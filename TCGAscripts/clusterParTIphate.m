%% ParTI pipeline for TCGA datasets
addpath ../ParTI/
origPath = pwd;
myQuantile = 0.3;%Metabric
%myQuantile = 0.4;%BLCA
%myQuantile = 0;
nArchetypes = 4;

global ForceNArchetypes; ForceNArchetypes = nArchetypes;

% Load the data into Matlab from a comma separated value (CSV) file
% The file is a purely numerical matrix, with patients as rows and genes as
% columns
geneExpression = dlmread('expMatrix.csv', ',');
% The file is formated as samples (i.e. patients) x genes. 
% We load gene names.
geneNames = importdata('geneListExp.list');

%% We expand the sample attributes by computing changes in GO category expression
% This section is optional. It makes it possible to determine broad gene 
% expression categories that are over-expressed in the vicinity of 
% archetypes. This is helpful to characterize the archetypes.
[GOExpression,GONames,~,GOcat2Genes] = MakeGOMatrix(geneExpression, geneNames, ...
                {'../ParTI/MSigDB/c2.cp.v4.0.symbols.gmt', '../ParTI/MSigDB/c5.all.v4.0.symbols.gmt'}, ...
                10);
GONames = regexprep(GONames, '_', ' ');
            
% GOExpression is a matrix of 2106 patients x 162 GO categories, and
% GONames contains the name of the GO categories.
% GOcat2Genes is a boolean matrix of genes x GO categories which
% indicates, for each category, which genes were used to compute it.
% In the next line, we expand this matrix so that it has as many columns as
% the number of continuous features (clinical + GO). Because clinical
% features are typically not directly based on specific genes, we add
% zeroes in the corresponding columns:
%GOcat2Genes=[zeros(size(GOcat2Genes,1),size(contAttr,2)),GOcat2Genes];

% We won't expand the continuous clinical features with GO-based continuous
% features but instead examine them in a separate analysis
%contAttrNames = [contAttrNames, GONames];
%contAttr = [contAttr, GOExpression];


%% Select genes
minExpr = quantile(mean(geneExpression,1), myQuantile);
selGenes = find(mean(geneExpression,1) >= minExpr);
geneExpression = geneExpression(:,selGenes);
geneNames = geneNames(selGenes,:);
cell2csv('geneNamesAfterExprFiltering.list', geneNames);

%% We import the sample attributes, i.e. the clinical data on patients
% These come in two kinds: 
% - discrete attributes, i.e. categorical data (citizenship, gender, cancer progression grade, ...)
% - continuous attributes, i.e. numerical data (weight, age, tumor volume, ...)
% We start by loading discrete attributes.
[discrAttrNames, discrAttr] = ...
    read_enriched_csv('discreteClinicalData_reOrdered.tsv', char(9));
%where discrAttr is a matrix of 2106 patients x 25 attributes. The names of
%the attributes are stored in discrAttrNames.

%Load continuous features
%[contAttrNames, contAttr] = ...
%   read_enriched_csv('continuousClinicalData_reOrdered.tsv', char(9));
if exist('continuousClinicalData_reOrdered_justData.csv', 'file') == 2
    contAttr = dlmread('continuousClinicalData_reOrdered_justData.csv', ',');
    contAttrNames = importdata('continuousClinicalData_reOrdered_featNames.list');
    contAttrNames = regexprep(contAttrNames, '_', ' ');
else
    contAttr = [];
    contAttrNames = [];
end

%% Finally, we substitute underscores '_' in variable names with spaces ' ' 
% to prevent the characters following underscores from appearing in indice
% position.
discrAttrNames = regexprep(discrAttrNames, '_', ' ');
%if exist('GOnames', 'var')
%    GONames = regexprep(GONames, '_', ' ');
%end

%% Remove normal tissues
featIdx = find(strcmp(discrAttrNames, 'sample type'));
if size(featIdx, 2) > 0
    noNormal = find(strcmp(discrAttr(:,featIdx), 'Primary Tumor') == 1);
    %noNormal = find(strcmp(discrAttr(:,featIdx), 'Solid Tissue Normal') == 1);
    geneExpression = geneExpression(noNormal,:);
    if exist('GOnames', 'var')
        GOExpression = GOExpression(noNormal,:);
    end
    discrAttr = discrAttr(noNormal,:);
    contAttr = contAttr(noNormal,:);
else
    noNormal = 1:size(discrAttr,1);
end

binSize=50/size(geneExpression,1); % 50 samples per bin by default
if size(geneExpression,1) * binSize > 100
    binSize = 100 / size(geneExpression,1);
end

%% Phate analysis
addpath ../phate
npca = 100;
pc = svdpca(geneExpression, npca, 'random');

% plot PCA
figure;
scatter3(pc(:,1), pc(:,2), pc(:,3), 5, 'b', 'filled');
%set(gca,'xticklabel',[]);
%set(gca,'yticklabel',[]);
%set(gca,'zticklabel',[]);
axis tight
title 'PCA before MAGIC'
xlabel 'PC1'
ylabel 'PC2'
zlabel 'PC3'

% Construct diffusion operator
k = 5; % nearest neighbors for adaptive bandwidth
a = 10; % decay of exponential kernel, a -> inf is uniform (unweighted) kernel
distfun = 'euclidean';
disp 'computing distances'
PDX = squareform(pdist(pc, distfun));
knnDST = sort(PDX,2);
disp 'computing kernel and operator'
epsilon = knnDST(:,k+1); % bandwidth(x) = distance to k-th neighbor of x
PDX = bsxfun(@rdivide,PDX,epsilon); % autotuning d(x,:) using epsilon(x)
GsKer = exp(-PDX.^a); % not really Gaussian kernel
GsKer = GsKer + GsKer'; % symmetrization
DiffDeg = diag(sum(GsKer,2)); % degrees
DiffOp = DiffDeg^(-1)*GsKer; % row stochastic
DiffOp_sp = sparse(DiffOp); % sparse version of operator

% find optimal t using VNE
DiffAff = DiffDeg^(-1/2)*GsKer*DiffDeg^(-1/2); % symmetric conjugate affinities
DiffAff = (DiffAff + DiffAff')/2; % clean up numerical inaccuracies to maintain symmetry
% Find the eigenvalues
disp 'Finding the eigenvalues'
[~,S] = svd(DiffAff); 
S = diag(S);
t_vec = 1:100;
disp 'Computing VNE'
H = nan(size(t_vec)); 
for I=1:length(t_vec)
    t = t_vec(I); 
    S_t=S.^t;
    P = S_t/sum(S_t);
    P=P(P>0);
    H(I) = -sum(P .* log(P)); % Entropy of eigenvalues
end
% Plot the entropy; choose a t in the flatter range after the 'knee' for generally best results
disp 'VNE plotted'
figure;
plot(t_vec,H,'*-')
xlabel('t')
ylabel('VNE')

% MAGIC
t = 10;
disp 'powering operator'
DiffOp_t = DiffOp^t;
disp 'imputing'
data_imputed = DiffOp_t * geneExpression;

% PCA after MAGC
npca = 100;
pc_magic = svdpca(data_imputed, npca, 'random');

% plot PCA MAGIC
figure;
scatter3(pc_magic(:,1), pc_magic(:,2), pc_magic(:,3), 5, 'b', 'filled');
%set(gca,'xticklabel',[]);
%set(gca,'yticklabel',[]);
%set(gca,'zticklabel',[]);
axis tight
title 'PCA after MAGIC'
xlabel 'PC1'
ylabel 'PC2'
zlabel 'PC3'

% PHATE
% t = 6;
distfun_mds = 'euclidean';
DiffOp_t = DiffOp^t;
disp 'potential recovery'
DiffOp_t(DiffOp_t<=eps)=eps;
DiffPot = -log(DiffOp_t);
disp(['MDS distfun: ' distfun_mds])
if strcmp(distfun_mds, 'euclidean')
    DiffPot_pca = svdpca(DiffPot, npca, 'random'); % to make pdist faster
end

% CMDS PHATE
D_DiffPot = squareform(pdist(DiffPot_pca, distfun_mds));
ndim = 3;
Y_cmds = randmds(D_DiffPot, ndim);

% plot Classic MDS PHATE
figure;
scatter3(Y_cmds(:,1), Y_cmds(:,2), Y_cmds(:,3), 5, 'b', 'filled');
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
set(gca,'zticklabel',[]);
axis tight
title 'CMDS PHATE'
xlabel 'MDS1'
ylabel 'MDS2'
zlabel 'MDS3'

% Metric MDS PHATE
ndim = 3;
opt = statset('display', 'iter');
Y_start = randmds(D_DiffPot, ndim);
Y_mmds = mdscale(D_DiffPot, ndim, 'options', opt, 'start', Y_start, 'Criterion', 'metricstress');

% plot MMDS PHATE
figure;
scatter3(Y_mmds(:,1), Y_mmds(:,2), Y_mmds(:,3), 5, 'b', 'filled');
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
set(gca,'zticklabel',[]);
axis tight
title 'MMDS PHATE'
xlabel 'MDS1'
ylabel 'MDS2'
zlabel 'MDS3'

%% replace gene expression
% geneExpression = data_imputed;
geneExpression = Y_mmds;

%% We are now ready to perform Pareto Task Inference.
% We use the Sisal algorithm (1), with up to 8 dimensions. We provide the
% discrete patient attributes, and ask ParTI to preliminary booleanize these
% attributes (0). We also pass continuous patient features. We pass a boolean 
% matrix specifiying which genes each continuous feature is baesd on (to be used
% in the leave-one-out procedure). 
% We specify that the enrichment analysis will be performed with a bin size 
% of 5%. Finally, the output of the the analysis will be stored in an
% Comma-Separated-Value text file, under the name 'Cancer_enrichmentAnalysis_*.csv'.

% noCtrls = find(~isnan(contAttr(:,5)));
% geneExpression = geneExpression(noCtrls,:);
% GOExpression = GOExpression(noCtrls,:);
% discrAttr = discrAttr(noCtrls,:);
% contAttr = contAttr(noCtrls,:);

cd ../ParTI
if exist(strcat(origPath, '/arcs_dims.tsv'), 'file') == 2
    fprintf('Reloading previously computed archetypes\n');
    load(strcat(origPath, '/arcs_dims.tsv'))
    arc = arcs_dims;
    load(strcat(origPath, '/arcsOrig_genes.tsv'))
    arcOrig = arcsOrig_genes;
else 
    %[arc, arcOrig, pc, coefs1] = ParTI_lite(geneExpression);
    [arc, arcOrig, pc, coefs1] = ParTI_lite(geneExpression, 1, ...
        10, [], [], 0, [], [], [], binSize, ...
        strcat(origPath, '/ParTI'));

    %save(strcat(origPath, '/arcs_dims.tsv'), 'arc', '-ascii')
    %save(strcat(origPath, '/arcsOrig_genes.tsv'), 'arcOrig', '-ascii')
    %save(strcat(origPath, '/pcsOrig_samplesXdims.tsv'), 'pc', '-ascii')
    %save(strcat(origPath, '/projOrig_varsXdims.tsv'), 'coefs1', '-ascii')
    
    save(strcat(origPath, '/arcs_dims.tsv'), 'arc', '-ascii')
    csvwrite(strcat(origPath, '/arcs_dims.csv'), arc)
    save(strcat(origPath, '/arcsOrig_genes.tsv'), 'arcOrig', '-ascii')
    csvwrite(strcat(origPath, '/arcsOrig_genes.csv'), arcOrig)
    %save(strcat(origPath, '/pcsOrig_samplesXdims.tsv'), 'pc', '-ascii')
    csvwrite(strcat(origPath, '/pcsOrig_samplesXdims.csv'), pc)
    %save(strcat(origPath, '/projOrig_varsXdims.tsv'), 'coefs1', '-ascii')
    csvwrite(strcat(origPath, '/projOrig_varsXdims.csv'), coefs1)
    %csvwrite(strcat(origPath, '/arcs_errs.csv'), errs)
    %csvwrite(strcat(origPath, '/arcs_pval.csv'), pval)
end

clear pc coefs1;

%% Genes
close all
ParTI_lite(geneExpression, 1, size(arcOrig,1), [], ...
    [], 0, geneNames, geneExpression, [], binSize, ...
    strcat(origPath, '/geneEnrichment'), arcOrig);

%% Clinical features
close all
ParTI_lite(geneExpression, 1, size(arcOrig,1), discrAttrNames, ...
    discrAttr, 0, contAttrNames, contAttr, [], binSize, ...
    strcat(origPath, '/clinicalEnrichment'), arcOrig);

%% MSigDB gene groups
clear discrAttr contAttr;
close all
ParTI_lite(geneExpression, 1, size(arcOrig,1), [], [], ...
    0, GONames, GOExpression, [], binSize, ...
    strcat(origPath, '/MSigDBenrichment'), arcOrig);

%% Mutations
clear GOExpression GONames GOcat2Genes;

mut = dlmread(strcat(origPath, '/mutMatrix_reOrdered_booleanized_justData.csv'), ',');
mutNames = importdata(strcat(origPath, '/mutMatrix_reOrdered_booleanized_geneNames.list'));
mut = mut(noNormal,:);

% Project mutations to gene groups
% posMut = find(cellfun('length', regexp(mutNames, '=1$')') > 0);
% [GOmut,GOmutNames,~,GOcat2mutGenes] = MakeGOMatrix(mut(:,posMut), strrep(mutNames(posMut), '=1', ''), ...
%                 {'MSigDB/c2.cp.v4.0.symbols.gmt', 'MSigDB/c5.all.v4.0.symbols.gmt'}, ...
%                 5);

close all
ParTI_lite(geneExpression, 1, size(arcOrig,1), mutNames, ...2
    mut, -1, [], [], [], binSize, ...
    strcat(origPath, '/mutEnrichment'), arcOrig);

% close all
% ParTI_lite(geneExpression, 1, ForceNArchetypes, [], ...
%     [], [], GOmutNames, GOmut, [], 0.1, ...
%     strcat(origPath, '/MSigDBmutEnrichment'), arcOrig);
% clear GOmut GOmutNames GOcat2mutGenes;

%% Copy number alterations
clear mut mutNames pc;

cop = dlmread(strcat(origPath, '/copMatrix_reOrdered_booleanized_justData.csv'), ',');
copNames = importdata(strcat(origPath, '/copMatrix_reOrdered_booleanized_geneNames.list'));
cop = cop(noNormal,:);

close all
ParTI_lite(geneExpression, 1, size(arcOrig,1), copNames, cop, ...
    0, [], [], [], binSize, ...
    strcat(origPath, '/cnvEnrichment'), arcOrig);


%% Finally, we perform the compete analysis, including randomization
% controls and archetype error estimation.
% close all
[arc, arcOrig, pc] = ParTI(geneExpression);

%[arc, arcOrig, pc, errs, pval] = ParTI(geneExpression, 1, 8, discrAttrNames, ...
%    discrAttr, 0, GONames, GOExpression, GOcat2Genes, 0.1, ...
%    strcat(origPath, '/clinicalMSigDBpValEnrichment'));
