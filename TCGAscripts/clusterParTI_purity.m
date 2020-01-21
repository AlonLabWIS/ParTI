%% ParTI pipeline for TCGA datasets
addpath(genpath('../ParTI/'))
origPath = pwd;
myQuantile = 0.0;%Metabric was .3
%nArchetypes = 4;

%global ForceNArchetypes; ForceNArchetypes = nArchetypes;

purity = 'low';

% Load the data into Matlab from a comma separated value (CSV) file
% The file is a purely numerical matrix, with patients as rows and genes as
% columns
geneExpression = dlmread(['expMatrix_', purity, 'Purity.csv'], ',');
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
if length(strfind(origPath, 'metabric')) == 1
    minExpr = quantile(mean(geneExpression,1), myQuantile);
    selGenes = find(mean(geneExpression,1) >= minExpr);
else
    genePctilesAvg = dlmread('../../genePctilesAvg.tab', ',');
    selGenes = find(genePctilesAvg >= myQuantile);   
end
geneExpression = geneExpression(:,selGenes);
geneNames = geneNames(selGenes,:);
cell2csv('geneNamesAfterExprFiltering.list', geneNames);

%% We import the sample attributes, i.e. the clinical data on patients
% These come in two kinds: 
% - discrete attributes, i.e. categorical data (citizenship, gender, cancer progression grade, ...)
% - continuous attributes, i.e. numerical data (weight, age, tumor volume, ...)
% We start by loading discrete attributes.
[discrAttrNames, discrAttr] = ...
    read_enriched_csv(['discreteClinicalData_reOrdered_', purity, 'Purity.tsv'], char(9));
%where discrAttr is a matrix of 2106 patients x 25 attributes. The names of
%the attributes are stored in discrAttrNames.

%Load continuous features
%[contAttrNames, contAttr] = ...
%   read_enriched_csv('continuousClinicalData_reOrdered.tsv', char(9));
if exist(['continuousClinicalData_reOrdered_justData_', purity, 'Purity.csv'], 'file') == 2
    contAttr = dlmread(['continuousClinicalData_reOrdered_justData_', purity, 'Purity.csv'], ',');
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
    if exist('GONames', 'var')
        GOExpression = GOExpression(noNormal,:);
    end
    discrAttr = discrAttr(noNormal,:);
    contAttr = contAttr(noNormal,:);
else
    noNormal = 1:size(discrAttr,1);
end

binSize=50/size(geneExpression,1); % 50 samples per bin by default
if binSize < .05
    binSize = .05; %5% bin size at least
end

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
    csvwrite(strcat(origPath, '/arcs_dims.csv'), arc)
    %save(strcat(origPath, '/arcsOrig_genes.tsv'), 'arcOrig', '-ascii')
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

%[arc, arcOrig, pc] = ParTI(geneExpression);

%[arc, arcOrig, pc, errs, pval] = ParTI(geneExpression, 1, 8, discrAttrNames, ...
%    discrAttr, 0, GONames, GOExpression, GOcat2Genes, 0.1, ...
%    strcat(origPath, '/clinicalMSigDBpValEnrichment'));
