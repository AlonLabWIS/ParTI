%% ParTI pipeline for TCGA datasetson
addpath ../ParTI/
origPath = pwd;
% myQuantile = 0.0;
% nArchetypes = 4;

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
                20);

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

binSize=.1; % 10% by default
if length(selGenes) * binSize > 100
    binSize = 100 / length(selGenes);
end

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
GONames = regexprep(GONames, '_', ' ');

%% Remove normal tissues
% featIdx = find(strcmp(discrAttrNames, 'sample_type'));
% noNormal = find(strcmp(discrAttr(:,featIdx), 'Solid Tissue Normal') == 0);
% geneExpression = geneExpression(noNormal,:);

%% We are now ready to perform Pareto Task Inference.
% We use the Sisal algorithm (1), with up to 8 dimensions. We provide the
% discrete patient attributes, and ask ParTI to preliminary booleanize these
% attributes (0). We also pass continuous patient features. We pass a boolean 
% matrix specifiying which genes each continuous feature is baesd on (to be used
% in the leave-one-out procedure). 
% We specify that the enrichment analysis will be performed with a bin size 
% of 5%. Finally, the output of the the analysis will be stored in an
% Comma-Separated-Value text file, under the name 'Cancer_enrichmentAnalysis_*.csv'.

cd ../ParTI
if exist(strcat(origPath, '/arcs_dims.tsv'), 'file') == 2
    fprintf('Reloading previously computed archetypes\n');
    load(strcat(origPath, '/arcs_dims.tsv'))
    arc = arcs_dims;
    load(strcat(origPath, '/arcsOrig_genes.tsv'))
    arcOrig = arcsOrig_genes;
else
    [arc, arcOrig, ~] = ParTI_lite(geneExpression);
    save(strcat(origPath, '/arcs_dims.tsv'), 'arc', '-ascii')
    save(strcat(origPath, '/arcsOrig_genes.tsv'), 'arcOrig', '-ascii')
end

%% Clinical features
close all
ParTI_lite(geneExpression, 1, ForceNArchetypes, discrAttrNames, ...
    discrAttr, 0, contAttrNames, contAttr, [], 0.1, ...
    strcat(origPath, '/clinicalEnrichment'), arcOrig);

%% MSigDB gene groups
clear discrAttr contAttr;
close all
ParTI_lite(geneExpression, 1, ForceNArchetypes, [], [], ...
    0, GONames, GOExpression, [], 0.1, ...
    strcat(origPath, '/MSigDBenrichment'), arcOrig);

%% Mutations
clear GOExpression GONames GOcat2Genes;

mut = dlmread(strcat(origPath, '/mutMatrix_reOrdered_booleanized_justData.csv'), ',');
mutNames = importdata(strcat(origPath, '/mutMatrix_reOrdered_booleanized_geneNames.list'));

% Project mutations to gene groups
posMut = find(cellfun('length', regexp(mutNames, '=1$')') > 0);
[GOmut,GOmutNames,~,GOcat2mutGenes] = MakeGOMatrix(mut(:,posMut), strrep(mutNames(posMut), '=1', ''), ...
                {'MSigDB/c2.cp.v4.0.symbols.gmt', 'MSigDB/c5.all.v4.0.symbols.gmt'}, ...
                5);

close all
ParTI_lite(geneExpression, 1, ForceNArchetypes, mutNames, ...
    mut, -1, [], [], [], 0.1, ...
    strcat(origPath, '/mutEnrichment'), arcOrig);

close all
ParTI_lite(geneExpression, 1, ForceNArchetypes, [], ...
    [], [], GOmutNames, GOmut, [], 0.1, ...
    strcat(origPath, '/MSigDBmutEnrichment'), arcOrig);
clear GOmut GOmutNames GOcat2mutGenes;

%% Copy number alterations
clear mut mutNames pc;

cop = dlmread(strcat(origPath, '/copMatrix_reOrdered_booleanized_justData.csv'), ',');
copNames = importdata(strcat(origPath, '/copMatrix_reOrdered_booleanized_geneNames.list'));

close all
ParTI_lite(geneExpression, 1, ForceNArchetypes, copNames, cop, ...
    0, [], [], [], 0.1, ...
    strcat(origPath, '/cnvEnrichment'), arcOrig);


%% Finally, we perform the compete analysis, including randomization
% controls and archetype error estimation.
% close all
% [arc, arcOrig, pc] = ParTI(geneExpression);

%[arc, arcOrig, pc, errs, pval] = ParTI(geneExpression, 1, 8, discrAttrNames, ...
%    discrAttr, 0, GONames, GOExpression, GOcat2Genes, 0.1, ...
%    strcat(origPath, '/clinicalMSigDBpValEnrichment'));
