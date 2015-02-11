%% Illustrating the ParTI method on gene expression in Cancer (RNAseq)
% This script an example of how to analyze gene expression data with the ParTI method.
% It can be used as a template to analyze other datasets.

% Load the data into Matlab from a comma separated value (CSV) file
% The file is a purely numerical matrix, with patients as rows and genes as
% columns
geneExpression = dlmread('Data/CancerRNAseq/CancerRNAseq_expression.csv', ',');
% The file is formated as samples (i.e. patients) x genes. 
% We load gene names.
geneNames = importdata('Data/CancerRNAseq/CancerRNAseq_geneNames.list');

%% We import the sample attributes, i.e. the clinical data on patients
% These come in two kinds: 
% - discrete attributes, i.e. categorical data (citizenship, gender, cancer progression grade, ...)
% - continuous attributes, i.e. numerical data (weight, age, tumor volume, ...)
% We start by loading a file with clinical attributes, both discrete and continuous.
[discrAttrNames, discrAttr] = ...
    read_enriched_csv('Data/CancerRNAseq/CancerRNAseq_discreteClinicalData.tsv', char(9));
%where discrAttr is a matrix of patients x attributes. The names of
%the attributes are stored in discrAttrNames.

% Continuous attributes are analysed using a different statistical procedure
% than discrete attributes. We therefore load a file with vectors 'discIdcs' and 
% 'contIdcs' containing the indices of discrete and continuous attributes 
% to be considered in the study:
load Data/CancerRNAseq/CancerRNAseq_featIdcs.mat;

% We extract continous attributes, discarding the sample ID
% (the order of the clinical records matches that of the expression data)
contAttrNames = discrAttrNames(:,contIdcs);
contAttr = discrAttr(:,contIdcs);
% Now we convert contAttr to a matrix of doubles
contAttr = str2double(contAttr);

% And select the remaining discrete attributes:
discrAttrNames = discrAttrNames(:,discIdcs);
discrAttr = discrAttr(:,discIdcs);

%% We expand the sample attributes by computing changes in GO category expression
% This section is optional. It makes it possible to determine broad gene 
% expression categories that are over-expressed in the vicinity of 
% archetypes. This is helpful to characterize the archetypes.
[GOExpression,GONames,~,GOcat2Genes] = MakeGOMatrix(geneExpression, geneNames, ...
                {'MSigDB/c2.cp.v4.0.symbols.gmt', 'MSigDB/c5.all.v4.0.symbols.gmt'}, ...
                100);
% GOExpression is a matrix of patients x GO categories, and
% GONames contains the name of the GO categories.
% GOcat2Genes is a boolean matrix of genes x GO categories which
% indicates, for each category, which genes were used to compute it.
% In the next line, we expand this matrix so that it has as many columns as
% the number of continuous features (clinical + GO). Because clinical
% features are typically not directly based on specific genes, we add
% zeroes in the corresponding columns:
GOcat2Genes=[zeros(size(GOcat2Genes,1),size(contAttr,2)),GOcat2Genes];
% and we expand the continuous clinical features with GO-based continuous
% features:
contAttrNames = [contAttrNames, GONames];
contAttr = [contAttr, GOExpression];

%% Finally, we substitute underscores '_' in variable names with spaces ' ' 
% to prevent the characters following underscores from appearing in indice
% position.
discrAttrNames = regexprep(discrAttrNames, '_', ' ');
contAttrNames = regexprep(contAttrNames, '_', ' ');

%% We are now ready to perform Pareto Task Inference.
% We use the Sisal algorithm (1), with up to 8 dimensions. We provide the
% discrete patient attributes, and ask ParTI to preliminary booleanize these
% attributes (0). We also pass continuous patient attributes. We pass a boolean 
% matrix specifiying which genes each continuous feature is baesd on (to be used
% in the leave-one-out procedure). We specify that the enrichment analysis 
% will be performed with a bin size of 5%. Finally, the output of the the 
% analysis will be stored in an Excel spreadsheet, under the name 
% 'CancerRNAseq_enrichmentAnalysis_*.csv'.
[arc, arcFinal, pc] = ParTI_lite(geneExpression, 1, 8, discrAttrNames, ...
    discrAttr, 0, contAttrNames, contAttr, [], 0.05, 'CancerRNAseq_enrichmentAnalysis');

%% Finally, we perform the compete analysis, including randomization
% controls and archetype error estimation.
[arc, arcOrig, pc, errs, pval] = ParTI(geneExpression, 1, 8, discrAttrNames, ...
    discrAttr, 0, contAttrNames, contAttr, GOcat2Genes, 0.05, 'CancerRNAseq_enrichmentAnalysis');
