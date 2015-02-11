%% Illustrating ParTI on gene expression in Mouse tissues (microarrays)
% This script an example of how to analyze gene expression data with ParTI.
% It can be used as a template to analyze other datasets.

% Load the data into Matlab from a comma separated value (CSV) file
% The file is a purely numerical matrix, with patients as rows and genes as
% columns
geneExpression = dlmread('Data/Mouse/Mouse_expression.tsv', char(9));
% The file is formated as samples (i.e. tissues) x genes. 
% We load gene names.
geneNames = importdata('Data/Mouse/Mouse_geneNames.list');

%% We import the sample attributes, i.e. the tissue types
[discrAttrNames, discrAttr] = ...
    read_enriched_csv('Data/Mouse/Mouse_tissueAnnotation.tsv', char(9));
% where discrAttr is a matrix of 89 tissue x 2 attributes. The names of
% the attributes are stored in discrAttrNames.
% The first colum of this table contains tissue names, which we won't need in
% the analysis, so we now discard them.
tissueNames = discrAttr(:,1);
discrAttr = discrAttr(:,2);
discrAttrNames = discrAttrNames(2);

% In this dataset, we did not consider continous attributes. If we had, we
% would load them now.
contAttrNames = [];
contAttr = [];
% [contAttrNames, contAttr] = ...
%    read_enriched_csv('Data/Mouse/Mouse_continuousTissueAnnotation.tsv', char(9));

%% We focus on samples with homogenous cell composition
isHomogenous = find(~strcmp(discrAttr(:,1), 'other (heterogenous)'));
geneExpression = geneExpression(isHomogenous,:);
discrAttr = discrAttr(isHomogenous,:);
tissueNames = tissueNames(isHomogenous,1);

%% We expand the sample attributes by computing changes in GO category expression
% This section is optional. It makes it possible to determine broad gene 
% expression categories that are over-expressed in the vicinity of 
% archetypes. This is helpful to characterize the archetypes.
[GOExpression,GONames,~,GOcat2Genes] = MakeGOMatrix(geneExpression, geneNames, ...
                {'MSigDB/c2.cp.v4.0.symbols.gmt', 'MSigDB/c5.all.v4.0.symbols.gmt'}, ...
                100);
% GOExpression is a matrix of samples x GO categories, and
% GONames contains the name of the GO categories. Finally, we expand the 
% continuous attributes
% GOcat2Genes is a boolean matrix of genes x GO categories which
% indicates, for each category, which genes were used to compute it.
GOcat2Genes=[zeros(size(GOcat2Genes,1),size(contAttr,2)),GOcat2Genes];
contAttrNames = [contAttrNames, GONames];
contAttr = [contAttr, GOExpression];

%% We are now ready to perform Pareto Task Inference.
% We use the MVSA algorithm (2), with up to 8 dimensions. We provide the
% discrete patient attributes, and ask ParTI to preliminary booleanize these
% attributes (0). We also pass continuous patient attributes. We also pass 
% continuous patient features. We pass a boolean matrix specifiying which 
% genes each continuous feature is baesd on (to be used in the 
% leave-one-out procedure). We specify that the enrichment analysis will 
% be performed with a bin size of 20%. 
% Finally, the output of the the analysis will be stored in an Excel
% spreadsheet, under the name 'Mouse_enrichmentAnalysis_*.csv'.
[arc, arcOrig, pc] = ParTI_lite(geneExpression, 2, 8, discrAttrNames, ...
     discrAttr, 0, contAttrNames, contAttr, [], 0.2, 'Mouse_enrichmentAnalysis');

%% Finally, we perform the compete analysis, including randomization
% controls and archetype error estimation.
[arc, arcOrig, pc, errs, pval] = ParTI(geneExpression, 2, 8, discrAttrNames, ...
     discrAttr, 0, contAttrNames, contAttr, GOcat2Genes, 0.2, 'Mouse_enrichmentAnalysis');
