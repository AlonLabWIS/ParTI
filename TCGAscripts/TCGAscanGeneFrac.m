%% ParTI pipeline for TCGA datasetson

% Default setting

%myQuantile = .5;
%nArchetypes = 4;

%%
addpath ../ParTI/
origPath = pwd;
% disp(origPath);

% Load the data into Matlab from a comma separated value (CSV) file
% The file is a purely numerical matrix, with patients as rows and genes as
% columns
geneExpression = dlmread('expMatrix.csv', ',');
% The file is formated as samples (i.e. patients) x genes. 
% We load gene names.
% geneNames = importdata('geneListExp.list');

[discrAttrNames, discrAttr] = ...
    read_enriched_csv('discreteClinicalData_reOrdered.tsv', char(9));

%% Select genes
% hist(reshape(geneExpression, 1, numel(geneExpression)),30);
% [f,x] = ecdf(reshape(geneExpression, 1, numel(geneExpression))); plot(x,f);
% clear f x;
% minExpr = 2;
% minExpr = 8;
minExpr = quantile(mean(geneExpression,1), myQuantile);
selGenes = find(mean(geneExpression,1) > minExpr);
geneExpression = geneExpression(:,selGenes);
% geneNames = geneNames(selGenes,:);

% [arc, arcOrig, pc] = ParTI_lite(geneExpression);

%% Fill in number of desired archetypesand quit after computing p-value
global ForceNArchetypes; ForceNArchetypes = nArchetypes;
global abortAfterPval; abortAfterPval = 1;
% global lowIterations; lowIterations = 1;

%% Remove normal tissues
featIdx = find(strcmp(discrAttrNames, 'sample_type'));
noNormal = find(strcmp(discrAttr(:,featIdx), 'Solid Tissue Normal') == 0);
geneExpression = geneExpression(noNormal,:);

%% Finally, we perform the compete analysis, including randomization
% controls and archetype error estimation.
% close all
cd ../ParTI
[arc, arcOrig, pc, errs, pval] = ParTI(geneExpression);
