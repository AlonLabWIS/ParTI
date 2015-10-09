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

cd ../ParTI

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

%% Finally, we perform the compete analysis, including randomization
% controls and archetype error estimation.
% close all
[arc, arcOrig, pc, errs, pval] = ParTI(geneExpression);
