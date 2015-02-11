%% Illustrating ParTI on a synthetic dataset
% This script illustrates ParTI on a synthetic, 2D triangular data 
% set of 100 points. Each point has features a single discrete attribute.
% This script gives a quick feeling for the ParTI method.
%
% Load the data into Matlab from a comma separated value (CSV) file
% The file is a purely numerical matrix, with samples as rows and traits as
% columns
data = dlmread('Data/Synthetic/Synthetic_traits.csv', ',');
% The file is formated as samples (i.e. points) x traits. 

%% We import the sample attributes.
% In this case, it is just a discrete variable that can values of A, B or C
% for each point.
[discrAttrNames, discrAttr] = ...
    read_enriched_csv('Data/Synthetic/Synthetic_attributes.csv', ',');
% where discrAttr is a matrix of 100 points x 1 attributes. The names of
% the attributes are stored in discrAttrNames.

% In this dataset, we did not consider continous attributes. If we had, we
% would load them now.
contAttrNames = [];
contAttr = [];
% [contAttrNames, contAttr] = ...
%    read_enriched_csv('Data/Synthetic/Synthetic_continuousTissueAnnotation.tsv', char(9));

%% We are now ready to perform Pareto Task Inference (light version).
% We use the SISAL algorithm (1), with up to 2 dimensions because the
% dataset has only two dimensions. We provide the discrete attributes,
% and ask ParTI to preliminary booleanize these attributes (0). We
% specify that the enrichment analysis will be performed with a bin size
% of 20%. Finally, the output of the the analysis will be stored in an
% Excel spreadsheet, under the name 'Synthetic_enrichmentAnalysis_*.csv'.
%
% NOTE: If you use this script as a template to analyze more complex
% datasets, remember to increase the dimension parameter (2 in this
% example) so that you can observe the point at which the explained
% variance curve saturates.
[arc, arcOrig, pc] = ParTI_lite(data, 1, 2, discrAttrNames, ...
    discrAttr, 0, contAttrNames, contAttr, [], 0.2, 'Synthetic_enrichmentAnalysis');
