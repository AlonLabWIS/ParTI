function []=calculateEnrichment(DataPoints,Archetypes,DiscFeatName,EnMatDis,ContFeatName,EnMatCont,binSize,OutputFileName,numIter,algNum,GOcat2Genes,DataPointsOrg,evalPmax)

%Inputs
% 1. Data points is the values of different traits (e.g. expression
% level of genes) - each sample is a row, each trait (gene) is a column 
% 2. Archetypes is a matrix with different archetypes in rows
% 3. DiscFeatName - List of the Discrete Feature labels (row)
% 4. EnMatDis is a boolean matrix with each Discrete feature as a column, each row is a
%  datapoint (order of datapoints and their features should match, order of columns and feature labels should match).
% 5. ContFeatName - List of the Continuous Feature labels
% 6. EnMatCont is a real matrix with each Continuous feature as a column, each row is a
%  datapoint (order of datapoints and their features should match).
% 7. binSize is the value of the percent each bin should have in the
%    enrichment calculation
% 8. OutputFileName is the name of the excel file that saves all enrichment
% 9. numIter is the number of iterations for the algorithm to find a minimal
% bounding simplex
% 10. algNum chooses the algorithm to run to find the bounding simplex (see
% findMinSimplex)
% 11. GOcat2Genes is a boolean matrix (genes x continuous features) that 
% indicates, for each continuous feature, which genes this feature is based 
% on. This matrix, used for the leave-on-out analysis, can be readily 
% computed by the MakeGOMatrix function.
% 12. DataPointsOrg is the original data passed to ParTI prior to PCA, which 
% is used in the leave-one-out analysis.
% 13. evalPmax determines whether to compute the probability that enrichment is maximal in the first bin. Setting this parameter to 0 skips this computation to save time.
% The output is separate CSV files for discrete and continuous features 
% with enrichment statistics.

% Initialization of parameters
minPmax = .5; %minimal probability to be maximally enriched in the bin closest to an archetype that we require to report the category as significantly enriched at the archetype

if binSize <= 0 || binSize > .5
    [~, numDataPoints] = size(DataPoints);
    defaultBinsN = round(numDataPoints / 10);
    if defaultBinsN < 2
        defaultBinsN = 2;
    end
    binSize = 1 / defaultBinsN;
    fprintf('binSize should be stricly larger than 0 and smaller or equal to 0.5!\nFor now, we will set it to %f\n', binSize);
end

% For each Archetype sort the points by their distance to the archetype
orderedIndices = sortDataByDistance(DataPoints,Archetypes);
disp('Finished sorting data points.');
    %orderedIndices list of indexes + distances
   
% Calculate enrichment and significance of discrete features (P/Q + p-value of HG test+
% + p-value of maximally enriched + Benjamini-Hochberg)
if size(DiscFeatName) == 0
   DiscFeatName = cellfun(@(x)['Feature ',num2str(x)],...
       num2cell(1:size(EnMatDis,2)),'UniformOutput',false);
end
discreteTable = DiscreteEnrichment(orderedIndices,EnMatDis,binSize,DiscFeatName,OutputFileName,evalPmax);

disp('Finished computing discrete enrichments.');
    % Calculate enrichment (P/Q) + p-value
    
    % Calculate p-value for maximally enriched
    
    % Correct for multiple hypothesis testing with BH procedure

if size(ContFeatName) == 0
    ContFeatName = cellfun(@(x)['Feature ',num2str(x)],...
       num2cell(1:size(EnMatCont,2)),'UniformOutput',false);
end

continuousTable = ContinuousEnrichment(orderedIndices,EnMatCont,binSize);  
disp('Finished computing continuous enrichments.');

%calculate the significance after a leave-1-out procedure
if length(continuousTable)>1 && evalPmax && ~isempty(ContFeatName)
    disp('Now applying leave-one-out to verify robustness enrichment of GO categories');
    NArchetypes=size(Archetypes,1);
    
    sigGOs=continuousTable(:,6)&continuousTable(:,7);
    sigGOind=any(reshape(sigGOs,NArchetypes,length(sigGOs)/NArchetypes))';

    isSignificantL1o=GoCategoriesLeaveOneOut(sigGOind,GOcat2Genes,DataPointsOrg,NArchetypes,algNum,numIter,EnMatCont,binSize);
    isSigL1oSpread=repmat(isSignificantL1o,NArchetypes,1);
    continuousTable(:,6)=continuousTable(:,6)&isSigL1oSpread(:);
end

% Calculate enrichment and significance of continuous features (Median difference
% + p-value of Mann-Whitney + Benjamini-Hochberg)

    % Calculate enrichment (P/Q) + p-value
    
    % Calculate p-value for maximally enriched
    
    % Correct for multiple hypothesis testing with BH procedure


% GO analysis for gene expression

% Plot all features enrichment for all archetypes

% Create tables of enrichment in excel (sheet1-Discrete, sheet2-Continuous)
%Discrete
if size(discreteTable,2) >= 3 && sum(~isnan(discreteTable(:,3))) > 0 
	if evalPmax
		DiscreteTitles = {'archetype #', 'Feature Name',  'P value (Hypergeom.)',...
	        'Significant after Benjamini-Hochberg correction?','P maximal'};
	else
		DiscreteTitles = {'archetype #', 'Feature Name',  'P value (Hypergeom.)',...
			'Significant after Benjamini-Hochberg correction?','Is first bin maximal?'};  
		minPmax = 1; %we don't have Pmax, only a boolean saying if the frequency of the feature is highest in the bin closest to the archetype
	end
    ordDiscTable = sortrows(discreteTable,[1 -4 3]);
    ordDiscFeatureNames = DiscFeatName(ordDiscTable(:,2))';
    fullTable = [DiscreteTitles; num2cell(ordDiscTable)];
    fullTable(2:end,2) = ordDiscFeatureNames;
    cell2csv(sprintf('%s_discrete_All.csv', OutputFileName), fullTable);
%     xlswrite([OutputFileName '_All.xls'], DiscreteTitles     ,1,'A1:E1');
%     xlswrite([OutputFileName '_All.xls'], ordDiscTable       ,1,'A2');
%     xlswrite([OutputFileName '_All.xls'], ordDiscFeatureNames,1,'B2');

    filteredDiscTable=ordDiscTable(ordDiscTable(:,4)>0 & ordDiscTable(:,5)>=minPmax,:);
    FiltDiscFeatureNames = DiscFeatName(filteredDiscTable(:,2))';
    
    fullTable = [DiscreteTitles; num2cell(filteredDiscTable)];
    fullTable(2:end,2) = FiltDiscFeatureNames;
    cell2csv(sprintf('%s_discrete_significant.csv', OutputFileName), fullTable);
    fprintf('*** Wrote summary of enriched discrete features to %s_discrete_significant.csv\n', OutputFileName);
%     xlswrite([OutputFileName '_significant.xls'], DiscreteTitles          ,1,'A1:E1');
%     xlswrite([OutputFileName '_significant.xls'], filteredDiscTable       ,1,'A2');
%     xlswrite([OutputFileName '_significant.xls'], FiltDiscFeatureNames    ,1,'B2');
end

%Continuous
if size(continuousTable,2) >= 3 && sum(~isnan(continuousTable(:,3))) > 0

    ContinuousTitles = { 'archetype #','Feature Name', 'P value (Mann-Whitney)'...
        ,'Median Difference' ,'Mean Difference','Significant after Benjamini-Hochberg correction?','Is first bin maximal?'}; 
    ordContTable = sortrows(continuousTable,[1 -6 -7 -4]);
    ordContFeatureNames = ContFeatName(ordContTable(:,2))';
    
    fullTable = [ContinuousTitles; num2cell(ordContTable)];
    fullTable(2:end,2) = ordContFeatureNames;
    cell2csv(sprintf('%s_continuous_All.csv', OutputFileName), fullTable);
%     xlswrite([OutputFileName '_All.xls'], ContinuousTitles     ,2,'A1:F1');
%     xlswrite([OutputFileName '_All.xls'], ordContTable       ,2,'A2');
%     xlswrite([OutputFileName '_All.xls'], ordContFeatureNames,2,'B2');

    filteredContTable=ordContTable(ordContTable(:,4)>0 &ordContTable(:,6)>0&ordContTable(:,7)>0,:);
    FiltContFeatureNames = ContFeatName(filteredContTable(:,2))';
    fullTable = [ContinuousTitles; num2cell(filteredContTable)];
    fullTable(2:end,2) = FiltContFeatureNames;
    cell2csv(sprintf('%s_continuous_significant.csv', OutputFileName), fullTable);
    fprintf('*** Wrote summary of enriched continuous features to %s_continuous_significant.csv\n', OutputFileName);
%     xlswrite([OutputFileName '_significant.xls'], ContinuousTitles          ,2,'A1:F1');
%     xlswrite([OutputFileName '_significant.xls'], filteredContTable         ,2,'A2');
%     xlswrite([OutputFileName '_significant.xls'], FiltContFeatureNames      ,2,'B2');
end
