function []=calculateEnrichment_lite(DataPoints,Archetypes,DiscFeatName,EnMatDis,ContFeatName,EnMatCont,binSize,OutputFileName)

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
% The output is separate CSV files for discrete and continuous features 
% with enrichment statistics.

% Initialization of parameters
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
   DiscFeatName = cellfun(@(x)['Feature ',num2str(x)]...
       ,num2cell(1:size(EnMatDis,2)),'UniformOutput',false);
end
discreteTable = DiscreteEnrichment_lite(orderedIndices,EnMatDis,binSize,DiscFeatName,OutputFileName);


 disp('Finished computing discrete enrichments.');
    % Calculate enrichment (P/Q) + p-value
    
    % Calculate p-value for maximally enriched
    
    % Correct for multiple hypothesis testing with BH procedure

if size(ContFeatName) == 0
    ContFeatName = cellfun(@(x)['Feature ',num2str(x)]...
       ,num2cell(1:size(EnMatCont,2)),'UniformOutput',false);
end

continuousTable = ContinuousEnrichment(orderedIndices,EnMatCont,binSize);  
disp('Finished computing continuous enrichments.');
% Calculate enrichment and significance of continuous features (Median difference
% + p-value of Mann-Whitney + Benjamini-Hochberg)

    % Calculate enrichment (P/Q) + p-value
    
    % Calculate p-value for maximally enriched
    
    % Correct for multiple hypothesis testing with BH procedure


% GO analysis for gene expression

% Plot all features enrichment for all archetypes

% Create tables of enrichment in excel (sheet1-Discrete, sheet2-Continuous)
%Discrete
if sum(~isnan(discreteTable)) > 0 
    DiscreteTitles = {'archetype #', 'Feature Name',  'P value (Hypergeom.)',...
        'Significant after Benjamini-Hochberg correction?','Is first bin maximal?'};  
    ordDiscTable = sortrows(discreteTable,[1 -4 3]);
    ordDiscFeatureNames = DiscFeatName(ordDiscTable(:,2))';
    fullTable = [DiscreteTitles; num2cell(ordDiscTable)];
    fullTable(2:end,2) = ordDiscFeatureNames;
    cell2csv([OutputFileName '_discrete_lite_All.csv'], fullTable);  
%     xlswrite([OutputFileName '_lite_All.xls'], DiscreteTitles     ,1,'A1:E1');
%     xlswrite([OutputFileName '_lite_All.xls'], ordDiscTable       ,1,'A2');
%     xlswrite([OutputFileName '_lite_All.xls'], ordDiscFeatureNames,1,'B2');

    filteredDiscTable=ordDiscTable(ordDiscTable(:,4)>0 & ordDiscTable(:,5)== 1,:);
    FiltDiscFeatureNames = DiscFeatName(filteredDiscTable(:,2))';
    
    fullTable = [DiscreteTitles; num2cell(filteredDiscTable)];
    fullTable(2:end,2) = FiltDiscFeatureNames;
    cell2csv([OutputFileName '_discrete_lite_significant.csv'], fullTable);
%     xlswrite([OutputFileName '_lite_significant.xls'], DiscreteTitles          ,1,'A1:E1');
%     xlswrite([OutputFileName '_lite_significant.xls'], filteredDiscTable       ,1,'A2');
%     xlswrite([OutputFileName '_lite_significant.xls'], FiltDiscFeatureNames    ,1,'B2');
end

%Continuous
if sum(~isnan(continuousTable)) > 0

    ContinuousTitles = { 'archetype #','Feature Name', 'P value (Mann-Whitney)'...
        ,'Median Difference' ,'Mean Difference','Significant after Benjamini-Hochberg correction?','Is first bin maixmal?'}; 
    ordContTable = sortrows(continuousTable,[1 -6 -7 -4]);
    ordContFeatureNames = ContFeatName(ordContTable(:,2))';
    
    fullTable = [ContinuousTitles; num2cell(ordContTable)];
    fullTable(2:end,2) = ordContFeatureNames;
    cell2csv([OutputFileName '_continuous_lite_All.csv'], fullTable);
%     xlswrite([OutputFileName '_lite_All.xls'], ContinuousTitles     ,2,'A1:F1');
%     xlswrite([OutputFileName '_lite_All.xls'], ordContTable       ,2,'A2');
%     xlswrite([OutputFileName '_lite_All.xls'], ordContFeatureNames,2,'B2');

    filteredContTable=ordContTable(ordContTable(:,4)>0 &ordContTable(:,6)>0&ordContTable(:,7)>0,:);
    FiltContFeatureNames = ContFeatName(filteredContTable(:,2))';
    fullTable = [ContinuousTitles; num2cell(filteredContTable)];
    fullTable(2:end,2) = FiltContFeatureNames;
    cell2csv([OutputFileName '_continuous_lite_significant.csv'], fullTable);  
%     xlswrite([OutputFileName '_lite_significant.xls'], ContinuousTitles          ,2,'A1:F1');
%     xlswrite([OutputFileName '_lite_significant.xls'], filteredContTable         ,2,'A2');
%     xlswrite([OutputFileName '_lite_significant.xls'], FiltContFeatureNames      ,2,'B2');
end