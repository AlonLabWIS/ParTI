function a=GoCategoriesLeaveOneOut(GOidx,GOcat2Genes,DataPoints,NArchetypes,algNum,numIter,EnMatCont,binSize)

%% Inputs
% GOidx - the index of all significant GO categories
% GOcat2Genes - A matrix storing for each GO category (column) the list of
% genes that are related to that GO category 
% DataPoints - Data points is the values of different traits (e.g. expression
% level of genes) - each sample is a row, each trait (gene) is a column
% NArchetypes - Number of archetypes to look for in the new projected data
% set
% algNum - states which archetype algorythem to use: 
%    algNum=1 :> Sisal (default)
%    algNum=2 :> MVSA
%    algNum=3 :> MVES
%    algNum=4 :> SDVMM
%    algNum=5 :> PCHA
% numIter - The number of iterations to run the algorithm for minimal
% bounding simplex
% EnMatCont - is a real matrix with each Continuous feature as a column, each row is a
%  datapoint (order of datapoints and their features should match).
% binSize - the bin size
% 
% Output
% a boolean vector indicating whether enrichment is still significant after
% leave-one-out (1) or not (0)

%% Initialization of data points to the new dimensions
indGO=find(GOidx)';
a=zeros(1,length(GOidx));

indGOlen = length(indGO);
for idx=1:indGOlen
%for idx=1:5 %DEBUG
    i = indGO(idx);
    if mod(idx,round(indGOlen/20)) == 0
        fprintf('%.0f%% done\n', 100*idx/indGOlen);
    end

    % Getting the list of genes that are in a certain GO category
    indGo2Genes=GOcat2Genes(:,i);
    if any(indGo2Genes) %check if the continuous feature is a GO category, if it is not, keep it significant
        % dropping the dimensions that are part of the specific GO category
        DataPointsWO=DataPoints(:,~indGo2Genes);
        
        %% Do PCA on the data
        %     fprintf('Starting to perform PCA, for big data on slow computers this may take a while...\n');
        [~,scores1,~] = pca(DataPointsWO);
        
        DataPCA=scores1;
        
        %% find the archetypes of the reduced data minimal simplex
        ArchsMin=findMinSimplex(numIter,DataPCA,algNum,NArchetypes,1);
        ArchsMin=ArchsMin';
        %% Calculating the distance of each point from the archetypes
        % For each Archetype sort the points by their distance to the archetype
        orderedIndices = sortDataByDistance(DataPCA(:,1:size(ArchsMin,2)),ArchsMin);
%         disp('Finished sorting data points.');
        
        %% Calculate enrichment with the new ordering but with the real points
        GOcontinuousTable = ContinuousEnrichment(orderedIndices,EnMatCont(:,i),binSize);
%         disp('Finished computing continuous enrichments.');
        
        %% Compare Enrichment in both cases to see whether a GO category is still enriched
        % The sixth column of the table holds whether the GO category is enriched or not
        a(i)=GOcontinuousTable(6);
    else % if this is not a GO category, keep it significant
        a(i)=1; %
    end
end

 
