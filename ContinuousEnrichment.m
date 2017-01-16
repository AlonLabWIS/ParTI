function [table, pval, medianDifference, isSignificantAfterFDR] = ContinuousEnrichment(DataPointsInd,EnMatCont,binSize)

%Inputs
% 1. data point indexes sorted according to distances from archetypes (the 
%    first output of sortDataByDistane(DataPoints,Archetypes).
% 2. EnMatCont is a real matrix with each feature as a column, each row is a
%  datapoint (order of datapoints and their features should match, order of columns and feature labels should match).
% 3. binSize is the value of the percent each bin should have in the
%    enrichment calculation

% Initialazing
% significance threshold after FDR (Benjamini Hochberg). 
ThreshHoldBH = 0.1;

[Numarchs, numDataPoints] = size(DataPointsInd);
[numDataPoints2, numFeatures] = size(EnMatCont);
if(numDataPoints2 ~= numDataPoints) 
    table = NaN;
    return; 
end
numOfBins = round(1 / binSize);
%%
% numDataPoints = 105; numOfBins = 5; %just for testing
breakPoints = floor(linspace(0.5, numDataPoints + 0.5, numOfBins+1));
% diff(breakPoints) %just for testing 
numPointInBin = diff(breakPoints);
breakPoints = breakPoints(2:end);
%%
pval    = ones(Numarchs,numFeatures);
medianDifference  = zeros(Numarchs,numFeatures);
meanDifference  = zeros(Numarchs,numFeatures);
PoverQ  = zeros(Numarchs,numFeatures);

%Turn strings into doubles
if iscell(EnMatCont) & ischar(EnMatCont{1,1})
    EnMatCont = cell2mat(cellfun(@(y)str2num(y),EnMatCont,'UniformOutput',0));
end

%On Each Archetype
for arch = 1:Numarchs
   %Divide to bins
   tempEnrich =  EnMatCont(DataPointsInd(arch,:),:);    
   Binned = mat2cell(tempEnrich,numPointInBin, numFeatures);
   
   %EnPQ 
   P = cellfun(@(x)nanmedian(x),Binned,'UniformOutput',0);
   [~, Q] = nanmax(cell2mat(P));
   PoverQ(arch,:) = (Q == 1) ;
  
   tempEnrich =  EnMatCont(DataPointsInd(arch,:),:);
   %Divide to bins
   bin =  tempEnrich(1:(numPointInBin(1)),:);
   rest = tempEnrich((numPointInBin(1)+1):numDataPoints,:);

   %EnGen - calculate p-val by mann-whitney test for each feature 
   % Pval of first bin Vs. all data
   toCheck = find(sum(~isnan(bin)) > 0 & sum(~isnan(rest)) > 0);
   pval(arch,toCheck) = arrayfun(@(x) nanranksum(bin(:,x),rest(:,x)), toCheck); 
            
   %EnPQ - Calculate median  diferences
   medianDifference(arch,:) = (nanmedian(bin)-nanmedian(rest));%./(mean([(median(bin));(median(rest))])+10^-6);
   meanDifference(arch,:) = (nanmean(bin)-nanmean(rest));%./(mean([(mean(bin));(mean(rest))])+10^-6);
end
%FDR
isSignificantAfterFDR = mafdr(pval(:),'BHFDR',true)<=ThreshHoldBH; 
%WriteTable
%Cols: ArchetypeNum, Features, pvals, median and mean differences, Significant After FDR 

table =  [repmat((1:Numarchs)',numFeatures,1),...
          reshape(ones(Numarchs,1)*(1:numFeatures),Numarchs*numFeatures,1),...
          pval(:),...
          medianDifference(:),...
          meanDifference(:),...
          isSignificantAfterFDR,...
          PoverQ(:)];

    
    


