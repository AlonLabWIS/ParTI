function [arc, arcOrig, pc] = PAAM_lite(DataPoints,algNum,dim,DiscFeatName,EnMatDis,cols,ContFeatName,EnMatCont,binSize,OutputFileName)
%% Inputs
% 1. DataPoints, double matrix with the values of different traits 
% (the coordinates, e.g. expression level of genes). Each sample is a row, 
% each trait (e.g. each gene) is a column.
% 2. algNum is an integer that chooses the algorithm to find the simplex:
%    algNum=1 :> Sisal (default)
%    algNum=2 :> MVSA 
%    algNum=3 :> MVES
%    algNum=4 :> SDVMM 
%    algNum=5 :> PCHA 
% 3. dim, an integer representing the dimension up to which the Explained Variance (ESV) 
% should be calculated. ESV is defined as the average distance between the 
% data points and the best simplex defined by the archetypes found by
% PCHA. We recommend to start with a dimension between 5 - 10, and 
%  increase the dimension until an inflexion be spotted the ESV curve.
% 4. DiscFeatName, a cell array of strings listing the Discrete Feature labels (row)  
% 5. EnMatDis, a boolean matrix with each Discrete feature as a column, 
% and each datapoint as a row. Note that the rows of the datapoints matrix 
% should match the rows of EnMatDis. In addition, the order of the columns 
% in EnMatDis should match that of DiscFeatName. it can also accept a cell
% array of strings and transform it into a boolean matrix.   
% 6. cols, The integer index of the columns of EnMatDis to booleanize, in case
% EnMatDis is a string cell array rather than a boolean matrix. To booleanize all
% columns, set to 0.
% 7. ContFeatName, a cell array of strings listing the Continuous Feature labels.
% 8. EnMatCont, a double matrix with each Continuous feature as a column, 
% and each datapoint as a row. The rows of datapoints should match those in 
% EnMatCont.
% 9. binSize, a double variable ranging between 0-1, that is the fraction of the datapoints that should be grouped in a 
% single bin when calculating feature enrichments.
% 10. OutputFileName, a string that represent the name of the comma or tab delimmited file that saves all enrichment
% data. Several files will be created with different ending to specify the continuous enrichment, the discrete enrichment 
% the significant enriched features and all the features. 
%
% Output:
% arc, a double matrix of the coordinates of the archetypes in the space spanned by the
% principle components.
% arcOrig, a double matrix of the coordinates of the archetypes in the original space, defined 
% by datapoints. 
% pc is a double matrix that saves the pc of the data 

% Sisal is presented at Bioucas-Dias JM (2009) in First Workshop on Hyperspectral Image and Signal Processing: Evolution in Remote Sensing, 2009. WHISPERS �09, pp 1�4.
% MVSA is presented at Li J, Bioucas-Dias JM (2008) in Geoscience and Remote Sensing Symposium, 2008. IGARSS 2008. IEEE International, pp III � 250�III � 253.
% SDVMM and MVES are taken from http://mx.nthu.edu.tw/~tsunghan/Source%20codes.html
% PCHA is taken from http://www.mortenmorup.dk/index_files/Page327.htm


DataPointsSize = size(DataPoints);
%The following makes prevents the user from specifying more dimensions than
%the dataset has (Jean)
if ( DataPointsSize(2) < dim )
    dim = DataPointsSize(2);
    fprintf('Warning: the maximal dimension was set to %d, which is the total dimensionality of the data \n', dim);
end


if iscell(EnMatDis) || ~isequal(unique(EnMatDis),[0;1])
    if cols>0
        [EnMatDis, DiscFeatName] = DiscreteToBoolean(EnMatDis, DiscFeatName, cols);
    else
        [EnMatDis, DiscFeatName] = DiscreteToBoolean(EnMatDis, DiscFeatName);
    end
end
% Initializing the running algorithm parameters
global lowIterations;

%maxRuns=1000; % current value for the number of data randomization
numIter=50; % current value for the number of iterations to run the algorithm
if exist('lowIterations', 'var') && ~isempty(lowIterations)
    %maxRuns= 20; %current value for the number of data randomization
    numIter= 5 ; %current value for the number of iterations to run the algorithm
    fprintf('Warning! lowIterations flag set: will only run numIter = %d\n', numIter);
end

[pc , arc, arcOrig] = findArchetypes_lite(DataPoints,algNum,dim,OutputFileName,numIter);

calculateEnrichment_lite(pc(:,1:size(arc,2)),arc,DiscFeatName,EnMatDis,ContFeatName,EnMatCont,binSize,OutputFileName);

end

