function [arc, arcOrig, pc, coefs1] = ParTI_lite(DataPoints,algNum,dim,DiscFeatName,EnMatDis,cols,ContFeatName,EnMatCont,GOcat2Genes,binSize,OutputFileName,arcOrig)
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
% columns, set to 0. To booleanize no columns, set to -1.
% 7. ContFeatName, a cell array of strings listing the Continuous Feature labels.
% 8. EnMatCont, a double matrix with each Continuous feature as a column, 
% and each datapoint as a row. The rows of datapoints should match those in 
% EnMatCont.
% 9. GOcat2Genes, a parameter that is ignored by ParTI_lite(). See the documentation of ParTI().
% 10. binSize, a double variable ranging between 0-1, that is the fraction of the datapoints that should be grouped in a 
% single bin when calculating feature enrichments.
% 11. OutputFileName, a string that represent the name of the comma or tab delimmited file that saves all enrichment
% data. Several files will be created with different ending to specify the continuous enrichment, the discrete enrichment 
% the significant enriched features and all the features. 
% 12. arcOrig, a matrix of archetypes x traits. If this matrix is given,
% archetypes won't be computed and this matrix will be used instead. This
% allows to keep the ordering of archetypes when repeatingly running ParTI.
%
% Output:
% arc, a double matrix of the coordinates of the archetypes in the space spanned by the
% principle components.
% arcOrig, a double matrix of the coordinates of the archetypes in the original space, defined 
% by datapoints.
% pc is a double matrix that saves the pc of the data

% Sisal is presented at Bioucas-Dias JM (2009) in First Workshop on Hyperspectral Image and Signal Processing: Evolution in Remote Sensing, 2009. WHISPERS '09, pp 1-4.
% MVSA is presented at Li J, Bioucas-Dias JM (2008) in Geoscience and Remote Sensing Symposium, 2008. IGARSS 2008. IEEE International, pp III.250-III.253.
% SDVMM and MVES are taken from http://mx.nthu.edu.tw/~tsunghan/Source%20codes.html
% PCHA is taken from http://www.mortenmorup.dk/index_files/Page327.htm
if nargin<12
    arcOrig=[];
end
if nargin<11
    OutputFileName='ParTIliteOutputFile';
end
if nargin<10
    binSize=0.05;
end
if nargin<9
    GOcat2Genes=[];
end
if nargin<8
    EnMatCont=[];
end
if nargin<7
    ContFeatName=[];
end
if nargin<6
    cols=[];
end
if nargin<5
    EnMatDis=[];
end
if nargin<4
    DiscFeatName=[];
end
if nargin<3
    dim=10; %Default dimension is set to 10
end
if nargin<2
    algNum=1; %Default algorithm is Sisal (alg=1)
end
if nargin<1
    error('Too few arguments in the input, define at least the data points')
end

DataPointsSize = size(DataPoints);
%The following prevents the user from specifying more dimensions than
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

maxRuns=0; % number of data randomization, none for lite!
numIter=50; % number of times we search for a simplex on the same dataset
if exist('lowIterations', 'var') && ~isempty(lowIterations)
    maxRuns= 0; %current value for the number of data randomization, none for lite!
    numIter= 5 ; %current value for the number of iterations to run the algorithm
    fprintf('Warning! lowIterations flag set: will only run numIter = %d\n', numIter);
end

if size(arcOrig,1) == 0
    [pc, arc, arcOrig,~,~,coefs1] = findArchetypes(DataPoints,algNum,dim,OutputFileName,numIter,maxRuns);
else
    %We allow passing archetypes so that indexing can stay constant from
    %run to run
    fprintf('Will use archetypes passed as argument instead of determining them from scratch.\n');
    [coefs1,pc] = pca(DataPoints);
    nArch = size(arcOrig,1);
    trArcOrig = bsxfun(@minus,arcOrig,mean(DataPoints));
    
    arc = trArcOrig * coefs1(:,1:(nArch-1));
end

if size(EnMatDis,1) > 0 || size(EnMatCont, 1) > 0
    calculateEnrichment(pc(:,1:size(arc,2)),arc,DiscFeatName,EnMatDis,ContFeatName,EnMatCont,binSize,OutputFileName,[],[],[],[],maxRuns);
else
    fprintf('Skipping enrichment analysis as no features were provided.\n');
end

end

