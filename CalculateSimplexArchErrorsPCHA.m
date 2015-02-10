function minArchsTot=CalculateSimplexArchErrorsPCHA(DataPoints,NumArch,maxRuns,numIter)

% Calculating errors in archetypes location by sampling the data with
% replacement and calculating the archetypes on the 

% Inputs
% DataPoints is the data after PCA with the number of dimensions that
% matches NumArch
% NumArch is the number of archetypes
% NumIter - number of iterations to run the algorithm to check that it got
% the minimal value for this set of Data


dim=NumArch-1; % the dimension of the data
Nmembers=NumArch; % number of Archetypes
NumDataPoints=size(DataPoints,1);

minArchsTot=cell(1,maxRuns); %The archeypes found per all datasets
minRandArchVol=zeros(1,maxRuns); % The volume of each simplex found

for k=1:maxRuns
    if mod(k,round(maxRuns/10)) == 0
        fprintf('%.0f%% done\n', 100*k/maxRuns);
    end
    %Initialization of arrays that contain the archetypes and their volume
    minArchs=cell(1,numIter);
    VolArch1=zeros(1,numIter);
        
    %sample with replacement the data points
    indSampled=randsample(1:NumDataPoints,NumDataPoints,true); %create a vector with the indices to sample
    
    % calculate the archetypes for a specific dataset
    for runs=1:numIter
        delta = 0.5;
        U=1:size(DataPoints,1); % Entries in X used that is modelled by the AA model
        I=1:size(DataPoints,1); % Entries in X used to define archetypes
        [Archs,~,~,~,varexpl]=PCHA1(DataPoints(indSampled,1:dim)',Nmembers,I,U,delta);
       
        % calculating volumes
        % PCHA algorithm
        if ~isnan(Archs)
            Arch1Red=bsxfun(@minus,Archs,Archs(:,Nmembers));
            VolArch1(runs)=abs(det(Arch1Red(:,1:end-1))/factorial(Nmembers-1));
            minArchs{runs}=Archs(:,1:Nmembers);
        else
            VolArch1(runs)= NaN;
            minArchs{runs}=NaN;
        end
        
    end
    [minRandArchVol(k) IndminVol]=max(VolArch1);
    minArchsTot{k}=minArchs{IndminVol};
end

    
