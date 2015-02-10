function minRandArchRatio=CalculateSimplexTratiosPCHA(DataPoints,NumArch,maxRuns,numIter)

% DataPoints is the data after PCA
% NumArch is the number of archetypes
% maxRuns is the number of shuffling in the randomization bootstrapping
% process, default set to 10000
% numIter - number of iterations per dataset run

% Calculating P-values of triangles and tetrahedrons
dim=NumArch-1; % dimension
Nmembers=NumArch; % number of Archetypes
% maxRuns=10000; %times the data is shuffled to get a P-value (10000)
% numIter=50; % times that algorithm runs to get the best fit (50)

minRandArchVol=zeros(maxRuns,1);
minRandArchRatio=zeros(maxRuns,1);   

% Calculate the P-Value: running 10,000 times of shuffled data to get
% the tRatios for a specific simplex in D dimensions 

VolConvRand=zeros(1,maxRuns); % This will hold the volume of the convex hull

for m=1:maxRuns
    if mod(m,round(maxRuns/10)) == 0
        fprintf('%.0f%% done\n', 100*m/maxRuns);
    end    
%     if mod(m,10)==0
%         disp(m);
%     end
    SimplexRand1=zeros(size(DataPoints));
    % Shuffle the Sampled data -
    for i=1:size(DataPoints,2) % for each dimension
        shuffInd=randperm(size(DataPoints,1)); % shuffle the data values of each axis
        SimplexRand1(:,i)=DataPoints(shuffInd,i);
    end
    
   % [~ , VolConvRand(m)]=convhulln(SimplexRand1(:,1:Nmembers-1));
    VolConvRand(m) = ConvexHull(SimplexRand1(:,1:Nmembers-1));
    VolArchRand=zeros(1,numIter);
    RandDataRatios=zeros(1,numIter);
    for k=1:numIter
        % PCHA algorithm
        delta = 0;
        U=1:size(SimplexRand1,1); % Entries in X used that is modelled by the AA model
        I=1:size(SimplexRand1,1); % Entries in X used to define archetypes
        [Arch3Rand,~,~,~,varexpl]=PCHA1(SimplexRand1',Nmembers,I,U,delta);
        % Calculating the volume of the randomized simplex
        if ~isnan(Arch3Rand)
            ArchRandRed=bsxfun(@minus,Arch3Rand,Arch3Rand(:,Nmembers));
            VolArchRand(k)=abs(det(ArchRandRed(:,1:end-1))/factorial(Nmembers-1));
            RandDataRatios(k)=VolArchRand(k)./VolConvRand(m);
        else
            VolArchRand(k)= NaN;
            RandDataRatios(k)=NaN;
        end
    end
    
    minRandArchVol(m)=max(VolArchRand);
    minRandArchRatio(m)=max(RandDataRatios);
    
end