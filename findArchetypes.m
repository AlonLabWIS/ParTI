function [DataPCA,meanClstErrs,realArchs,ArchsErrors,PvalueRatio,coefs1]=findArchetypes(DataPoints,algNum,dim,OutputFileName,numIter,maxRuns)
%Inputs
% 1. Data points is the values of different traits (e.g. expression
% level of genes) - each sample is a row, each trait (gene) is a column
% 2. algNum is for choosing the algorithm to find the simplex:
%    algNum=1 :> Sisal (default)
%    algNum=2 :> MVSA
%    algNum=3 :> MVES
%    algNum=4 :> SDVMM
%    algNum=5 :> PCHA
%
% Sisal is presented at Bioucas-Dias JM (2009) in First Workshop on Hyperspectral Image and Signal Processing: Evolution in Remote Sensing, 2009. WHISPERS �09, pp 1�4.
% MVSA is presented at Li J, Bioucas-Dias JM (2008) in Geoscience and Remote Sensing Symposium, 2008. IGARSS 2008. IEEE International, pp III � 250�III � 253.
% SDVMM and MVES are taken from http://mx.nthu.edu.tw/~tsunghan/Source%20codes.html
% PCHA is taken from http://www.mortenmorup.dk/index_files/Page327.htm
% 3. dim is the dimension up to which dimension should the ESV be calculated
% 4. OutputFileName, a tag that will be used in the name of the files in which
% figures will be saved.
% 5. numIter, the number of iterations for the algorithm to find a minimal
% bounding simplex
% 6. maxRuns, how many bootstrapping iterations should be performed to
% estimate the errors on the archetypes
%
%Outputs
% 1. The data after PCA
% 2. Postition of the archetypes in the low-dim space of PCA
% 4. The calculated archetypes in the Original coordinates
% 3. Covariance matrix of the errors on the archetypes obtained by bootstrapping
% 5. The statistical significance of the simplex (p-value at the t-ratio test)

addpath(genpath(pwd)); %Add all subfolders of the current directory to run all the diff. algorithms
global ForceNArchetypes;
global abortAfterPval;
PvalueRatio = [];

% Initialize the parameters
if nargin<2
    algNum=1;
    dim=10;
    DimFig=3;
else if nargin<3
        dim=10;
        DimFig=3;
    else if nargin<4
            DimFig=3;
        end
    end
end


%% Do PCA on the data
fprintf('Starting to perform PCA, for big data on slow computers this may take a while...\n');
if exist('PCA') == 0
    error('This package requires the pca() function from the Matlab Statistics and Machine Learning Toolbox. Please install the Matlab Statistics and Machine Learning Toolbox.');
end
[coefs1,scores1,variances] = pca(DataPoints);
percent_explained = 100*cumsum(variances)/sum(variances);

figure;
plot(percent_explained(1:dim),'.-','linewidth',2,'MarkerSize',20);
title('Cumulative variability explained per principle component','fontsize',14);
xlabel('Dimension','fontsize',14);ylabel('% variability explained','fontsize',14);
if exist('savefig')
    savefig([OutputFileName,'_CumVarExpPCA.fig']);
end

figure;
plot(1:dim,diff([0, percent_explained(1:dim)']),'.-','linewidth',2,'MarkerSize',20);
title('Variability explained by each principle component','fontsize',14);
xlabel('Dimension','fontsize',14);ylabel('% variability explained','fontsize',14);

DataPCA=scores1;

%% Calculate ESV for dimensions 2-min(10,Dimension of the data)

DataPCA_centered=bsxfun(@minus,DataPCA,mean(DataPCA,1));

% disp('calculating ESV');

TotESV1=zeros(1,dim);
varexpl=zeros(1,dim);
%calculate ESV using the standard PCHA method.
fprintf('Calculating explained variance with PCHA (Morup M, Hansen KL, 2011)\n');

for indNmembers=1:dim
    Nmembers=indNmembers+1; % number of Archetypes is dimension+1
    
    delta = 0;
    U=1:size(DataPCA,1); % Entries in X used that is modelled by the AA model
    I=1:size(DataPCA,1); % Entries in X used to define archetypes
    
    [~,~,~,~,varexpl(indNmembers)]=PCHA1(DataPCA(:,1:dim)',Nmembers,I,U,delta);
    % calculate the archtypes with SDVMM
    %     [Arch12,~] = SDVMM(DataPCA_centered',Nmembers,0);
    %
    %     %initialize the ESV values
    %     ESV1=zeros(1,size(DataPCA_centered,1));
    %
    %     %calculate ESV values
    %     for i=1:size(DataPCA_centered,1)
    %         distPerp = norm(DataPCA_centered(i,Nmembers:end));
    %         distance1=distFromPoly(DataPCA_centered(i,1:Nmembers-1),Arch12(1:Nmembers-1,:));
    %         ESV1(i)=1-(sqrt(distance1^2+distPerp^2))/norm(DataPCA_centered(i,:));
    %     end
    
    %TotESV1 is the ESV value per dimension
    %    TotESV1(indNmembers)=1/length(ESV1)*sum(ESV1);
end
TotESV1 = varexpl*(sum(variances(1:dim))/sum(variances));
% plot the ESV curve to extract the desired dimension
figure;
plot(2:dim+1,100*TotESV1,'.-','linewidth',2,'MarkerSize',20);
title('ESV for different numbers of archetypes','fontsize',14);
xlabel('Number of archetypes','fontsize',14);ylabel('% variability explained','fontsize',14);
if exist('savefig') == 2
    savefig([OutputFileName,'_cumESV.fig']);
end

figure;
plot(2:dim+1,diff([0 100*TotESV1]),'.-','linewidth',2,'MarkerSize',20);
title('ESV for each archetype','fontsize',14);
xlabel('Number of archetypes','fontsize',14);ylabel('% variability explained','fontsize',14);
if exist('savefig')
    savefig([OutputFileName,'_ESV.fig']);
end

%% Get the desired dimension from the user
if exist('ForceNArchetypes','var') && ~isempty(ForceNArchetypes)
    fprintf('Warning! ForceNArchetypes preset in workspace to %d. Will now use that value.\n', ForceNArchetypes);
    NArchetypes=ForceNArchetypes;
else
    if dim > 2
        SuggNArchetypes = DimensionFinder(TotESV1) + 1;
        NArchetypes=input(['Elbow method suggests ' num2str(SuggNArchetypes) ...
            ' archetypes.\nPlease indicate the desired number of archetypes (or press enter for using the suggestion): ']);
        if isempty(NArchetypes)
            NArchetypes =  SuggNArchetypes;
        end
    else
        NArchetypes=input('Please indicate the desired number of archetypes: ');
    end
end

if NArchetypes<3 && algNum >=3
    msgbox('Only sisal and MVSA algorithms allow you to run less than 3 archetypes!','Error','error');
    error ('Only sisal and MVSA algorithms allow you to run less than 3 archetypes!');
end

%% Find the archetypes of the bounding simplex in d-dimensions
DataDim=size(DataPCA,2);
%We need to figure out how to generalize volumes to non-simplical polytopes
%before we allow running PCHA with NArchetypes>DataDim+1
%if (algNum~=5)
    if (NArchetypes>DataDim+1) %meaning NumArchetypes>dim+1
        fprintf('Warning! Number of Archetypes (%d) exceeds data dimensions (%d) + 1\nWe reset to %d archetypes.\n', ...
            NArchetypes, DataDim, DataDim+1);
        NArchetypes = DataDim+1;
        %return;
    end
    if (NArchetypes>DataDim) %meaning NumArchetypes=dim+1
        DataPCA=[DataPCA, ones(size(DataPCA,1),1)]; %embedding the data in a D+1 space
    end
%end

[ArchsMin,VolArchReal]=findMinSimplex(numIter,DataPCA,algNum,NArchetypes);
%Now we'll re-order archetypes according to their coordinate on the first
%PC. This should help achieve some reproducibility from ParTI run to ParTI
%run.
[~,ArchsOrder] = sort(ArchsMin(1,:));
ArchsMin = ArchsMin(:,ArchsOrder);
disp('finished finding the archetypes');

if NArchetypes < 4
    DimFig = 2;
else
    DimFig = 3;
end

Xeltot=cell(1,NArchetypes);
Yeltot=cell(1,NArchetypes);
Zeltot=cell(1,NArchetypes);

if maxRuns > 0
	%% Calculate the p-value from t-ratios

	% calculate the volume of the convex hull of the real dataNArchetypes-1
	ConHullVol = ConvexHull(DataPCA(:,1:min(NArchetypes-1,DataDim)));
	%[~ , ConHullVol]=convhulln(DataPCA(:,1:NArchetypes-1));

	% calculate the t-ratio for the real data
	tRatioReal=VolArchReal/ConHullVol;

	% run a function that gets the data and outputs the bootstraped shuffled data t-ratios
	% when high number of points exists, Sisal is recommended. If number of
	% points is relatively low, SDVMM is recommended.
	fprintf('Now computing t-ratios.\n');

    if NArchetypes > 2
        switch algNum
            case 1 %    algNum=1 :> Sisal (default)
    	        tRatioRand = CalculateSimplexTratiosSisal(DataPCA(:,1:NArchetypes),NArchetypes,maxRuns,numIter);
    	    case 2 %    algNum=2 :> MVSA
    	        tRatioRand = CalculateSimplexTratiosMVSA(DataPCA(:,1:NArchetypes),NArchetypes,maxRuns,numIter);
    	    case 3 %    algNum=3 :> MVES
    	        tRatioRand = CalculateSimplexTratiosMVES(DataPCA(:,1:NArchetypes-1),NArchetypes,maxRuns,numIter);
    	    case 4 %    algNum=4 :> SDVMM
    	        tRatioRand = CalculateSimplexTratiosSDVMM(DataPCA(:,1:NArchetypes-1),NArchetypes,maxRuns,numIter);
    	    case 5 %    algNum=5 :> PCHA
    	        tRatioRand = CalculateSimplexTratiosPCHA(DataPCA(:,1:min(NArchetypes-1,DataDim)),NArchetypes,maxRuns,numIter);
        end
    
    	if algNum<4 %for first three methods the t-ratio is larger than 1, and you want the minimal one
    	    PvalueRatio=sum(tRatioRand<tRatioReal)/maxRuns;
    	else %for last two methods (SDVMM, PCHA) t-ratio is smaller than 1, and you want the maximal one
    	    PvalueRatio=sum(tRatioRand>tRatioReal)/maxRuns;
        end
    else
        %if there is only 2 archetypes, compute the ratio of first two
        %eigenvalues as a statistic for how line-like the data is.
        %it is necessary to use a different statistic compared to the 3
        %archetype case because the volume of a line is 0.
        
        asymReal = variances(1) / variances(2);

        asymsRand=zeros(1,maxRuns); % This will hold the volume of the convex hull
        for i=1:maxRuns
            if mod(i,round(maxRuns/10)) == 0
                fprintf('%.0f%% done\n', 100*i/maxRuns);
            end
            
            SimplexRand=zeros(size(DataPoints));
            % Shuffle the Sampled data -
            for j=1:size(DataPoints,2) % for each dimension
                shuffInd=randperm(size(DataPoints,1)); % shuffle the data values of each axis
                SimplexRand(:,j)=DataPoints(shuffInd,j);
            end
            [~,~,variancesRand] = pca(SimplexRand);
            asymsRand(i) = variancesRand(1)/variancesRand(2);
        end
        PvalueRatio = sum(asymsRand>asymReal)/maxRuns;
    end

	fprintf('The significance of %d archetypes has p-value of: %2.5f \n',NArchetypes,PvalueRatio);
end

if maxRuns > 0 && ~(exist('abortAfterPval','var') && ~isempty(abortAfterPval))
	%% Calculate errors in archetypes (by bootstrapping)
    fprintf('Now calculating errors on the archetypes.\n');
    switch algNum
        case 1 %    algNum=1 :> Sisal (default)
            ArchsErrors = CalculateSimplexArchErrorsSisal(DataPCA(:,1:NArchetypes),NArchetypes,maxRuns,numIter);
        case 2 %    algNum=2 :> MVSA
            ArchsErrors = CalculateSimplexArchErrorsMVSA(DataPCA(:,1:NArchetypes),NArchetypes,maxRuns,numIter);
        case 3 %    algNum=3 :> MVES
            ArchsErrors = CalculateSimplexArchErrorsMVES(DataPCA(:,1:NArchetypes-1),NArchetypes,maxRuns,numIter);
        case 4 %    algNum=4 :> SDVMM
            ArchsErrors = CalculateSimplexArchErrorsSDVMM(DataPCA(:,1:NArchetypes-1),NArchetypes,maxRuns,numIter);
        case 5 %    algNum=5 :> PCHA
            ArchsErrors = CalculateSimplexArchErrorsPCHA(DataPCA(:,1:min(NArchetypes-1,DataDim)),NArchetypes,maxRuns,numIter);
    end

    % create the error clouds per archetype
    % cluster the archtypes to clouds
    ArchsErrorsMat=cell2mat(ArchsErrors)';
    switch NArchetypes
        case 1
            clusteredArchsErrorInd = ones(size(ArchsErrorsMat,1),1);
        case 2
            clusteredArchsErrorInd = (ArchsErrorsMat > 0) + 1;
        otherwise
            clusteredArchsErrorInd = kmeans(ArchsErrorsMat,NArchetypes,'distance','cosine','replicates',10);        
    end

    meanClstErrs=zeros(NArchetypes,NArchetypes-1);
    if NArchetypes < 3
        meanClstErrs(:,2) = zeros(size(meanClstErrs,1),1);
        if NArchetypes < 2
            meanClstErrs(:,1) = zeros(size(meanClstErrs,1),1);
        end
    end
    archReordering = zeros(1,NArchetypes);
    for l=1:NArchetypes
        clusteredArchsError=ArchsErrorsMat(clusteredArchsErrorInd==l,:);
        if NArchetypes < 3
            clusteredArchsError(l,2) = 0;
            if NArchetypes < 2
                clusteredArchsError(l,1) = 0;
            end
        end
        meanClstErrs(l,:)=mean(clusteredArchsError);
        repArchCoord = repmat(meanClstErrs(l,:), NArchetypes, 1);
        [~,archReordering(l)] = min(sum((repArchCoord - ArchsMin').^2,2));
    end

    if length(unique(archReordering)) < NArchetypes
        archReordering = 1:NArchetypes;
        disp('Warning: could not align archetypes. Archetype order is random and will change if you rerun ParTI again.');
    end
    %[~,archReordering]=sort(archReordering);
    %[archReordering(clusteredArchsErrorInd)',clusteredArchsErrorInd]
    clusteredArchsErrorIndAligned = archReordering(clusteredArchsErrorInd)';

    El1=zeros(1,NArchetypes);
    El2=zeros(1,NArchetypes);
    phi=zeros(1,NArchetypes);
    Coeff2d = cell(1,NArchetypes);
    RotEllipsoidArch = cell(1,NArchetypes);
    for l=1:NArchetypes
        clusteredArchsError=ArchsErrorsMat(clusteredArchsErrorIndAligned==l,:);
        if NArchetypes < 3
            clusteredArchsError(l,2) = 0;
            if NArchetypes < 2
                clusteredArchsError(l,1) = 0;
            end
        end
        meanClstErrs(l,:)=mean(clusteredArchsError);

        % remove the mean of each column - move the errors to zero
        clstArchErrMeanless=bsxfun(@minus,clusteredArchsError(:,1:DimFig),meanClstErrs(l,1:DimFig));
        % calculating the axes of the principal components
        [Coeff2d{l},~,loadings2d]=pca(clstArchErrMeanless(:,1:2));
        El1(l) = loadings2d(1)^(1/2);
        El2(l) = loadings2d(2)^(1/2);

        if DimFig >= 3
            [Coeff,~,loadings]=pca(clstArchErrMeanless);
            % generate the ellipsoid
            [Xel,Yel,Zel]=ellipsoid(0,0,0,loadings(1)^(1/2),loadings(2)^(1/2),loadings(3)^(1/2),25);
            % move the ellipsoid to the archtype location and rotate the ellipsoid
            % to its principal axes
            RotEllipsoidArch{l}=arrayfun(@(x,y,z) Coeff*[x,y,z]',Xel,Yel,Zel,'uniformoutput',0);
            RotEllipMat=cell2mat(RotEllipsoidArch{l});
            Xeltot{l}=meanClstErrs(l,1)+RotEllipMat(1:3:end,:);
            Yeltot{l}=meanClstErrs(l,2)+RotEllipMat(2:3:end,:);
            Zeltot{l}=meanClstErrs(l,3)+RotEllipMat(3:3:end,:);
        end
    end
    if NArchetypes < 3
        meanlessTemp = meanClstErrs(:,1:NArchetypes-1)*(coefs1(:,1:NArchetypes-1)');
    else
        meanlessTemp = meanClstErrs*(coefs1(:,1:NArchetypes-1)');
    end
    realArchs = bsxfun(@plus,meanlessTemp,mean(DataPoints));

    ArchsErrors = cell(1,NArchetypes);
    %vol = abs(det(bsxfun(@minus,meanClstErrs(1:end-1,:),meanClstErrs(end,:)))...
    %    /factorial(NArchetypes-1));
    for l=1:NArchetypes
        clusteredArchsError=ArchsErrorsMat(clusteredArchsErrorIndAligned==l,:);
        if NArchetypes < 3
            clusteredArchsError(l,2) = 0;
            if NArchetypes < 2
                clusteredArchsError(l,1) = 0;
            end
        end
         %[~,~,loadingsArc]=pca(clusteredArchsError);
         %volarc = exp(mean(log(loadingsArc)));
         ArchsErrors{l} = cov(clusteredArchsError); %volarc / power(vol,1/(NArchetypes-1));
    end

    disp('finished finding the archetypes error distribution'); 
else
	meanClstErrs = ArchsMin';
	meanlessTemp = meanClstErrs*(coefs1(:,1:NArchetypes-1)');
	realArchs = bsxfun(@plus,meanlessTemp,mean(DataPoints));
	ArchsErrors = [];
	%PvalueRatio = [];
    % At this point, we have all we need to return for the 'lite' version of this function
    
    %We'll just compute 'dot' sizes for the archetypes and we're done:
    tmp=meanClstErrs;
    if NArchetypes < 3
        tmp(:,2) = zeros(size(tmp,1),1);
        if NArchetypes < 2
           tmp(:,1) = zeros(size(tmp,1),1);
        end
    end
    p = 0.02 * norm(tmp);
    for l=1:NArchetypes
        % generate the ellipsoid
        [Xel,Yel,Zel]= sphere;
        % move the ellipsoid to the archtype location and rotate the ellipsoid
        % to its principal axes
        RotEllipsoidArch=arrayfun(@(x,y,z) [p * x,p* y,p* z]',Xel,Yel,Zel,'uniformoutput',0);
        RotEllipMat=cell2mat(RotEllipsoidArch);
        Xeltot{l}=tmp(l,1)+RotEllipMat(1:3:end,:);
        Yeltot{l}=tmp(l,2)+RotEllipMat(2:3:end,:);

        if DimFig >= 3
            Zeltot{l}=tmp(l,3)+RotEllipMat(3:3:end,:);
        end
    end
end

%% Jean> Quick hack so we can work with two archetypes
if NArchetypes == 2
    meanClstErrs(:,2) = 0;
end
%% End hack

% plotting the data in the first 2 PC's
styleel={'-r','-g','-b','-m','-y','-c','-k','--r','--b','--g','--k','--m','--c','--y'};
style={'.r','.g','.b','.m','.y','.c','.k','or','ob','og','ok','om','oc','oy'};
cmap=[1 0 0;
    0 1 0;
    0 0 1;
    1 0 1;
    1 1 0;
    0 1 1;
    0 0 0;];

figure;
% plot the data points
plot(DataPCA(:,1),DataPCA(:,2),'.k');
hold on;
% plot the archetypes in 2d
for arcCol = 1:NArchetypes
	if maxRuns > 0 && ~(exist('abortAfterPval','var') && ~isempty(abortAfterPval))
        %Only show errors if we actually asked to compute them
	    ellipse(meanClstErrs(arcCol,1),meanClstErrs(arcCol,2),El1(arcCol),...
	        El2(arcCol),Coeff2d{arcCol},styleel{mod(arcCol-1,14)+1});
	else
		plot(meanClstErrs(arcCol,1),meanClstErrs(arcCol,2),style{mod(arcCol-1,14)+1},'markersize',35);
	end
    text(meanClstErrs(arcCol,1),meanClstErrs(arcCol,2),['    ', num2str(arcCol)],'FontSize',15);
end
axis equal
xlabel('PC1','fontsize',14);ylabel('PC2','fontsize',14);
if exist('savefig')
    savefig([OutputFileName,'_ArchsIn2D.fig']);
end

if (DimFig == 3)
    % plotting the data in the first 3 PC's
    figure;
    % plot the data points
    plot3(DataPCA(:,1),DataPCA(:,2),DataPCA(:,3),'.k');
    hold on;
    % plot the error cloud per archetype
    for i=1:NArchetypes
        surf(Xeltot{i},Yeltot{i},Zeltot{i},'EdgeColor', 'none');
        colormap(cmap(mod(i-1,7)+1,:));
        freezeColors;
        text(meanClstErrs(i,1),meanClstErrs(i,2),meanClstErrs(i,3),['   ',num2str(i)],'fontsize',15);
    end
    axis equal
    box on
    xlabel('PC1','fontsize',14);
	ylabel('PC2','fontsize',14);
    zlabel('PC3','fontsize',14);
    if exist('savefig')
        savefig([OutputFileName,'_ArchsIn3D.fig']);
    end
end

% Remove the embbed ones before we return the PCA projected data
DataPCA = DataPCA(:,1:size(scores1,2));

end


