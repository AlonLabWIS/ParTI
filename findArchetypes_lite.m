function [DataPCA,ArchsFinal,realArchs]=findArchetypes_lite(DataPoints,algNum,dim,OutputFileName,numIter)

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
%
%Outputs 
% 1. The data after PCA
% 2. The calculated archetypes in the PCA coordinates
% 3. The calculated archetypes in the Original coordinates 

addpath(genpath(pwd)); %Add all subfolders of the current directory to run all the diff. algorithms
global ForceNArchetypes;

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
if exist('princomp') == 0
    error('This package requires the princomp() function from the Matlab Statistical Toolbox. Please install the Matlab Statistical Toolbox or provide an implementation of the princomp() function.');
end
[coefs1,scores1,variances] = princomp(DataPoints);
percent_explained = 100*cumsum(variances)/sum(variances);

figure;
plot(percent_explained(1:dim),'.-','linewidth',2,'MarkerSize',20);
title('Cumulative variability explained per principle component','fontsize',14);
xlabel('Dimension','fontsize',14);ylabel('% variability explained','fontsize',14);
if exist('savefig')
    savefig([OutputFileName,'1_CumVarExpPCA.fig']);
end

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
title('ESV for different dimensions','fontsize',14);
xlabel('Number of Archetypes','fontsize',14);ylabel('% variability explained','fontsize',14);
if exist('savefig')
    savefig([OutputFileName,'2_ESV.fig']);
end

%% Get the desired dimension from the user
if exist('ForceNArchetypes','var') && ~isempty(ForceNArchetypes)
    fprintf('Warning! ForceNArchetypes preset in workspace to %d. Will now use that value.\n', ForceNArchetypes);
    NArchetypes=ForceNArchetypes;
else
    NArchetypes=input('please indicate the desired number of archetypes: ');
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

disp('finished finding the archetypes');


ArchsFinal = ArchsMin';
meanlessTemp = ArchsFinal*(coefs1(:,1:NArchetypes-1)');
realArchs = bsxfun(@plus,meanlessTemp,mean(DataPoints));

%% Calculate errors in archetypes (by bootstrapping)

if NArchetypes < 4
    DimFig = 2;
else
    DimFig = 3;
end

Xeltot=cell(1,NArchetypes);
Yeltot=cell(1,NArchetypes);
Zeltot=cell(1,NArchetypes);

meanClstErrs=ArchsMin';
 
 if NArchetypes < 3
        meanClstErrs(:,2) = zeros(size(meanClstErrs,1),1);
        if NArchetypes < 2
           meanClstErrs(:,1) = zeros(size(meanClstErrs,1),1);
        end
 end
  p = 0.02 * norm(meanClstErrs);
for l=1:NArchetypes
    % generate the ellipsoid
   [Xel,Yel,Zel]= sphere;
    % move the ellipsoid to the archtype location and rotate the ellipsoid
    % to its principal axes
    RotEllipsoidArch=arrayfun(@(x,y,z) [p * x,p* y,p* z]',Xel,Yel,Zel,'uniformoutput',0);
    RotEllipMat=cell2mat(RotEllipsoidArch);
    Xeltot{l}=meanClstErrs(l,1)+RotEllipMat(1:3:end,:);
    Yeltot{l}=meanClstErrs(l,2)+RotEllipMat(2:3:end,:);
    
    if DimFig >= 3
        Zeltot{l}=meanClstErrs(l,3)+RotEllipMat(3:3:end,:);
    end
end

%%disp('finished finding the archetypes error distribution');

%% Plot the Data + Archetypes in 3 first PC's + Error clouds of the archetypes
% switch DimFig
%     case 2
        % plotting the data in the first 2 PC's
        figure;
        % plot the data points
 style={'.r','.g','.b','.m','.y','.c','.k','or','ob','og','ok','om','oc','oy'};

        plot(DataPCA(:,1),DataPCA(:,2),'.k');
        hold on;
        % plot the archetypes in 2d
        for arcCol = 1:NArchetypes
         plot(meanClstErrs(arcCol,1),meanClstErrs(arcCol,2),style{mod(arcCol-1,14)+1},'markersize',35);
         text(meanClstErrs(arcCol,1),meanClstErrs(arcCol,2),['   ', num2str(arcCol)]); 
        end
        axis equal
        xlabel('PC1','fontsize',14);ylabel('PC2','fontsize',14); 
	if exist('savefig')
            savefig([OutputFileName,'3_ArchsIn2D.fig']);
	end

   if (DimFig == 3)
        cmap=[1 0 0; 
              0 1 0; 
              0 0 1;
              1 0 1;
              1 1 0;
              0 1 1;
              0 0 0;];
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
         xlabel('PC1','fontsize',14);ylabel('PC2','fontsize',14); 
         zlabel('PC3','fontsize',14); 
	if exist('savefig')
            savefig([OutputFileName,'4_ArchsIn3D.fig']);
	end
    end
end


