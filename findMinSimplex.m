function [ArchsMin,VolArchReal]=findMinSimplex(numIter,DataPCA,algNum,NArchetypes,silent)

% numIter - number of iterations for the algorithm to find a minimal
% bounding simplex
% DataPCA - Data points after PCA
% algNum - which algorithm to run in order to find the bounding simplex:
%    algNum=1 :> Sisal (default)
%    algNum=2 :> MVSA
%    algNum=3 :> MVES
%    algNum=4 :> SDVMM
%    algNum=5 :> PCHA
% NArchetypes - Number of archetypes to find (dim+1)
% silent - no information on the used algorithm will be printed out if 1
%
% Output
% 1. Best fitting simplex
% 2. Simplex volume

if nargin<5
    silent=0;
end

minArchsIter=cell(1,3*numIter);
VolArch=zeros(1,3*numIter);
DataDim=size(DataPCA,2);

%% Jean> I feel like this section was already implemented in findArchetypes.m
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
%% End

switch algNum
    case 1 %    algNum=1 :> Sisal (default)
        if ~silent
            fprintf('Calculating archetypes positions with SISAL (Bioucas-Dias JM, 2009)\n');
        end
        
        for i=1:3*numIter
            [Archs,~, ~, ~] = sisal(DataPCA(:,1:NArchetypes)',NArchetypes,'VERBOSE',0);
            if ~isnan(Archs)
                %calculate the volume of the simplex
                Arch1Red=bsxfun(@minus,Archs,Archs(:,NArchetypes));
                VolArch(i)=abs(det(Arch1Red(1:end-1,1:end-1))/factorial(NArchetypes-1));
                %save the archetypes
                minArchsIter{i}=Archs(1:end-1,1:NArchetypes);
            else
                VolArch(i)= NaN;
                minArchsIter{i}=NaN;
            end
        end
        [VolArchReal IndminVol]=min(VolArch); %find the minimal volume simplex
        ArchsMin=minArchsIter{IndminVol}; %get the minimal archetypes that form this simplex
        
    case 2 %    algNum=2 :> MVSA
        if ~silent
            fprintf('Calculating archetypes positions with MVSA (Li J, Bioucas-Dias JM, 2008)\n');
        end
        
        for i=1:3*numIter
            [Archs,~, ~, ~] = mvsa(DataPCA(:,1:NArchetypes)',NArchetypes,'verbose',0);
            if ~isnan(Archs)
                %calculate the volume of the simplex
                Arch1Red=bsxfun(@minus,Archs,Archs(:,NArchetypes));
                VolArch(i)=abs(det(Arch1Red(1:end-1,1:end-1))/factorial(NArchetypes-1));
                %save the archetypes
                minArchsIter{i}=Archs(1:end-1,1:NArchetypes);
            else
                VolArch(i)= NaN;
                minArchsIter{i}=NaN;
            end
        end
        [VolArchReal IndminVol]=min(VolArch); %find the minimal volume simplex
        ArchsMin=minArchsIter{IndminVol}; %get the minimal archetypes that form this simplex
        
        
    case 3 %    algNum=3 :> MVES
        if ~silent
            fprintf('Calculating archetypes positions with MVES (Chan T-H, Chi C-Y, Huang Y-M and Ma W-K, 2009)\n');
        end
        
        for i=1:3*numIter
            [Archs,~, ~, ~] = MVES(DataPCA(:,1:NArchetypes-1)',NArchetypes,0);
            if ~isnan(Archs)
                %calculate the volume of the simplex
                Arch1Red=bsxfun(@minus,Archs,Archs(:,NArchetypes));
                VolArch(i)=abs(det(Arch1Red(:,1:end-1))/factorial(NArchetypes-1));
                %save the archetypes
                minArchsIter{i}=Archs;
            else
                VolArch(i)= NaN;
                minArchsIter{i}=NaN;
            end
        end
        [VolArchReal IndminVol]=min(VolArch); %find the minimal volume simplex
        ArchsMin=minArchsIter{IndminVol}; %get the minimal archetypes that form this simplex
        
    case 4 %    algNum=4 :> SDVMM
        if ~silent
            fprintf('Calculating archetypes positions with SDVMM (Chan T-H, Chi C-Y, Huang Y-M and Ma W-K, 2009)\n');
        end
        
        for i=1:3*numIter
            [Archs,~] = SDVMM(DataPCA(:,1:NArchetypes-1)',NArchetypes,0);
            if ~isnan(Archs)
                %calculate the volume of the simplex
                Arch1Red=bsxfun(@minus,Archs,Archs(:,NArchetypes));
                VolArch(i)=abs(det(Arch1Red(:,1:end-1))/factorial(NArchetypes-1));
                %save the archetypes
                minArchsIter{i}=Archs;
            else
                VolArch(i)= NaN;
                minArchsIter{i}=NaN;
            end
        end
        [VolArchReal IndminVol]=max(VolArch); %find the minimal volume simplex
        ArchsMin=minArchsIter{IndminVol}; %get the minimal archetypes that form this simplex
        
        
    case 5 %    algNum=5 :> PCHA
        if ~silent
            fprintf('Calculating archetypes positions with PCHA (Morup M & Hansen LK, 2011)\n');
        end
        
        for i=1:3*numIter
            delta = 0;
            U=1:size(DataPCA,1); % Entries in X used that is modelled by the AA model
            I=1:size(DataPCA,1); % Entries in X used to define archetypes
            [Archs,~,~,~,varexpl]=PCHA1(DataPCA(:,1:min(NArchetypes-1,DataDim))',NArchetypes,I,U,delta);
            if ~isnan(Archs)
                %calculate the volume of the simplex
                Arch1Red=bsxfun(@minus,Archs,Archs(:,NArchetypes));
                %This formula is correct for a simplex
                VolArch(i)=abs(det(Arch1Red(:,1:end-1))/factorial(NArchetypes-1));
                %save the archetypes
                minArchsIter{i}=Archs(1:end,1:NArchetypes);
            else
                VolArch(i)= NaN;
                minArchsIter{i}=NaN;
            end
        end
        [VolArchReal IndminVol]=max(VolArch); %find the minimal volume simplex
        ArchsMin=minArchsIter{IndminVol}; %get the minimal archetypes that form this simplex
end

% We already report the following in 'findArchetypes.m', so we don't need it here:
% if ~silent
%    disp('finished finding the archetypes');
% end
