function [GOmatrix,GOfullNames,nGenesPerGO,GOcat2Genes]=MakeGOMatrix(DataPoints,GeneNames,MSigDBFiles,minGenes)
% This function transforms a matrix of samples x genes matrix to an
% samples x MSigDB matrix, for use in continuous enrichment analysis.
% In the resulting matrix, each column represent an MSigDB category, and
% expression values represent the average expression of genes belonging to
% that MSigDB category.
% 
% Arguments
%  DataPoints: samples x genes
%  GeneNames: a list of IDs that are known in MSigDB, in the same order as
%             the columns of DataPoints
%  MSigDBFiles: a list of files mapping MSigDB categories to genes names, in MSigDB format
%  minGenes: only consider MSigDB categories with at least minGenes genes
%
% Output
%  GOmatrix: transformed matrix, samples x GO categories
%  GOfullNames: the names of the MSigDB categories (= columns of GOmatrix)
%  nGenesPerGO: how many genes we averaged on (= columns of GOmatrix)
%  GOcat2Genes: a genes x MSigDB categories matrix indicating, for each MSigDB 
%               category which genes are part of this go category

%% Initialize
% load('Data/Cancer/CancerGenesNames.mat');
% load('Data/Cancer/CancerDataPoints.mat');
% GeneNames = Top_mRNA_genes;
% clear Top_mRNA_genes;

%% Function starts here
%MSigDBFiles = {'MSigDB/c2.cp.v4.0.symbols.gmt', 'MSigDB/c5.all.v4.0.symbols.gmt'};
% minGenes = 5;

%% Load MSigDB files
GOtoGenes = cell(1, 10000);
GOfullNames = cell(1, 10000);
fileIdx = 1;
nGOs = 0;
for fileIdx = 1:length(MSigDBFiles)
    DBFile = MSigDBFiles{fileIdx};
    GOfileId = fopen(DBFile, 'r');
    while ~feof(GOfileId)
        line = fgetl(GOfileId);
        nFields = length(find(line == char(9)))+1;
        myCols = cell(1,nFields);
        %Split line into tokens
        remain = line;
        for fieldIdx = 1:nFields
            [str, remain] = strtok(remain, char(9));
            %if isempty(str), break; end
            myCols{fieldIdx} = str;
        end
        if length(myCols)-2 < minGenes
            continue;
        end
        nGOs = nGOs + 1;
        GOtoGenes{nGOs} = myCols(3:end);
        GOfullNames{nGOs} = myCols{1};
    end
    fclose(GOfileId);
end
GOtoGenes = GOtoGenes(1:nGOs);
GOfullNames = GOfullNames(1:nGOs);

%% Now compute GO matrix from gene matrix
DataPointsDims = size(DataPoints);
GOmatrix = zeros(DataPointsDims(1), nGOs);
GOcat2Genes=zeros(DataPointsDims(2), nGOs);
GOmatrix(:,:) = NaN;
nGenesPerGO = zeros(1, nGOs);

GOidx = 1;
for GOidx = 1:nGOs
    myGenes = GOtoGenes{GOidx};    
    matchingIdcs = zeros(1, length(myGenes));
    matchingIdcs(1:end) = NaN;
    
    geneIdx = 1;
    for geneIdx = 1:length(myGenes)
        matchVec = strcmp(myGenes{geneIdx}, GeneNames);
        if ( sum(matchVec) == 0 ) %this gene is not in the dataset
            continue;
        end
        dataGeneIdx = find(matchVec);
        if length(dataGeneIdx) > 1
            error(sprintf('Could not compute GOmatrix because gene named %s appears to occure more than once in the data.', GeneNames{dataGeneIdx(1)}));
        end
        matchingIdcs(geneIdx) = dataGeneIdx;
    end
    matchingIdcs = matchingIdcs(~isnan(matchingIdcs));
    GOcat2Genes(matchingIdcs,GOidx)=1;
    if isempty(matchingIdcs) %this GO category matches no gene in the dataset
        continue;
    end
    if length(matchingIdcs) < minGenes
        continue;
    end
    %Average gene expression
    GOmatrix(:,GOidx) = mean(DataPoints(:,matchingIdcs), 2);
    nGenesPerGO(GOidx) = length(matchingIdcs);
end

%% Post-processing
keepGOs = sum(isnan(GOmatrix),1) == 0;

GOmatrix = GOmatrix(:,keepGOs);
GOfullNames = GOfullNames(keepGOs);
nGenesPerGO = nGenesPerGO(keepGOs);
GOcat2Genes=GOcat2Genes(:,keepGOs);

if isempty(GOfullNames)
    fprintf('No gene category of the database contained at least %d genes:\nTry lowering the minGenes parameter.\n', ...
            minGenes);
end