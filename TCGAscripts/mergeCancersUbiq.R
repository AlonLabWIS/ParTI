rm(list=ls())

library(tidyverse)

cancerIDs <- read.table("../TCGA_frac_nArchs.tab", h=F, as.is=T)[,1]
cancerIDs <- setdiff(cancerIDs, c("SYNT", "STAD"))
#STAD has a different gene set

## Merge gene expression data from different cancers, focusing on
## ubiquitous genes

cancerID <- cancerIDs[1]

##################################################
## Merge gene expression

## genePctiles <- 
##     lapply(cancerIDs, function(cancerID) {
##         cat(sprintf("Processing %s\n", cancerID))
##         expr <- read_tsv(sprintf("../%s_UCSC/expMatrix.tsv", cancerID),
##                          col_names=F)
##         expr <- as.data.frame(expr)
##         rownames(expr) <- expr[,1]
##         expr <- expr[,-1]

##         genes <- read_tsv(sprintf("../%s_UCSC/geneListExp.list", cancerID),
##                           col_names=F)
##         colnames(expr) <- unlist(genes[,1])
        
##         clinProps <-
##             read_tsv(sprintf("../%s_UCSC/discreteClinicalData_reOrdered.tsv",
##                              cancerID))
##         isNormal <- clinProps[,"sample_type"] == "Solid Tissue Normal"
##         isPrimary <- clinProps[,"sample_type"] == "Primary Tumor"
##         pctiles <- apply(expr[isPrimary,], 1, rank) / ncol(expr)
##         medPctiles <- apply(pctiles, 1, median)
##         names(medPctiles) <- colnames(expr)
        
##         return(enframe(medPctiles))
##     })
## save(genePctiles, file="genePctiles.rda")
load("genePctiles.rda")
i <- 2
allGenePctiles <- genePctiles[[1]]
for (i in 2:length(genePctiles)) {
    allGenePctiles <- inner_join(allGenePctiles, genePctiles[[i]], by="name")
}
## sum((allGenePctiles %>%
##     group_by(name) %>% summarize_all(mean) %>% .[,"name"] %>% unlist
##     %>% table)>1)
sum((allGenePctiles %>% .[,"name"] %>% unlist %>% table)>1) == 0 #are names unique?

allGenePctiles <- as.data.frame(allGenePctiles)
rownames(allGenePctiles) <- allGenePctiles[,1]
allGenePctiles <- allGenePctiles[,-1]
genePctilesAvg <- apply(allGenePctiles, 1, median)

## Write out percentiles across cancer types; gene are ordered as
## in expression data
write.table(genePctilesAvg, file="genePctilesAvg.tab", quote=F,
            row.names=F, col.names=F)

##################################################

load("allExpr.rda")
allExpr <- allExpr[cancerIDs]
## Gene names are already aligned: md5sum *_UCSC/geneListExp.list
allExprMat <-
    matrix(NA,
           sum(sapply(allExpr, function(x) { nrow(x) })),
           ncol(allExpr[[1]]))
n <- 1;
for (i in 1:length(allExpr)) {
    allExprMat[seq(n, n+nrow(allExpr[[i]])-1),] <- allExpr[[i]]
    n <- n + nrow(allExpr[[i]]);
}

## round(allExprMat[1:5,1:5], 4)
## write.table(round(allExprMat, 4), file="expMatrix.csv", quote=F,
##             row.names=F, col.names=F, sep=",")
write_csv(as_tibble(round(allExprMat, 4)), "expMatrix.csv", col_names=F)

rm(allExpr); gc();
dim(allExprMat)

cancerID <- cancerIDs[1]
cancerLabel <- ## only for primary tumors
    unlist(sapply(cancerIDs, function(cancerID) {
        clinProps <-
            read_tsv(sprintf("../%s_UCSC/discreteClinicalData_reOrdered.tsv",
                             cancerID))
        rep(cancerID, sum(clinProps[,"sample_type"] == "Primary Tumor"))
    }))

isPT <- ## for all tumors
    unlist(sapply(cancerIDs, function(cancerID) {
        clinProps <-
            read_tsv(sprintf("../%s_UCSC/discreteClinicalData_reOrdered.tsv", cancerID))
        clinProps[,"sample_type"] == "Primary Tumor"
    }))

avgExprPT <- apply(allExprMat[isPT,], 2, function(x) { mean(x) })
length(apply(allExprMat[isPT,], 1, function(x) { 1 }))
## only for primary tumors:
allExprMat0 <- t(apply(allExprMat[isPT,], 1, function(x) { x - avgExprPT }))
cl <- "BRCA"
SDbyCancer <- sapply(unique(cancerLabel), function(cl) {
    ## sd(allExprMat0[cancerLabel == cl,])
    dcmp <- svd(allExprMat0[cancerLabel == cl,], nu=3, nv=3)
    sqrt(sum(dcmp$d[1:4]^2))
})
allExprMat0 <- 
    t(sapply(1:nrow(allExprMat0), function(i) {
        allExprMat0[i,]
        ## 1000 * allExprMat0[i,] / SDbyCancer[cancerLabel[i]]
    }))
dim(allExprMat0)

boxplot(apply(allExprMat0, 1, sd) ~ cancerLabel)

decomp <- svd(allExprMat0, nu=3, nv=3)

barplot((decomp$d^2 / sum(decomp$d^2))[1:10], main="% variance explained")

dim(decomp$v)
exprProj <- allExprMat0 %*% decomp$v
library(LSD); heatscatter(exprProj[,1], exprProj[,2], add.cont=T)

plot(exprProj[,1], exprProj[,2], col=as.factor(cancerLabel))
legend("bottomright", levels(as.factor(cancerLabel)),
       col=1:length(levels(as.factor(cancerLabel))), pch=20)
plot(exprProj[,2], exprProj[,3], col=as.factor(cancerLabel))
legend("bottomright", levels(as.factor(cancerLabel)),
       col=1:length(levels(as.factor(cancerLabel))), pch=20)
library(rgl); plot3d(exprProj[,1], exprProj[,2], exprProj[,3])
plot3d(exprProj[,1], exprProj[,2], exprProj[,3],
       col=as.numeric(as.factor(cancerLabel)))

cancerID <- "LGG"
cancerID <- "LUSC"
cancerID <- "LIHC"
cancerID <- "BLCA"
cancerID <- "BRCA"
cancerID <- "COAD"

EVs <-
    sapply(cancerIDs, function(cancerID) {
        sel <- cancerLabel == cancerID;
        plot3d(exprProj[sel,1], exprProj[sel,2], exprProj[sel,3], alpha=.8)
        points3d(exprProj[!sel,1], exprProj[!sel,2], exprProj[!sel,3],
                 col="grey", alpha=.5)
        EVs <- eigen(cov(exprProj[sel,1:3]))$values
        ## EVs <- cumsum(EVs)
        return(EVs)
    })
ggplot(as_tibble(t(EVs)) %>%
       mutate(cancer=cancerIDs) %>%
       gather(EVid, EV, -cancer)) +
    geom_line(aes(x=EVid, y=EV, col=cancer, group=cancer), size=1.5) +
    ylim(c(0,3200))


##################################################
## Merge treatments, M status and cancer types

## cancerID <- cancerIDs[1]
## cancerID <- "PRAD"

## allClinical <- 
##     lapply(cancerIDs, function(cancerID) {
##         clinical <-
##             read.table(sprintf("../%s_UCSC/discreteClinicalData_reOrdered.tsv",
##                                cancerID),
##                        h=T, as.is=T, sep="\t")
##         Midx <- grep("pathologic_M", colnames(clinical));
##         if ( length(Midx) == 0 ) {
##             Mvec <- rep(NA, nrow(clinical))
##         } else {
##             Mvec <- clinical[,Midx]
##         }
##         cType <- rep(cancerID, nrow(clinical));
##         clinical <-
##             data.frame(
##                 cancer_type=cType,
##                 pathologic_M=Mvec,
##                 clinical[,grep("^treat", colnames(clinical))])
##         return(clinical)
##     })

## sapply(allClinical, ncol)
## allClinMat <- allClinical[[1]]
## for (i in 2:length(allClinical)) {
##     allClinMat <- rbind(allClinMat, allClinical[[i]])
## }

## ## The number of patients with gene expression should match that with
## ## clinical record:
## nrow(allExprMat) == nrow(allClinMat)

## write.table(allClinMat, file="discreteClinicalData_reOrdered.tsv",
##             quote=F, sep="\t", row.names=F);

##################################################
## Merge all clinical features

featTypes <- c("discrete", "continuous");
featType <- featTypes[2]

for ( featType in featTypes ) {
    featsByCancer <- sapply(cancerIDs, function(cancerID) {
        clinical <-
            read.table(sprintf("../%s_UCSC/%sClinicalData_reOrdered.tsv",
                               cancerID, featType),
                       h=T, as.is=T, sep="\t")
        colnames(clinical)
    })
    featCount <- table(unlist(featsByCancer))

    hist(featCount, seq(0, max(featCount)));

    if ( featType == "discrete" ) {
        featCount <- c(featCount, "cancer_type"=0)
    }

    allClinMat <- matrix(NA, 0, length(featCount));
    colnames(allClinMat) <- names(featCount)

    cancerID <- cancerIDs[1]
    cancerID <- "PRAD"

    for ( cancerID in cancerIDs ) {
        clinical <-
            read.table(sprintf("../%s_UCSC/%sClinicalData_reOrdered.tsv",
                               cancerID, featType),
                       h=T, as.is=T, sep="\t")
        if ( featType == "discrete" ) {
            cType <- rep(cancerID, nrow(clinical));
            clinical <-
                data.frame(
                    cancer_type=cType,
                    clinical);
        }
        clinicalPadded <- matrix(NA, nrow(clinical), length(featCount));
        colnames(clinicalPadded) <- names(featCount);
        clinicalPadded <- as.data.frame(clinicalPadded);
        clinicalPadded[,colnames(clinical)] <- clinical

        allClinMat <- rbind(allClinMat, clinicalPadded)
    }

    ## For features that are cancer specific (one occurrence only), we
    ## tag the feature name with the cancer type
    toRename <- "psa_value";
    for ( toRename in names(featCount)[which(featCount == 1)] ) {
        doesFeatBelongTo <- sapply(featsByCancer, function(x) {
            any(x == toRename) })
        belongsTo <- names(doesFeatBelongTo)[doesFeatBelongTo]
        colnames(allClinMat)[colnames(allClinMat) == toRename] <-
            sprintf("%s.%s", belongsTo, toRename)
    }

    if ( featType == "discrete" ) { barplot(table(allClinMat[,"cancer_type"]), las=3)}
    ## The number of patients with gene expression should match that with
    ## clinical record:
    if ( nrow(allExprMat) != nrow(allClinMat) ) {
        stop("Mismatch between number of patients with gene expression and treatments")
    }
    write.table(allClinMat,
                file=sprintf("%sClinicalData_reOrdered.tsv", featType),
                quote=F, sep="\t", row.names=F);

    if ( featType == "continuous" ) {
        write.table(allClinMat,
                    file=sprintf("%sClinicalData_reOrdered_justData.csv", featType),
                    quote=F, row.names=F, col.names=F, sep=",", na="NaN");
        write.table(colnames(allClinMat),
                    file=sprintf("%sClinicalData_reOrdered_featNames.list", featType),
                    quote=F, row.names=F, col.names=F, sep=",");
    }
}

##################################################
## Merge patientIDs

patientIDs <- sapply(cancerIDs, function(cancerID) {
    patients <-
        read_tsv(sprintf("../%s_UCSC/patientIDs.list",
                         cancerID), col_names=F)
    as.character(unlist(patients))
})

if ( sum(sapply(patientIDs, length)) != nrow(allExprMat) ) {
    stop("number of patients does not match number of tumors")
}
write_tsv(as.data.frame(unlist(patientIDs)), "patientIDs.list", col_names=F)

##################################################
## Merge all mutations

mutsByCancer <- sapply(cancerIDs, function(cancerID) {
    mutations <-
        read.table(sprintf("../%s_UCSC/mutMatrix_reOrdered_booleanized_geneNames.list",
                           cancerID), h=F, as.is=T)
    mutations <- mutations[grep("=1$", mutations[,1]),1]
})
mutsCount <- table(unlist(mutsByCancer))

hist(mutsCount, seq(0, max(mutsCount)));

## Let's keep all mutations re-occuring it at least 5 cancer subtypes,
## which gives us 500 mutations to work with and see if ParTI eats
## it.
mutsCount <- mutsCount[mutsCount >= 7]
length(mutsCount)

## Let's add commonly mutated genes
library(tidyverse); library(stringr)
keptSNVs <- 
    unique(
        c(names(mutsCount),
          read_csv(str_c("~/work/cancerTaskAtlas/",
                         "SandersNatGenetics2013/",
                         "Ciriello2013_ng.2762-S2.csv")) %>% 
          filter(Type == "MUTATION") %>%
          transmute(mutation=str_c(`Altered Locus/Gene`, "=1")) %>% 
          as.data.frame() %>% .[,1]
          )
    )
length(keptSNVs)

allMutMat <- matrix(NA, 0, length(keptSNVs));
colnames(allMutMat) <- keptSNVs

cancerID <- cancerIDs[1]
cancerID <- "BLCA"

for ( cancerID in cancerIDs ) {
    cat(sprintf("Importing mutations from %s\n", cancerID))
    mutations <-
        read.table(sprintf("../%s_UCSC/mutMatrix_reOrdered_booleanized_justData.csv",
                           cancerID), h=F, as.is=T, sep=",")
    colnames(mutations) <-
        read.table(sprintf("../%s_UCSC/mutMatrix_reOrdered_booleanized_geneNames.list",
                           cancerID), h=F, as.is=T, sep=",")[,1]

    mutationsPadded <- matrix(NA, nrow(mutations), length(keptSNVs));
    colnames(mutationsPadded) <- keptSNVs;
    mutationsPadded <- as.data.frame(mutationsPadded);

    toKeepMuts <- intersect(colnames(mutations), keptSNVs)
    mutationsPadded[,toKeepMuts] <- mutations[,toKeepMuts]

    allMutMat <- rbind(allMutMat, mutationsPadded)
}

## The number of patients with gene expression should match that with
## mutation record:
if ( nrow(allExprMat) != nrow(allMutMat) ) {
    stop("Mismatch between number of patients with gene expression and treatments")
}
write.table(allMutMat,
            file="mutMatrix_reOrdered_booleanized_justData.csv",
            quote=F, sep=",", row.names=F, col.names=F, na="NaN");
write.table(colnames(allMutMat),        #
            file="mutMatrix_reOrdered_booleanized_geneNames.list",
            quote=F, sep=",", row.names=F, col.names=F);

##################################################
## Merge treatments

## WARNING this code seems flawed somewhere:
## TCGA-DK-AA6S-01 and TCGA-EL-A3ZS-11 after merging are wrong!

## treatPerCancer <-
##     lapply(cancerIDs, function(cancerID) {
##         discrete <-
##             read_tsv(
##                 sprintf("../%s_UCSC/treatments_reOrdered.tab",
##                         cancerID),
##                 col_names=c("patient",
##                             as.character(unlist(
##                                 read_tsv(
##                                     sprintf("../%s_UCSC/treatments_reOrdered.tab",
##                                             cancerID),
##                                     col_names=F, comment="TCGA")))),
##                 skip=1)
##     })

## write_tsv(bind_rows(treatPerCancer), "treatTab.tsv")

