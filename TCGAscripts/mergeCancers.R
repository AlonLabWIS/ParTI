rm(list=ls())

cancerIDs <- read.table("../TCGA_frac_nArchs.tab", h=F, as.is=T)[,1]
cancerIDs <- setdiff(cancerIDs, "SYNT")

cancerID <- cancerIDs[1]

##################################################
## Merge gene expression

## pb <- txtProgressBar(style=3);
## allExpr <- 
##     sapply(cancerIDs, function(cancerID) {
##         expr <- read.table(sprintf("../%s_UCSC/expMatrix.tsv", cancerID),
##                            h=F, as.is=T, row.names=1)
##         expr0 <- apply(expr, 2, function(x) { x - mean(x) } )
##         setTxtProgressBar(pb, which(cancerID == cancerIDs) / length(cancerIDs))
##         return(expr0);
##     })
## save(allExpr, file="allExpr.rda");
load("allExpr.rda");

allExprMat <-
    matrix(NA, sum(sapply(allExpr, function(x) { nrow(x) })), ncol(allExpr[[1]]))
n <- 1;
for (i in 1:length(allExpr)) {
    allExprMat[seq(n, n+nrow(allExpr[[i]])-1),] <- allExpr[[i]]
    n <- n + nrow(allExpr[[i]]);
}

## round(allExprMat[1:5,1:5], 4)
write.table(round(allExprMat, 4), file="expMatrix.csv", quote=F,
            row.names=F, col.names=F, sep=",")

rm(allExpr); gc();
decomp <- svd(allExprMat, nu=3, nv=3)


barplot((decomp$d^2 / sum(decomp$d^2))[1:10], main="% variance explained")

dim(decomp$v)
exprProj <- allExprMat %*% decomp$v

library(LSD); heatscatter(exprProj[,1], exprProj[,2], add.cont=T)
library(rgl); plot3d(exprProj[,1], exprProj[,2], exprProj[,3])

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

allMutMat <- matrix(NA, 0, length(mutsCount));
colnames(allMutMat) <- names(mutsCount)

cancerID <- cancerIDs[1]
cancerID <- "PRAD"

for ( cancerID in cancerIDs ) {
    cat(sprintf("Importing mutations from %s\n", cancerID))
    mutations <-
        read.table(sprintf("../%s_UCSC/mutMatrix_reOrdered_booleanized_justData.csv",
                           cancerID), h=F, as.is=T, sep=",")
    colnames(mutations) <-
        read.table(sprintf("../%s_UCSC/mutMatrix_reOrdered_booleanized_geneNames.list",
                           cancerID), h=F, as.is=T, sep=",")[,1]

    mutationsPadded <- matrix(NA, nrow(mutations), length(mutsCount));
    colnames(mutationsPadded) <- names(mutsCount);
    mutationsPadded <- as.data.frame(mutationsPadded);

    toKeepMuts <- intersect(colnames(mutations), names(mutsCount))
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
write.table(colnames(allMutMat),
            file="mutMatrix_reOrdered_booleanized_geneNames.list",
            quote=F, sep=",", row.names=F, col.names=F);

