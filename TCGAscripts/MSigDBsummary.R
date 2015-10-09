rm(list=ls())

myFiles <-
    c("BRCA_UCSC/MSigDBenrichment_continuous_significant.csv",
      ## "KIRC_UCSC/MSigDBenrichment_continuous_significant.csv",
      "LUAD_UCSC/MSigDBenrichment_continuous_significant.csv",
      ## "OV_UCSC/MSigDBenrichment_continuous_significant.csv",
      "SKCM_UCSC/MSigDBenrichment_continuous_significant.csv",
      "UCEC_UCSC/MSigDBenrichment_continuous_significant.csv")
labs <- gsub("_UCSC.*$", "", myFiles)

myFile <- myFiles[1]
enrichRes <- sapply(myFiles, function(myFile) {
    myLab <- gsub("_UCSC.*$", "", myFile)
    myTab <- read.csv(myFile, as.is=T)
    ## myTab <- data.frame(myLab, myTab[myTab[,3] < 1e-12,1:2])
    myTab <- data.frame(myLab, myTab[myTab[,3] < 1,1:2])
    return(myTab);
})
colnames(enrichRes) <- labs;

allMSigDBcat <- unlist(sapply(labs, function(lab) { enrichRes[["Feature.Name",lab]] }))
allMSigDBcat <- rev(sort(table(allMSigDBcat)))
plot(allMSigDBcat,
     xlab="# categories", ylab="% enrichment covered", type="l",
     log="x")
abline(v=100, lty=2)
plot(100 * cumsum(allMSigDBcat) / sum(allMSigDBcat),
     xlab="# categories", ylab="% enrichment covered", type="l")

myCats <- names(allMSigDBcat[1:400]);
MSumMat <- matrix(NA, length(myCats), 4*length(labs))
rownames(MSumMat) <- myCats
colnames(MSumMat) <- paste(rep(labs, each=4), rep(1:4, 4))
lab <- labs[1]
i <- 1;
for (i in 1:length(labs)) {
    lab <- labs[i]
    j <- 1;
    for (j in 1:4) {
        sel <- enrichRes[["archetype..",lab]] == j;
        MSumMat[,(i-1)*4 + j] <-
            as.numeric(
                sapply(myCats, function(x) {
                    sum(enrichRes[["Feature.Name", lab]][sel] == x) > 0
                })
                )
    }
}
## heatmap(MSumMat, scale="none", margins=c(5, 30))
dim(unique(MSumMat))
pdf("MSumMat.pdf", height=15, width=10)
heatmap(unique(MSumMat), scale="none", margins=c(1, 30))
dev.off();

library(ade4)
dudi1 <- dudi.coa(unique(MSumMat));
pdf("MSumMatCOA.pdf", height=4, width=6);
s.label(dudi1$co);
dev.off();
dudi1$co[order(dudi1$co[,1]),]
dudi1$li[order(dudi1$li[,1]),]
## Axis 1: mitosis & gene expression (-) vs differentiation (+)
dudi1$li[order(dudi1$li[,2]),]
## Axis 2: TGF beta signaling (-) vs protein synthesis & metabolism (+)
dudi1$li[order(dudi1$li[,3]),]
## Axis 3: Pol II (-) vs Pol I (+); mainly there to explain BRCA.1 (LumB)
s.label(dudi1$co[,c(1,3)])

## image(t(MSumMat - .5) %*% (MSumMat -.5) / nrow(MSumMat))
## heatmap(t(MSumMat - .5) %*% (MSumMat -.5) / nrow(MSumMat))

##################################################
## Let's try something else

myCats <- names(allMSigDBcat);
MSumMat <- matrix(NA, length(myCats), 4*length(labs))
rownames(MSumMat) <- myCats
colnames(MSumMat) <- paste(rep(labs, each=4), rep(1:4, 4))
lab <- labs[1]
i <- 1;
for (i in 1:length(labs)) {
    lab <- labs[i]
    j <- 1;
    for (j in 1:4) {
        sel <- enrichRes[["archetype..",lab]] == j;
        MSumMat[,(i-1)*4 + j] <-
            as.numeric(
                sapply(myCats, function(x) {
                    sum(enrichRes[["Feature.Name", lab]][sel] == x) > 0
                })
                )
    }
}

selCats <- c();

heatmap(MSumMat, scale="none", margins=c(5, 30))

library(ade4)
dudi1 <- dudi.coa(MSumMat);
s.label(dudi1$co)
dudi1$co[order(dudi1$co[,1]),]
dim(dudi1$li)
