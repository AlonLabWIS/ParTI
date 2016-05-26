rm(list=ls())

myFiles <-
    c("LUAD_UCSC/MSigDBenrichment_continuous_significant.csv",
      ## "SKCM_UCSC/MSigDBenrichment_continuous_significant.csv",
      ## "UCEC_UCSC/MSigDBenrichment_continuous_significant.csv",
      "BRCA_UCSC/MSigDBenrichment_continuous_significant.csv",
      ## "KIRC_UCSC/MSigDBenrichment_continuous_significant.csv",
      ## "LUSC_UCSC/MSigDBenrichment_continuous_significant.csv",
      ## "PRAD_UCSC/MSigDBenrichment_continuous_significant.csv",
      "COAD_UCSC/MSigDBenrichment_continuous_significant.csv",
      "GBM_UCSC/MSigDBenrichment_continuous_significant.csv",
      "LGG_UCSC/MSigDBenrichment_continuous_significant.csv",
      ## "THCA_UCSC/MSigDBenrichment_continuous_significant.csv",
      "BLCA_UCSC/MSigDBenrichment_continuous_significant.csv",
      "LIHC_UCSC/MSigDBenrichment_continuous_significant.csv");
      

labs <- gsub("_UCSC.*$", "", myFiles)

myFile <- myFiles[1]
enrichRes <- sapply(myFiles, function(myFile) {
    cancerType <- gsub("_UCSC.*$", "", myFile)
    myTab <- read.csv(myFile, as.is=T)
    ## myTab <- data.frame(cancerType, myTab[myTab[,3] < 1e-12,1:2])
    myTab <- data.frame(cancerType, myTab[myTab[,3] < 1,1:2])
    return(myTab);
})
colnames(enrichRes) <- labs;

allMSigDBcat <- unlist(sapply(labs, function(lab) { enrichRes[["Feature.Name",lab]] }))
allMSigDBcat <- rev(sort(table(allMSigDBcat)))
plot(allMSigDBcat / length(labs),
     xlab="# categories", ylab="occurences per dataset", type="l",
     log="x")
abline(v=1000, lty=2)
plot(100 * cumsum(allMSigDBcat) / sum(allMSigDBcat),
     xlab="# categories", ylab="% enrichment covered", type="l")
allArchNs <-
    sapply(colnames(enrichRes),
                  function(can) {
                      unique(unlist(enrichRes[2,can])) })

myCats <- names(allMSigDBcat[1:1000]);
MSumMat <- matrix(NA, length(myCats), sum(sapply(allArchNs, length)))
rownames(MSumMat) <- myCats
colnames(MSumMat) <- 
    unlist(sapply(names(allArchNs),
                  function(can) {
                      paste(can, allArchNs[[can]], sep=".")  }))


i <- 1;
for (i in 1:length(labs)) {
    lab <- labs[i]
    j <- 1;
    for (j in allArchNs[[lab]]) {
        sel <- enrichRes[["archetype..",lab]] == j;
        offset <- 0;
        if ( i > 1 ) {
            offset <-
                sum(
                    sapply(seq(1, which(labs == labs[i-1])),
                           function(k) { length(allArchNs[[k]]) })
                    )
        }
        MSumMat[,offset + j] <-
            as.numeric(
                sapply(myCats, function(x) {
                    sum(enrichRes[["Feature.Name", lab]][sel] == x) > 0
                })
                )
    }
}

write.csv(MSumMat, file="MSumMat.csv")
          
## heatmap(MSumMat, scale="none", margins=c(5, 30))
dim(unique(MSumMat))
MSumMat0 <- MSumMat;
rownames(MSumMat0) <- c();
pdf("MSumMat.pdf", height=15, width=10)
## heatmap(unique(MSumMat0), scale="none", margins=c(1, 30))
h <- heatmap(MSumMat0, scale="none", margins=c(5, 5), keep.dendro=T)
dev.off();

myClusters <- list(
    h$Rowv[[2]][[2]][[2]],
    h$Rowv[[2]][[2]][[1]],
    h$Rowv[[2]][[1]],
    h$Rowv[[1]][[1]],
    h$Rowv[[1]][[2]]);

getLeaves <- function(y) {
    unlist(dendrapply(y, function(x) {
        if ( length(attr(x, 'label')) ) {
            return(attr(x, 'label'));
        } }))
}

## A: misc, among which cell adhesion
rownames(MSumMat)[getLeaves(myClusters[[1]])] 
## B: synapses & neurons, differentiation, MTOR signaling
rownames(MSumMat)[getLeaves(myClusters[[2]])]
## C: immune response, regulation of apoptosis, transcription, DNA
## damage checkpoint, degradation of ECM
rownames(MSumMat)[getLeaves(myClusters[[3]])]
## D1: gene expression (transcription & translation), DNA repair,
## mitosis
rownames(MSumMat)[getLeaves(myClusters[[4]])]
## D2: lipoprotein biosynthesis, respiration, regulation of mitosis,
## ribosome & translation, tRNAs, WNT signaling, DNA synthesis
rownames(MSumMat)[getLeaves(myClusters[[5]])]



heatmap(matrix(c(0,0,1), nrow=2), scale="none")
library(ade4)
## dudi1 <- dudi.coa(unique(MSumMat));
dudi1 <- dudi.coa(MSumMat);
plot(100 * cumsum(dudi1$eig) / sum(dudi1$eig),
     xlim=c(1,10), ylim=c(0,100), type="l",
     xlab="# components", ylab="% variance explained")
pdf("MSumMatCOA.pdf", height=4, width=6);
s.label(dudi1$co);
dev.off();
dudi1$co[order(dudi1$co[,1]),]
dudi1$li[order(dudi1$li[,1]),]
## Axis 1: mitosis & gene expression (-) vs differentiation (+)
dudi1$li[order(dudi1$li[,2]),2]
## Axis 2: TGF beta signaling (-) vs protein synthesis & metabolism (+)
dudi1$li[order(dudi1$li[,3]),]
## Axis 3: Pol II (-) vs Pol I (+); mainly there to explain BRCA.1 (LumB)
s.label(dudi1$co[,c(1,3)])

## Focus on the genes with best projetion
catLen <- apply(dudi1$li, 1, function(x) { sqrt(sum(x^2)) });
plot(ecdf(catLen))
cutOff <- 2.3;
abline(v=cutOff, lty=2)
ecdf(catLen)(cutOff)
catCoords <- dudi1$li[catLen>cutOff,];
library(rgl)
rownames(catCoords) <- sub("(KEGG|REACTOME|PID|BIOCARTA) ", "", rownames(catCoords))
plot3d(catCoords)
lims <- range(catCoords);
decorate3d(xlim=lims, ylim=lims, zlim=lims, aspect=T)
s.arrow(catCoords[,1:2])
rownames(catCoords)[catCoords[,2]>0]
rownames(catCoords)[catCoords[,2]<0]

x11();
plot3d(dudi1$li)
lims <- range(dudi1$co);
## decorate3d(xlim=lims, ylim=lims, zlim=lims, aspect=T)
ax3 <- dudi1$li[order(dudi1$li[,3]),]
ax3 <- ax3[ax3[,1]>0,]
ax3 <- ax3[abs(ax3[,2])<.1,]
write.csv(ax3, file="PCloadings.csv")

x11();
plot3d(dudi1$co, type="n")
lims <- range(dudi1$co);
## decorate3d(xlim=lims, ylim=lims, zlim=lims, aspect=T)
text3d(dudi1$co, texts=rownames(dudi1$co))

toHighlight <- c("MATRIX", "DOUBLE STRAND BREAK REPAIR",
                 "TRANSCRIPTION INITIATION FROM RNA POLYMERASE II PROMOTER",
                 "REACTOME TRANSLATION",
                 "RESPIRATORY ELECTRON TRANSPORT$",
                 "REACTOME REGULATION OF APOPTOSIS",
                 "REGULATION OF IMMUNE RESPONSE",
                 "PROGRAMMED CELL DEATH");
dudi1$li[grep("VEGF", rownames(dudi1$li)),]
dudi1$li[grep("REGULATION OF ANGIOGENESIS", rownames(dudi1$li)),]
dudi1$li[grep("RESPIRAT", rownames(dudi1$li)),]
dudi1$li[grep("TRANSLATION", rownames(dudi1$li)),]
dudi1$li[grep("TGF", rownames(dudi1$li)),]
dudi1$li[grep("MATRIX", rownames(dudi1$li)),]
dudi1$li[grep("IMMUN", rownames(dudi1$li)),]
## dudi1$li[grep("TCA CYCLE", rownames(dudi1$li)),]
dudi1$li[grep("REACTOME CYTOKINE SIGNALING IN IMMUNE SYSTEM", rownames(dudi1$li)),]
dudi1$li[grep("PROGRAMMED CELL DEATH", rownames(dudi1$li)),]


## Next, let's project chosen functions on the 3D PCA of samples
myIdcs <- unlist(sapply(toHighlight, function(x) { grep(x, rownames(dudi1$li)) }))

plot3d(dudi1$li, type="n")
lims <- range(dudi1$li);
## decorate3d(xlim=lims, ylim=lims, zlim=lims, aspect=T)
text3d(dudi1$li[myIdcs,1:3], texts=rownames(dudi1$li)[myIdcs])

pdf("MCAaxes.pdf", height=15, width=30);
par(mfrow=c(1,2))
s.arrow(dudi1$li[myIdcs,1:3])
text(dudi1$co[,1], dudi1$co[,2], labels=rownames(dudi1$co))
s.arrow(dudi1$li[myIdcs,1:3], xax=1, yax=3)
text(dudi1$co[,1], dudi1$co[,3], labels=rownames(dudi1$co))
dev.off();

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
