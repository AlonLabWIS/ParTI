rm(list=ls())

library(fdrtool)
library(gdata)
library(gplots)
library(ade4)
library(rgl)

## TCGAfracNarchs <- read.table("TCGA_frac_nArchs.tab", as.is=T, row.names=1)
## colnames(TCGAfracNarchs) <- c("quantile", "nArchetypes");
## cancerIDs <- rownames(TCGAfracNarchs);

## cancerID <- cancerIDs[2]
## cancerID <- "GBM"
## treatPerCancer <-
##     sapply(cancerIDs, function(cancerID) {
##         discrete <-
##             read.table(
##                 sprintf("%s_UCSC/treatTab.tsv",
##                         cancerID), as.is=T, h=T, sep="\t")[,c(-1, -2)]
##         apply(discrete, 2, function(x) { mean(x, na.rm=T) })
##     })

## exprPerCancer <-
##     sapply(cancerIDs, function(cancerID) {
##         avgExpr <-
##             read.table(
##                 sprintf("%s_UCSC/avgExpr.tsv",
##                         cancerID), as.is=T, h=F, sep="\t")[,1]
##         geneList <-
##             read.table(
##                 sprintf("%s_UCSC/geneListExp.list",
##                         cancerID), as.is=T, h=F, sep="\t")[,1]
##         names(avgExpr) <- geneList;
##         return(avgExpr)
##     })

## ## exprCancers <- 
## ##     lapply(c(cancerIDs, "ALL"), function(cancerID) {
## ##         cat(paste("Loading", cancerID, "expression\n"))
## ##         geneList <-
## ##             read.table(
## ##                 sprintf("%s_UCSC/geneListExp.list",
## ##                         cancerID), as.is=T, h=F, sep="\t")[,1]
## ##         expr <- read.csv(sprintf("%s_UCSC/expMatrix.csv", cancerID),
## ##                          as.is=T, h=F)
## ##         colnames(expr) <- geneList;
## ##         return(expr)
## ##     })
## ## names(exprCancers) <- c(cancerIDs, "ALL")
## ## save(exprCancers, file="exprCancers.rda")
## load("exprCancers.rda")

## ## archPerCancerList <-
## ##     sapply(c(cancerIDs, "ALL"), function(cancerID) {
## ##         cat(paste("Loading", cancerID, "archetypes\n"))
## ##         arcsOrig <-
## ##             t(read.table(
## ##                 sprintf("%s_UCSC/arcsOrig_genes.tsv",
## ##                         cancerID), as.is=T, h=F))
## ##         geneList <-
## ##             read.table(
## ##                 sprintf("%s_UCSC/geneNamesAfterExprFiltering.list",
## ##                         cancerID), as.is=T, h=F, sep="\t")[,1]
## ##         rownames(arcsOrig) <- geneList;
## ##         return(arcsOrig)
## ##     })
## ## save(archPerCancerList, file="archPerCancerList.rda")
## load("archPerCancerList.rda")

## ## Keep only genes that are present in all the archetypes
## geneTab <- table(unlist(sapply(archPerCancerList, function(x) { rownames(x) })))
## ## These genes are found in all archetypes of all cancer
## ubiqGenes <- names(geneTab)[geneTab == (length(cancerIDs)+1)] 

## alnArchs <-
##     matrix(NA, length(ubiqGenes), 
##            sum(sapply(archPerCancerList, function(a) { ncol(a) })));
## rownames(alnArchs) <- ubiqGenes;
## n <- 0;
## for (i in 1:length(archPerCancerList)) {
##     ## Normalize out abs gene expression inside each cancer
##     ## archPerCancerList[[i]] <-
##     ##     t(apply(archPerCancerList[[i]], 1, function(x) {x - mean(x)}))
##     for (j in 1:ncol(archPerCancerList[[i]]) ) {
##         n <- n + 1;
##         alnArchs[,n] <- archPerCancerList[[i]][ubiqGenes,j]
##     }
## }

## a <- names(archPerCancerList)[1]
## colnames(alnArchs) <-
##     unlist(sapply(names(archPerCancerList), function(a) {
##         sprintf("%s.%d",
##                 rep(a, ncol(archPerCancerList[[a]])),
##                 1:ncol(archPerCancerList[[a]])
##                 )
##     }))

## ## Compute optimal cut-off on patent response
## respRanking <- c("Clinical Progressive Disease",
##                  "Stable Disease", "Partial Response",
##                  "Complete Response")

## pdf("respCutOffs.pdf", height=3*4, width=4*4)
## par(mfrow=c(3,4))
## cancerID <- "GBM"
## respCutOffs <-
##     sapply(cancerIDs, function(cancerID) {
##         cat(paste(cancerID,"\n"))
        
##         treatments <- 
##             read.table(
##                 sprintf("%s_UCSC/discreteClinicalData_reOrdered.tsv",
##                         cancerID), as.is=T, h=T, sep="\t")
##         treatments <- treatments[,c(which("treat.patientBestResp" == colnames(treatments)),
##                                     grep("treat.target.", colnames(treatments)))]
##         colnames(treatments) <- gsub("^treat.target.", "", colnames(treatments))
##         colnames(treatments)[1] <- "best.resp";

##         best.respVec <-
##             sapply(treatments[,"best.resp"], function(x) {
##                 if ( is.na(x) ) { return(NA) }
##                 which(respRanking == x)
##             })
##         treatments[,"best.resp"] <- best.respVec;

##         if ( ! any(!is.na(best.respVec)) ) { return(NA) }
        
##         bp <- barplot(table(as.numeric(best.respVec)), main=cancerID)
##         respCutOff <- median(best.respVec, na.rm=T);
##         respCutOff <- respCutOff +
##             which.min(abs(c(sum(best.respVec < respCutOff - .5, na.rm=T) / sum(!is.na(best.respVec)),
##                             sum(best.respVec < respCutOff + .5, na.rm=T) /
##                             sum(!is.na(best.respVec))) - .5)) - 1.5
##         ## if ( respCutOff == 1 ) {
##         ##     respCutOff <- 1.5
##         ## } else if ( respCutOff == 4 ) {
##         ##     respCutOff <- 3.5
##         ## }
##         abline(v=(respCutOff - 1)/(nrow(bp)-1) * diff(range(bp)) + min(bp), lty=2)
##         return(respCutOff)
##     })
## dev.off();

## save.image("drugVsCancers_init.rda")

## Done Init

load("drugVsCancers_init.rda")

##################################################

drugTargets <- read.csv("drugTargets.csv", as.is=T)
drugTargets[is.na(drugTargets)] <- 0
head(drugTargets)

treatPerCancerMat <-
    sapply(which(names(treatPerCancer) != "SYNT"), function(i) {
        treatPerCancer[[i]]
    })
colnames(treatPerCancerMat) <- setdiff(names(treatPerCancer), "SYNT")
treatPerCancer <- treatPerCancerMat;
rm(treatPerCancerMat)

dim(treatPerCancer)
which(apply(treatPerCancer, 2, sum) == 0)
treatPerCancer <- treatPerCancer[apply(treatPerCancer, 1, sum) > 0,]

heatmap.2(treatPerCancer, scale="none",
          col = heat.colors(256),
          hclustfun=function(x) { hclust(x, method="ward.D") })

treatPerCancer0 <-
    t(apply(treatPerCancer, 1, function(x) { (x - median(x)) / sd(x) }))
sdmd <- function(x) { mean(abs(x - median(x))) }
sum(apply(treatPerCancer0, 1, sd) < 1e-10)

heatmap.2(treatPerCancer0, scale="none",
          col = heat.colors(256),
          hclustfun=function(x) { hclust(x, method="ward.D") })

treatPerCancer0["dabrafenib",]

heatmap.2(log(treatPerCancer + 1), scale="none",
          col = heat.colors(256),
          hclustfun=function(x) { hclust(x, method="ward.D") })

treatPerCancer0["radiotherapy",]


hist(log10(apply(treatPerCancer, 1, mean)),
     xlab="log10 average usage frequency")
isFrequent <- apply(treatPerCancer, 1, mean) > .01
heatmap(treatPerCancer[isFrequent,], scale="none",
        hclustfun=function(x) { hclust(x, method="ward.D") })

heatmap.2(treatPerCancer[isFrequent,], scale="none",
          col = heat.colors(256),
          hclustfun=function(x) { hclust(x, method="ward.D") })

##################################################
## Same analysis by mechanism

d <- "radiotherapy"
allTherapies <- rownames(treatPerCancer);
allTherapies <- allTherapies[-grep("^target", allTherapies)]
alnDrugTargets <-
    sapply(allTherapies, function(d) {
        sel <-
            gsub("^([0-9])", "X\\1",
                 gsub("[ -]", ".", drugTargets[,1])) == d;
        if ( sum(sel) != 1 ) { cat(sprintf("Trouble with %s\n", d)) }
        return(as.numeric(drugTargets[sel,c(-1, -2)]))
    })

treatPerCancerS <-
    alnDrugTargets %*%
    treatPerCancer[-grep("^target", rownames(treatPerCancer)),]
rownames(treatPerCancerS) <- colnames(drugTargets)[c(-1,-2)]

heatmap.2(treatPerCancerS, scale="none",
          col = heat.colors(256), margins = c(5, 15),
          hclustfun=function(x) { hclust(x, method="ward.D") })

treatPerCancerS0 <-
    t(apply(treatPerCancerS, 1,
            function(x) { (x - median(x)) ## / sd(x)
                      }))

pdf("treatPerCancerS0.pdf", height=5, width=7.5);
heatmap.2(treatPerCancerS0, scale="none", trace="none",
          col = heat.colors, margins = c(5, 15),
          hclustfun=function(x) { hclust(x, method="ward.D") })
dev.off()

## Take out average gene expression to focus on changes between cancers
exprPerCancer0 <- t(apply(exprPerCancer, 1, function(x) { x - mean(x) }))

minExpr <- rev(sort(apply(exprPerCancer, 1, mean)))[100]
sel <- apply(exprPerCancer, 1, mean) >= minExpr
sel <- rev(order(apply(exprPerCancer[,setdiff(colnames(exprPerCancer), "SYNT")], 1, sd)))[1:1000]
topExprPerCancer <-
    t(apply(exprPerCancer[sel,setdiff(colnames(topExprPerCancer), "SYNT")],
          1, function(x) { x - mean(x) }))

pdf("exprClustering.pdf", height=5, width=5);
heatmap.2(topExprPerCancer, trace="none", scale="none")
dev.off();

##################################################

pdf("distTreatExpr.pdf", height=8, width=8)
par(mfrow=c(2,2))

par(mfrow=c(2,2))
dTreat <- as.matrix(dist(t(treatPerCancerS0)))
dTreatN <- dTreat / median(dTreat)
image(dTreatN, main="Treatments")

dExpr <- as.matrix(dist(t(topExprPerCancer)))
dExprN <- dExpr / median(dExpr)
image(dExprN, main="Expression")

## plot(dExprN, dTreatN)
## cor.test(dExprN, dTreatN)

rndTEPR <-
    t(apply(topExprPerCancer, 1,
            function(x) { sample(x) }))
dRndTEPR <- as.matrix(dist(t(rndTEPR)))
dRndTEPRN <- dRndTEPR / median(dExpr)
image(dRndTEPRN, main="Random expr.")

## Normalization factor should be the same as is Expression?

rndDists <-
    sapply(1:1000, function(i) {
        rndTEPR <-
            t(apply(topExprPerCancer, 1,
                    function(x) { sample(x) }))
        dRndTEPR <- as.matrix(dist(t(rndTEPR)))
        dRndTEPRN <- dRndTEPR / median(dExpr)
        sum((dTreatN - dRndTEPRN)^2)
    })

initDist <- sum((dTreatN - dExprN)^2)
hist(rndDists, 20,
     xlim=range(c(rndDists, initDist)),
     xlab="distance between matrices")
abline(v=initDist, lwd=2, col="red")

dev.off();

##################################################

pdf("covTreatExpr.pdf", height=8, width=8)
par(mfrow=c(2,2))

tmp <- t(treatPerCancerS0) %*% treatPerCancerS0;
tmp2 <- tmp / 
    sqrt(as.matrix(apply(treatPerCancerS0^2, 2, sum)) %*%
         t(as.matrix(apply(treatPerCancerS0^2, 2, sum))))
image(tmp2, col = heat.colors(256), breaks=seq(-1, 1, len=257),
      main="Treatments")

image(cor(topExprPerCancer), col = heat.colors(256),
      breaks=seq(-1, 1, len=257), main="Expression")

rndTEPR <-
    t(apply(topExprPerCancer, 1,
            function(x) { sample(x) }))
image(cor(rndTEPR), col = heat.colors(256),
      breaks=seq(-1, 1, len=257), main="Expression")

rndDists <-
    sapply(1:1000, function(i) {
        rndTEPR <-
            t(apply(topExprPerCancer, 1,
                    function(x) { sample(x) }))
        sum((tmp2 - cor(rndTEPR))^2)
    })

initDist <- sum((tmp2 - cor(topExprPerCancer))^2)
hist(rndDists, 20,
     xlim=range(c(rndDists, initDist)),
     xlab="distance between matrices")
abline(v=initDist, lwd=2, col="red")
dev.off();

par(mfrow=c(1,2))
plot(hclust(dist(t(exprPerCancer0)), method="ward.D"), main="Expression")
plot(hclust(dist(t(treatPerCancerS0)), method="ward.D"), main="Treatments")

##################################################

## Now study archetypes instead of cancers

## Take out average gene expression to focus on changes between cancers
sel <- rev(order(apply(alnArchs[,-grep("SYNT|ALL", colnames(alnArchs))], 1, sd)))[1:1000]
## sel <- 1:nrow(alnArchs)
topExprArchs <-
    t(apply(alnArchs[sel,-grep("SYNT|ALL", colnames(alnArchs))], 1,
            function(x) { x - mean(x) }))

pdf("topExprArchsClustering.pdf", height=7, width=7);
## bitmap("topExprArchsClustering.png", height=10, width=10);
heatmap.2(topExprArchs, scale = "none", trace="none", col=heat.colors)
dev.off();

image(as.matrix(dist(t(topExprArchs))))

## Cluster using all genes, just dendrogram

## Is the distance between cancer centroids large compared to the
## distance between archetypes and their centroid?
head(alnArchs[sel,])
head(exprPerCancer[ubiqGenes,][sel,])
## archCentr <- alnArchs[sel,-grep("SYNT", colnames(alnArchs))]
archCentr <-
    cbind(alnArchs[sel,-grep("SYNT|ALL", colnames(alnArchs))],
          exprPerCancer[ubiqGenes,setdiff(colnames(exprPerCancer),
                                          c("SYNT", "ALL"))][sel,])

dudi1 <- dudi.pca(t(archCentr), scale=F, scannf=F, nf=3)
barplot(dudi1$eig / sum(dudi1$eig))

colVec <-
    c(rep(20, ncol(alnArchs[,-grep("SYNT|ALL", colnames(alnArchs))])),
      rep(2, ncol(exprPerCancer[,setdiff(colnames(exprPerCancer), c("SYNT", "ALL"))])))
cTypeFac <- as.factor(gsub("\\.[0-9]$", "", colnames(archCentr)));
shapeVec <- as.numeric(cTypeFac)

pdf("archCentrPCA.pdf", height=6, width=6);
plot(dudi1$li[,1:2], pch=shapeVec, col=colVec)
legend("bottomright", c("archetype", "centroid"), fill=unique(colVec))
legend("bottomleft", levels(cTypeFac), col="grey", pch=min(shapeVec):max(shapeVec))
dev.off();

shapeVec3 <- rep(NA, length(shapeVec))
shapeVec3[shapeVec == 1] <- "p"
shapeVec3[shapeVec == 2] <- "s"
library(rgl)
plot3d(dudi1$li[shapeVec == 1,1:3],
       type="p", col=colVec[shapeVec == 1])
plot3d(dudi1$li[shapeVec == 2,1:3],
       type="s", add=T, col=colVec[shapeVec == 2])

##################################################

## Create a treatment profile for each archetype.
## We need:
## - archetypes
## - gene expression in different tumors
## - treatment of each tumor

dim(alnArchs)

cancerID <- cancerIDs[1]
cancerID <- "LGG";
cancerID <- cancerIDs[6]
EVsL <- list();
load("allCancers_EVs.rda")

## treatPerCancerArch <- list();
## for ( cancerID in setdiff(cancerIDs, "SYNT") ) {
##     cat(paste(cancerID,"\n"))
##     geneList <-
##         read.table(
##             sprintf("%s_UCSC/geneListExp.list",
##                     cancerID), as.is=T, h=F, sep="\t")[,1]
##     expr <- read.csv(sprintf("%s_UCSC/expMatrix.csv", cancerID),
##                      as.is=T, h=F)
##     colnames(expr) <- geneList;
    
##     ## exprU <- t(expr[,ubiqGenes])
##     minExpr <- quantile(apply(expr, 2, mean), TCGAfracNarchs[cancerID,"quantile"]);
##     exprU <- t(expr[,apply(expr, 2, mean) > minExpr])
    
##     ## discreteClinicalData_reOrdered_withTreatment.tsv, select
##     ## only treat.target
##     treatments <- 
##         read.table(
##             sprintf("%s_UCSC/discreteClinicalData_reOrdered_withTreatment.tsv",
##                     cancerID), as.is=T, h=T, sep="\t")
##     treatments <- treatments[,grep("treat.target.", colnames(treatments))]
##     colnames(treatments) <- gsub("^treat.target.", "", colnames(treatments))

##     myArchs <- archPerCancerList[[cancerID]]
##     if ( sum(rownames(myArchs) != rownames(exprU)) > 0 ) {
##         stop("genes in archetypes and samples are mis-aligned.\n")
##     }
     
##     centerVec <- apply(exprU, 1, mean)
##     exprU0 <-
##         t(sapply(1:nrow(exprU), function(i) {
##             exprU[i,] - centerVec[i] }))
##     archs0 <-
##         t(sapply(1:nrow(exprU), function(i) {
##             myArchs[i,] - centerVec[i]
##         }))

##     nPCs <- 4;
##     if ( is.null(EVsL[[cancerID]]) ) {
##         ## EVsL[[cancerID]] <- eigen(cov(t(exprU0)))$vectors[,1:nPCs]
##         EVsL[[cancerID]] <- svd(exprU0, nu=nPCs, nv=nPCs)$u[,1:nPCs];
##         save(EVsL, file="allCancers_EVs.rda")
##     } else {
##         cat("Reusing previously computed eigenvectors.\n")
##     }
##     EVs <- EVsL[[cancerID]]

##     exprProj <- t(t(EVs) %*% exprU0)
##     archProj <- t(t(EVs) %*% archs0)
    
##     plot(exprProj[,1], exprProj[,2],
##          xlim=range(c(exprProj[,1], archProj[,1])),
##          ylim=range(c(exprProj[,2], archProj[,2])))
##     points(archProj, pch=20, col="red", cex=2)

##     ## Take the topPct closest to each archetype: use only half
##     ## the data
##     topPct <- 1 / (ncol(archs0)+1)
    
##     j <- 1;
##     archTreatProfile <-
##         sapply(1:nrow(archProj), function(j) {
##             dists <- sapply(1:nrow(exprProj), function(i) {
##                 sqrt(sum((exprProj[i,] - archProj[j,])^2))
##             })
##             closestPts <-
##                 order(dists)[1:round(topPct*length(dists))];
##             apply(treatments[closestPts,], 2, function(x) {
##                 mean(x, na.rm=T) })
##         })

##     treatPerCancerArch[[cancerID]] <- archTreatProfile
## }
## save(treatPerCancerArch, file="treatPerCancerArch.rda")

load("treatPerCancerArch.rda")
alnTreats <-
    matrix(NA, nrow(treatPerCancerArch[[1]]),
           sum(sapply(archPerCancerList, function(a) { ncol(a) })) -
           ncol(archPerCancerList[["SYNT"]]) -
           ncol(archPerCancerList[["ALL"]])
           );
rownames(alnTreats) <- rownames(treatPerCancerArch[[1]]);

n <- 0;
for (i in 1:length(treatPerCancerArch)) {
    for (j in 1:ncol(treatPerCancerArch[[i]]) ) {
        n <- n + 1;
        alnTreats[,n] <- treatPerCancerArch[[i]][,j]
    }
}

colnames(alnTreats) <-
    unlist(sapply(setdiff(cancerIDs, "SYNT"), function(a) {
        idx <- which(cancerIDs == a)
        sprintf("%s.%d",
                rep(a, ncol(treatPerCancerArch[[idx]])),
                1:ncol(treatPerCancerArch[[idx]])
                )
    }))

alnTreats0 <- t(apply(alnTreats, 1, function(x) { x - mean(x, na.rm=T) }))
alnTreats0[is.na(alnTreats0)] <- 0;

pdf("treatPerArch.pdf", height=6, width=12);
heatmap.2(alnTreats0[apply(alnTreats, 1, mean) > .1,], scale="none", trace="none",
          col = heat.colors, margins = c(5, 20),
          hclustfun=function(x) { hclust(x, method="ward.D") })
dev.off()

#####

pdf("distTreatExprArch.pdf", height=8, width=8)
par(mfrow=c(2,2))

dTreat <- as.matrix(dist(t(alnTreats0)))
dTreatN <- dTreat / median(dTreat)
image(dTreatN, main="Treatments")

dExpr <- as.matrix(dist(t(topExprArchs)))
dExprN <- dExpr / median(dExpr)
image(dExprN, main="Expression")

## plot(dExprN, dTreatN)
## cor.test(dExprN, dTreatN)

rndTEPR <-
    t(apply(topExprArchs, 1, function(x) { sample(x) }))
dRndTEPR <- as.matrix(dist(t(rndTEPR)))
dRndTEPRN <- dRndTEPR / median(dExpr)
image(dRndTEPRN, main="Random expr.")

## Normalization factor should be the same as is Expression?

rndDists <-
    sapply(1:1000, function(i) {
        rndTEPR <-
            t(apply(topExprArchs, 1,
                    function(x) { sample(x) }))
        dRndTEPR <- as.matrix(dist(t(rndTEPR)))
        dRndTEPRN <- dRndTEPR / median(dExpr)
        sum((dTreatN - dRndTEPRN)^2)
    })

initDist <- sum((dTreatN - dExprN)^2)
hist(rndDists, 20,
     xlim=range(c(rndDists, initDist)),
     xlab="distance between matrices")
abline(v=initDist, lwd=2, col="red")

dev.off();

##################################################

## Is there a specific interaction between treatments and archetypes
## associated to outcome?

cancerID <- cancerIDs[1]
cancerID <- "LGG";
cancerID <- cancerIDs[6]
cancerID <- "BRCA" #not a good cancer to find new drug-arch
                   #combinations: therapy is already well tuned
cancerID <- "COAD"
cancerID <- "BLCA";
cancerID <- "SYNT";
EVsL <- list();
load("allCancers_EVs.rda")

showEverything <- F;
responsesL <- list();

## treatPerCancerArch <- list();

pdf("treatArchResponse.pdf", height=8, width=8);
for ( cancerID in names(respCutOffs[!is.na(respCutOffs)]) ) {
    cat(paste(cancerID,"\n"))
    par(mfrow=c(2,2))

    expr <- exprCancers[[cancerID]];
    minExpr <- quantile(apply(expr, 2, mean), TCGAfracNarchs[cancerID,"quantile"]);
    exprU <- t(expr[,apply(expr, 2, mean) > minExpr])
    
    ## discreteClinicalData_reOrdered_withTreatment.tsv, select
    ## only treat.target
    treatments <- 
        read.table(
            sprintf("%s_UCSC/discreteClinicalData_reOrdered.tsv",
                    cancerID), as.is=T, h=T, sep="\t")
    treatments <- treatments[,c(which("treat.patientBestResp" == colnames(treatments)),
                                grep("treat.target.", colnames(treatments)))]
    colnames(treatments) <- gsub("^treat.target.", "", colnames(treatments))
    colnames(treatments)[1] <- "best.resp";
    colnames(treatments)

    best.respVec <-
        sapply(treatments[,"best.resp"], function(x) {
            if ( is.na(x) ) { return(NA) }
            which(respRanking == x)
        })
    treatments[,"best.resp"] <- best.respVec;
    
    myArchs <- archPerCancerList[[cancerID]]
    if ( sum(rownames(myArchs) != rownames(exprU)) > 0 ) {
        stop("genes in archetypes and samples are mis-aligned.\n")
    }
    
    centerVec <- apply(exprU, 1, mean)
    exprU0 <-
        t(sapply(1:nrow(exprU), function(i) {
            exprU[i,] - centerVec[i] }))
    archs0 <-
        t(sapply(1:nrow(exprU), function(i) {
            myArchs[i,] - centerVec[i]
        }))

    nPCs <- 4;
    ## if ( is.null(EVsL[[cancerID]]) ) {
        ## EVsL[[cancerID]] <- eigen(cov(t(exprU0)), symmetric=T)$vectors[,1:nPCs]
        EVsL[[cancerID]] <- svd(exprU0, nu=nPCs, nv=nPCs)$u[,1:nPCs];
        save(EVsL, file="allCancers_EVs.rda")
    ## } else {
    ##     cat("Reusing previously computed eigenvectors.\n")
    ## }
    EVs <- EVsL[[cancerID]]

    exprProj <- t(t(EVs) %*% exprU0)
    archProj <- t(t(EVs) %*% archs0)
    
    plot(exprProj[,1], exprProj[,2],
         xlim=range(c(exprProj[,1], archProj[,1])),
         ylim=range(c(exprProj[,2], archProj[,2])),
         main=cancerID)
    points(archProj, pch=20, col="blue", cex=2)

    ## Take the topPct closest to each archetype: use only half
    ## the data
    ## topPct <- .5  / ncol(archs0)
    topPct <- 1 / (ncol(archs0)+1)

    j <- 1;
    closestPtss <-
        sapply(1:nrow(archProj), function(j) {
            dists <- sapply(1:nrow(exprProj), function(i) {
                sqrt(sum((exprProj[i,] - archProj[j,])^2))
            })
            closestPts <-
                order(dists)[1:round(topPct*length(dists))];
            ## apply(treatments[closestPts,], 2, function(x) {
            ##     mean(x, na.rm=T) })
            return(closestPts)
        })
    archFac <- rep(NA, nrow(exprProj))
    for (i in 1:ncol(closestPtss)) {
        ## FIXME the last archetypes snatched patients from other
        ## archetypes
        archFac[closestPtss[,i]] <- i;
    }
    archFac <- as.factor(archFac);

    minOcc <- 10;
    treatOcc <- sort(apply(treatments[,-1], 2, function(x) { sum(x, na.rm=T) }))
    keptTreats <- names(which(treatOcc > minOcc))
    if ( length(keptTreats) == 0 ) { next; }
    treatFact <- as.data.frame(treatments[,keptTreats])
    for ( i in 1:ncol(treatFact) ) {
        treatFact[,i] <- as.factor(treatFact[,i])
    }
    ## table(treatments[,"best.resp"])
    colnames(treatFact) <- keptTreats

    respCutOff <- respCutOffs[cancerID]
    
    preGlmData <-
        data.frame(
            ## best.resp=as.numeric(treatments[,"best.resp"] !=
            ##     "Clinical Progressive Disease"),
            best.resp=as.numeric(treatments[,"best.resp"] > respCutOff),
            treatFact,
            arch=archFac)
    noNA <- !apply(is.na(preGlmData[,-ncol(preGlmData)]), 1, any) ## & Mstatus == "M0"
    ## Previously: glmData <- glmData[apply(is.na(glmData[,-ncol(glmData)]), 1, sum) == 0,]
    glmData <- preGlmData[noNA,]

    summary(glmData)
    i <- "nucleotide.depletion";
    i <- "DNA.damaging";
    i <- "angiogenesis.signalling.inhibitor";
    i <- "microtubule";

    responses <- lapply(colnames(glmData)[2:(ncol(glmData)-1)], function(i) {
        cat(paste(i, "\n"))
        tgd <- glmData[glmData[,i] == 1, c("best.resp","arch")]

        logORSE <- matrix(NA, 2, nrow(archProj));
        rownames(logORSE) <- c("logOR", "SE");
        colnames(logORSE) <- 1:ncol(logORSE)
        
        if ( sum(dim(table(tgd)) < 2) > 0 ) { return(logORSE); }
        cTab <- table(tgd) #+1 #FIXME
        keptArch <- which(apply(cTab, 2, sum) > 0);
        cTab <- cTab[,keptArch]
        ## p <- chisq.test(cTab, simulate.p.value=F)$p.value
        
        ## if ( p < .1 ) {
        if ( sum(cTab) >= minOcc ) {
            cTabP <- cTab + 1;
            logORSE <-
                sapply(1:ncol(cTabP), function(j) {
                    ORs <- c(cTabP[2,j] / cTabP[1,j],
                             (sum(tgd[,"best.resp"] == 1) + 1) /
                             (sum(tgd[,"best.resp"] == 0) + 1));
                    ## retMat <- matrix(NA, 2, ncol(cTab));
                    ## colnames(retMat) <- colnames(cTab);
                    ## rownames(retMat) <- c("logOR", "SE")
                    ## retMat["logOR",] <- log(ORs[1] / ORs[2]);
                    ## retMat["SE",] <-
                    ##     sqrt(sum(1/cTabP[,j]) +
                    ##          sum(1/c(sum(tgd[,"best.resp"] == 1) + 1,
                    ##                  sum(tgd[,"best.resp"] == 0) + 1)))
                    return(
                        c(logOR=log(ORs[1] / ORs[2]),
                          SE=sqrt(sum(1/cTabP[,j]) +
                              sum(1/c(sum(tgd[,"best.resp"] == 1) + 1,
                                      sum(tgd[,"best.resp"] == 0) +
                                      1))))
                        )
                    ## sqrt(sum(1/cTabP[,j]) + sum(1/apply(cTabP[,-j],
                    ## 1, sum))))
                })
            colnames(logORSE) <- keptArch;
            
            if ( any(abs(logORSE["logOR",]) -
                     1.96 * logORSE["SE",] > 0) || showEverything ) { #FIXME
                bp <- barplot(logORSE["logOR",],
                              names=paste("A", colnames(logORSE), sep=""),
                              main=sprintf("%s", i),
                              ylab=sprintf("log OR (n=%d)", nrow(tgd)))
                sapply(1:ncol(logORSE), function(j) {
                    arrows(bp[j,1], logORSE["logOR",j] +
                           1.96 * logORSE["SE",j],
                           bp[j,1], logORSE["logOR",j] -
                           1.96 * logORSE["SE",j],
                           angle=90, code=3, length=0.05)
                })
            }
        }
        
        return(logORSE);
    })
    names(responses) <- colnames(glmData)[2:(ncol(glmData)-1)]
    responsesL[[cancerID]] <- responses;
}
dev.off();

logORs <- matrix(NA, ncol(glmData) - 2, ncol(alnArchs))
rownames(logORs) <- colnames(glmData)[2:(ncol(glmData)-1)];
colnames(logORs) <- colnames(alnArchs)
    
cancerID <- cancerIDs[1]
treat <- colnames(glmData)[2:(ncol(glmData)-1)][1]
for (cancerID in cancerIDs) {
    cat(paste(cancerID, "\n"));
    nArchs <- length(grep(cancerID, colnames(alnArchs)));
    for (idx in 1:nArchs) {
        cat(paste(idx, "\n"));
        for (treat in colnames(glmData)[2:(ncol(glmData)-1)]) {
            ## cat(paste(treat, "\n"));
            if ( is.null(responsesL[[cancerID]][[treat]]) ) { next; }
            sel <-
                which(names(
                    responsesL[[cancerID]][[treat]]["logOR",]
                    ) == idx)
            if ( length(sel) == 0 ) { next; }
            logORs[treat, sprintf("%s.%d", cancerID, idx)] <-
                responsesL[[cancerID]][[treat]]["logOR",sel]
        }
    }
}

alnResponses <- logORs;
selTreats <- names(rev(sort(apply(!is.na(alnResponses), 1, sum)))[1:4])
## selArchs <- which(apply(!is.na(alnResponses), 2, sum) >= 3)
selArchs <-
    intersect(
        which(apply(!is.na(alnResponses), 2, sum) >= 1),
        setdiff(1:ncol(alnResponses), grep("SYNT", colnames(alnResponses)))
        )

## alnResponses0 <- t(apply(alnResponses, 1, function(x) { x - mean(x, na.rm=T) }))
alnResponses0 <- alnResponses[selTreats,selArchs]

alnResponses0NA <- alnResponses0;
alnResponses0NA[is.na(alnResponses0NA)] <- 0;

pdf("respPerArch.pdf", height=6, width=8);
heatmap.2(alnResponses0NA, scale="none", trace="none",
          col = heat.colors, margins = c(5, 12.5),
          hclustfun=function(x) { hclust(x, method="ward.D") })
dev.off()


##################################################

## Take out average gene expression to focus on changes between cancers
sel <- rev(order(apply(alnArchs, 1, sd)))[1:1000]
## sel <- 1:nrow(alnArchs)
topExprArchs <-
    t(apply(alnArchs[sel,], 1, function(x) { x - mean(x) }))

pdf("distRespExprArch.pdf", height=8, width=8)
par(mfrow=c(2,2))

dTreat <- as.matrix(dist(t(alnResponses0)))
dTreatN <- dTreat / median(dTreat)
image(dTreatN, main="Responses")

dExpr <- as.matrix(dist(t(topExprArchs[,selArchs])))
dExprN <- dExpr / median(dExpr)
image(dExprN, main="Expression")

## plot(dExprN, dTreatN)
## cor.test(dExprN, dTreatN)

## rndTEPR <-
##     t(apply(topExprArchs, 1, function(x) { sample(x) }))
## dRndTEPR <- as.matrix(dist(t(rndTEPR)))
## dRndTEPRN <- dRndTEPR / median(dRndTEPR)
## image(dRndTEPRN, main="Random expr.")

## Normalization factor should be the same as is Expression?

## rndDists <-
##     sapply(1:1000, function(i) {
##         rndTEPR <-
##             t(apply(topExprArchs, 1,
##                     function(x) { sample(x) }))
##         dRndTEPR <- as.matrix(dist(t(rndTEPR)))
##         dRndTEPRN <- dRndTEPR / median(dRndTEPR)
##         sum((dTreatN - dRndTEPRN)^2)
##     })

rnddExprN <- dExprN;
upperTriangle(rnddExprN, diag=F) <- 
    sample(upperTriangle(dExprN, diag=F))
lowerTriangle(rnddExprN) <- NA; ## upperTriangle(rnddExprN)
image(rnddExprN, main="Random expr.")

rndDists <-
    sapply(1:1000, function(i) {
        rnddExprN <- 
            sample(upperTriangle(dExprN, diag=F))
        sum((upperTriangle(dTreatN, diag=F) -
             rnddExprN)^2)
    })

initDist <- sum((upperTriangle(dTreatN) - upperTriangle(dExprN))^2)
hist(rndDists, 20,
     xlim=range(c(rndDists, initDist)),
     xlab="distance between matrices")
abline(v=initDist, lwd=2, col="red")
dev.off();

##################################################

## We test whether good responses are associated to specific
## archetypes (theta test)

cancerID <- "BLCA";
cancerID <- "BRCA";
cancerID <- "GBM";
cancerID <- "LGG";
cancerID <- "SYNT";
cancerID <- "ALL";

EVsL <- list();
load("allCancers_EVs.rda")
showEverything <- F;

testArchResp <- function(glmData, thetas,
                         what=c("p-value", "direction"), minOcc=5) {
    sapply(colnames(glmData)[2:ncol(glmData)], function(i) {
        cat(paste(i, "\n"))
        tgd <- glmData[glmData[,i] == 1, "best.resp"]
        thisTheta <- thetas[glmData[,i] == 1,]
        if ( sum(tgd == 0) > minOcc/2 && sum(tgd == 1) > minOcc/2 ) {
            ## plot(exprProj[,1], exprProj[,2],
            ##      xlim=range(c(exprProj[,1], archProj[,1])),
            ##      ylim=range(c(exprProj[,2], archProj[,2])),
            ##      main=cancerID)
            ## points(archProj, pch=20, col="blue", cex=2)
            
            ## points(exprProj[which(noNA)[which(glmData[,i] == 1)[tgd == 1]],1],
            ##        exprProj[which(noNA)[which(glmData[,i] == 1)[tgd == 1]],2],
            ##        col="green", pch=20);
            ## points(exprProj[which(noNA)[which(glmData[,i] == 1)[tgd == 0]],1],
            ##        exprProj[which(noNA)[which(glmData[,i] == 1)[tgd == 0]],2],
            ##        col="red", pch=20);

            lapply(1:ncol(thetas), function(arcIdc) {
                if ( what == "p-value" ) {
                    return(
                        wilcox.test(thisTheta[tgd == 1, arcIdc],
                                    thisTheta[tgd == 0, arcIdc])$p.value
                        )
                } else {
                    return(
                        (-(median(thisTheta[tgd == 1, arcIdc]) -
                           median(thisTheta[tgd == 0, arcIdc])))
                        )
                }
            })
        } else {
            return(lapply(1:ncol(thetas), function(arcIdc) {
                return(NA)
            }))
        }
    })
}

pdf("treatArchResponse.pdf", height=8, width=8);
pVals <-
    sapply(names(respCutOffs)[!is.na(respCutOffs)], function(cancerID) {
        cat(paste(cancerID,"\n"))
        par(mfrow=c(2,2))
        
        expr <- exprCancers[[cancerID]];
        myQuantile <- TCGAfracNarchs[cancerID,"quantile"];
        if ( is.na(myQuantile) || myQuantile == 0 ) {
            exprU <- t(expr);
        } else {
            meanGeneExpr <- apply(expr, 2, mean);
            minExpr <- quantile(meanGeneExpr, myQuantile);
            exprU <- t(expr[,meanGeneExpr >= minExpr])
        }
        rm(expr);
        
        ## discreteClinicalData_reOrdered_withTreatment.tsv, select
        ## only treat.target
        treatmentsFull <- 
            read.table(
                sprintf("%s_UCSC/discreteClinicalData_reOrdered.tsv",
                        cancerID), as.is=T, h=T, sep="\t")
        if ( any(colnames(treatmentsFull) == "pathologic_M") ) {
            Mstatus <- treatmentsFull[,"pathologic_M"]
            ## Nstatus <- treatments[,"pathologic_N"]
        } else {
            cat(sprintf("No metastasis information for %s\n", cancerID))
            Mstatus <- rep("M0", nrow(treatmentsFull))
        }
        treatments <- treatmentsFull[,c(which("treat.patientBestResp" == colnames(treatmentsFull)),
                                    grep("treat.target.", colnames(treatmentsFull)))]
        colnames(treatments) <- gsub("^treat.target.", "", colnames(treatments))
        colnames(treatments)[1] <- "best.resp";

        best.respVec <-
            sapply(treatments[,"best.resp"], function(x) {
                if ( is.na(x) ) { return(NA) }
                which(respRanking == x)
            })
        treatments[,"best.resp"] <- best.respVec;
        
        myArchs <- archPerCancerList[[cancerID]]
        if ( sum(rownames(myArchs) != rownames(exprU)) > 0 ) {
            stop("genes in archetypes and samples are mis-aligned.\n")
        }
        
        centerVec <- apply(exprU, 1, mean)
        exprU0 <-
            t(sapply(1:nrow(exprU), function(i) {
                exprU[i,] - centerVec[i] }))
        archs0 <-
            t(sapply(1:nrow(exprU), function(i) {
                myArchs[i,] - centerVec[i]
            }))
        rm(exprU); gc();

        nPCsInit <- 10;
        nPCs <- 4;
        if ( is.null(EVsL[[cancerID]]) ) {
            EVsL[[cancerID]] <- svd(exprU0, nu=nPCsInit, nv=nPCsInit)$u[,1:nPCsInit];
            save(EVsL, file="allCancers_EVs.rda")
        } else {
            cat("Reusing previously computed eigenvectors.\n")
        }
        EVs <- EVsL[[cancerID]][,1:nPCs];

        exprProjInit <- t(t(EVsL[[cancerID]]) %*% exprU0)
        write.table(exprProjInit,
                    file=sprintf("%s_samples_10PCs.csv", cancerID),
                    quote=F, row.names=F, col.names=F, sep=",");
        rm(exprProjInit); gc();

        archProjInit <- t(t(EVsL[[cancerID]]) %*% archs0)
        write.table(archProjInit,
                    file=sprintf("%s_archetypes_10PCs.csv", cancerID),
                    quote=F, row.names=F, col.names=F, sep=",");
        rm(archProjInit);
        
        exprProj <- t(t(EVs) %*% exprU0)
        archProj <- t(t(EVs) %*% archs0)
        
        plot(exprProj[,1], exprProj[,2],
             xlim=range(c(exprProj[,1], archProj[,1])),
             ylim=range(c(exprProj[,2], archProj[,2])),
             main=cancerID)
        points(archProj, pch=20, col="blue", cex=2)

        plot3d(exprProj[,1], exprProj[,2], exprProj[,3],
               size=3, type="p", alpha=.5,
               col="red", xlab="", ylab="", zlab="")
        points3d(archProj[,1], archProj[,2], archProj[,3],
                 col="blue", size=8)

        ## Compute distance from each tumor to each archetype
        distTtoA <-
            t(apply(exprProj[,1:3], 1, function(tumor) {
                apply(archProj[,1:3], 1, function(arch) {
                    sqrt(sum((tumor - arch)^2))
                })
            }))
        rownames(distTtoA) <-
            read.table(sprintf("%s_UCSC/patientIDs.list", cancerID))[,1]
        write.csv(distTtoA,
                  file=sprintf("%s_dist_tumor_to_arch.csv", cancerID))
        
        isPRAD <- treatmentsFull[,"cancer_type"] == "BRCA"
        plot3d(exprProj[isPRAD,1], exprProj[isPRAD,2],
               exprProj[isPRAD,3], size=.5, type="s",
               col="red", xlab="", ylab="", zlab="")
        points3d(archProj[,1], archProj[,2], archProj[,3],
                 col="blue", size=8)
        rgl.points(exprProj[!isPRAD,1], exprProj[!isPRAD,2],
                   exprProj[!isPRAD,3], size=2, col="black", alpha=.25)
        
        ## Move the last archetype to 0,0:
        exprProjA <- t(apply(exprProj, 1, function(x) { x - archProj[nrow(archProj),] }))
        archProjA <-
            t(apply(archProj, 1,
                    function(x) { x - archProj[nrow(archProj),] }))[-nrow(archProj),]
        ## archetypes x dimensions
        thetas <- exprProjA[,1:nrow(archProjA)] %*% solve(archProjA[,1:nrow(archProjA)])
        thetas <- cbind(thetas, apply(thetas, 1, function(x) { 1 - sum(x) }))

        minOcc <- 10;
        if ( cancerID == "ALL" ) { minOcc <- 100; }
        treatOcc <- sort(apply(treatments[,-1], 2, function(x) { sum(x, na.rm=T) }))
        keptTreats <- names(which(treatOcc > minOcc))
        if ( length(keptTreats) == 0 ) { return(NULL) }
        treatFact <- as.data.frame(treatments[,keptTreats])
        for ( i in 1:ncol(treatFact) ) {
            treatFact[,i] <- as.factor(treatFact[,i])
        }
        colnames(treatFact) <- keptTreats
        
        table(treatments[,"best.resp"])

        if ( cancerID == "ALL" ) {
            ## x <- unique(treatmentsFull[,"cancer_type"])[1]
            ## isBestResp <- 
            ##     unlist(
            ##         sapply(unique(treatmentsFull[,"cancer_type"]), function(x) {
            ##             ## sum(treatmentsFull[,"cancer_type"] == x)
            ##             treatments[treatmentsFull[,"cancer_type"] == x,
            ##                        "best.resp"] > respCutOffs[x]
            ##         })
            ##         )
            isBestResp <- treatments[,"best.resp"] > median(respCutOffs, na.rm=T)
        } else {
            respCutOff <- respCutOffs[cancerID]
            isBestResp <- treatments[,"best.resp"] > respCutOff;
        }
        preGlmData <-
            data.frame(
                ## best.resp=as.numeric(treatments[,"best.resp"] !=
                ##     "Clinical Progressive Disease"),
                best.resp=isBestResp,
                ## best.resp=as.numeric(treatments[,"best.resp"] == "Complete Response"),
                treatFact)
        ## glmData <- glmData[apply(is.na(glmData), 1, sum) == 0,]
        hist(apply(is.na(preGlmData), 1, sum))
        cat(sprintf("%d patients with known response.\n",
                    sum(!is.na(preGlmData[,1]))))
        noNA <- !apply(is.na(preGlmData), 1, any) ## & Mstatus == "M0"
        glmData <- preGlmData[noNA,]
        if ( nrow(glmData) < minOcc ) {
            return(NULL);
        }
        apply(glmData[,-1], 2, function(x) { sum(as.numeric(x)) })
        summary(glmData)
        i <- "nucleotide.depletion";
        i <- "DNA.damaging";
        i <- "angiogenesis.signalling.inhibitor";
        i <- "microtubule";
        i <- "growth.signalling.inhibitor";
        i <- "topoisomerase.inhibitor";
        i <- "other";

        dir <- testArchResp(glmData, thetas[noNA,], what="direction", minOcc=minOcc)
        p <- testArchResp(glmData, thetas[noNA,], "p-value", minOcc=minOcc)
        ## q <- as.numeric(p) * sum(!is.na(p)) #MTC arch x drugs
        q <- as.numeric(p) * sum(!apply(is.na(p), 2, any)) #MTC drugs
        ## q <- as.numeric(p) * nrow(p) #MTC arch
        q[q>1] <- 1

        q <- matrix(q, nrow(p), ncol(p))
        colnames(q) <- colnames(p);
        return(list(dir=dir, p=p, q=q))
    })
dev.off();

cancerID <- cancerIDs[1]
cancerID <- "SYNT"
cancerID <- "BLCA"
cancerID <- "GBM"
## sapply(colnames(pVals), function(cancerID) {
sapply(names(pVals), function(cancerID) {
    cat(paste(cancerID, "\n"))

    ## qMat <- pVals[["p", cancerID]]
    qMat <- pVals[[cancerID]]
    if ( is.null(qMat) ) { return() }
    
    i <- 1;
    sapply(1:nrow(qMat[["p"]]), function(i) {
        j <- 1;
        sapply(1:ncol(qMat[["p"]]), function(j) {
            if ( !is.na(qMat[["p"]][i,j]) && qMat[["p"]][i,j] < .15 ) {
                ## browser();
                ## cat(sprintf("%s: arch %d - %s (%.2f)\n",
                ##             cancerID, i, colnames(qMat)[j],
                ##             pVals[["q", cancerID]][i,j]))
                cat(sprintf("%s: arch %d - %s (dir=%.2f, p=%.3f)\n",
                            cancerID, i, colnames(qMat[["p"]])[j],
                            pVals[[cancerID]][["dir"]][i,j],
                            pVals[[cancerID]][["p"]][i,j]
                            ))
            }

        })
        return()
    })
    return()
})
