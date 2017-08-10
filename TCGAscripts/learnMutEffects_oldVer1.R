rm(list=ls()); gc();

isMetabric <- length(grep("metabric", getwd())) == 1

geneNames <- read.table("geneListExp.list", as.is=T, h=F, sep="\t")[,1]
geneNamesFilt <- read.table("geneNamesAfterExprFiltering.list",
                            as.is=T, h=F, sep="\t")[,1]
geneExpression <- read.csv("expMatrix.csv", as.is=T, h=F)
colnames(geneExpression) <- geneNames;
rm(geneNames);
tmp <- sapply(geneNamesFilt, function(g) {
    geneExpression[,g]
})
geneExpression <- tmp
dim(geneExpression)

mutations <- read.csv("mutMatrix_reOrdered_booleanized_justData.csv",
                      as.is=T, h=F)
mutGeneNames <-
    read.table("mutMatrix_reOrdered_booleanized_geneNames.list",
               h=F, as.is=T)[,1]
colnames(mutations) <- mutGeneNames;
rm(mutGeneNames)
dim(mutations)
mutations <- as.matrix(mutations)

if ( isMetabric ) {
    isHealthy <- apply(is.nan(mutations), 1, sum) == ncol(mutations)
} else {
    ## NaNcol <- grep("=NaN", colnames(mutations));
    ## isHealthy <- as.logical(apply(mutations[,NaNcol], 1, mean))
    ## mutations <- mutations[,setdiff(1:ncol(mutations), NaNcol)]
    ## rm(NaNcol)
    discClin <- read.table("discreteClinicalData_reOrdered.tsv",
                           as.is=T, sep="\t", h=T)
    isHealthy <- discClin[,"sample_type"] == "Solid Tissue Normal"
    rm(discClin)
}

CNAs <- read.csv("copMatrix_reOrdered_booleanized_justData.csv",
                 as.is=T, h=F)
copGeneNames <-
    read.table("copMatrix_reOrdered_booleanized_geneNames.list",
               h=F, as.is=T)[,1]
colnames(CNAs) <- copGeneNames
x <- copGeneNames[1]
copTab <- t(sapply(copGeneNames, function(x) {
    strsplit(x, "=")[[1]]
}))
## if ( isMetabric ) {
##     CNAgenes <- ## from Cis-acting CNAs
##         as.character(
##             read.table("/home/jhausser/Desktop/projects/cancerTaskAtlas/brca_metabric/top1000cisActingGenes.list",
##                        as.is=T, h=F, sep="\t")[,1])
## } else {
##     CNAgenes <- ## from Chris Sanders
##         as.character(unlist(
##             sapply(
##                 read.csv("/home/jhausser/work/cancerTaskAtlas/Ciriello2013_ng.2762-S2.csv",
##                          as.is=T)[,"Genes.in.GISTIC.peak"],
##                 function(x) {
##                     strsplit(x, ",")[[1]]
##                 })))
## }
CNAgenes <- ## from Cis-acting CNAs
    as.character(
        read.table("/home/jhausser/Desktop/projects/cancerTaskAtlas/brca_metabric/top1000cisActingGenes.list",
                   as.is=T, h=F, sep="\t")[,1])

CNAgenes <- unique(CNAgenes)
x <- CNAgenes[1]
sel <- sort(unique(unlist(sapply(CNAgenes, function(x) {
    grep(sub("$", "=", x), copGeneNames)
}))))
CNAsSmall <- CNAs[,sel]
colnames(CNAsSmall) <- copGeneNames[sel]
plot(rev(sort(apply(CNAsSmall, 2, mean)))); abline(v=600, lty=2)
colnames(CNAsSmall) <- sprintf("%s_CNA", colnames(CNAsSmall));
rm(copGeneNames)
CNAs <- CNAsSmall;
rm(CNAsSmall)
## CNAs <- CNAs[,rev(order(apply(CNAs, 2, mean)))]
## plot(apply(CNAs, 2, mean)); abline(v=4000, lty=2)
## CNAs <- CNAs[,1:4000]

healthyProfile <- apply(geneExpression[isHealthy,], 2, median)

gE0 <- apply(geneExpression, 1, function(x) { x - mean(x) })
rm(geneExpression)
## gE00 <- apply(gE0, 1, function(x) { x - mean(x) })
gE00 <- apply(gE0, 2, function(x) { x - healthyProfile } )
rm(gE0); gc();

library(MASS)
M <- t(mutations)
rm(mutations)
Mp <- rbind(M, t(CNAs))
## Mp <- M
Mp <- t(CNAs);

rm(CNAs)
E <- gE00
rm(gE00)

## save.image("learnMutEffects.rda");
rm(list=ls()); gc();
load("learnMutEffects.rda");

if ( isMetabric ) {
    E <- E[,apply(is.na(M), 2, sum) == 0]
    M <- M[,apply(is.na(M), 2, sum) == 0]
}
## Mp <- M

## 50 is the tikhonov factor (lambda)
## X <-
##     t(
##     (solve( M %*% t(M) + 50 * diag(1, nrow=dim(M)[1]) ) ) %*%
##     M %*% t(E)
##     )
X <-
    t(
    (solve( Mp %*% t(Mp) + 50 * diag(1, nrow=dim(Mp)[1]) ) ) %*%
    Mp %*% t(E)
    )
hist(X, seq(-max(abs(X)), max(abs(X)), len=40))
hist(apply(X, 2, sd), 20)

dim(E)
dim(Mp)
dim(diag(1 / apply(Mp, 1, sum)))

## At least 10 examples so vectors aren't too noisy
commonMuts <- intersect(which(apply(Mp, 1, sum) >= 10),
                        setdiff(1:nrow(Mp), grep("=NaN_", rownames(Mp))))
## Xd <-
##     E %*% t(Mp) %*% diag(1 / apply(Mp, 1, sum)) -
##     E %*% t(1-Mp) %*% diag(1 / apply(1-Mp, 1, sum))
Xd <-
    E %*% t(Mp[commonMuts,]) %*% diag(1 / apply(Mp[commonMuts,], 1, sum)) -
    E %*% t(1-Mp[commonMuts,]) %*% diag(1 / apply(1-Mp[commonMuts,], 1, sum))
colnames(Xd) <- rownames(Mp)[commonMuts]
dim(Xd)
Xd[1:5,1:5]
hist(apply(Mp, 1, sum))
sum(apply(Mp, 1, sum) >= 10) / nrow(Mp)
dim(Mp)

## ## At least 10 examples so vectors aren't too noisy
## commonMuts <- intersect(which(apply(Mp, 1, sum) >= 10),
##                         setdiff(1:nrow(Mp), grep("=NaN_", rownames(Mp))))
## m <- commonMuts[1]
## Xd <- sapply(commonMuts, function(m) {
##     cat(sprintf("%d\n", m))
##     clab <- sub("=.*_CNA", "=0_CNA", rownames(Mp)[m])
##     ## Mp[m,]
##     E %*% matrix(Mp[m,], ncol=1) / sum(Mp[m,]) -
##         E %*% matrix(Mp[clab,], ncol=1) / sum(Mp[clab,])
## })
## colnames(Xd) <- rownames(Mp)[commonMuts]
## Xd[1:5,1:5]

## x <- rownames(Mp)[1]
## Xd <-
##     sapply(## rownames(Mp)
##         1:nrow(Mp), function(x) {
##             cat(sprintf("%s\n", x))
##             apply(E[,as.logical(Mp[x,])], 1, mean) -
##                 apply(E[,!as.logical(Mp[x,])], 1, mean)
##         })
## Xd <- mcapply(Mp, 1, function(x) {
##     sel <- as.logical(x)
##     if ( sum(sel) < 10 ) { return(rep(NA, nrow(E))) }
##     cat(sprintf("%d\n", sum(sel)))
##     apply(E[,sel], 1, mean) - apply(E[,!sel], 1, mean)
## })
## Xd[1:5,1:5]

## save(Xd, file="Xdelta.rda");
load("Xdelta.rda");
X <- Xd;

## X <- X[,commonMuts]

## These are the mutations with strongest impact on gene expression
rev(sort(apply(X, 2, sd)))[1:10]

library(ade4)
dudi1 <- dudi.pca(t(X), scale=F, center=F, scannf=F, nf=3)
cumsum(dudi1$eig)[1:5] / sum(dudi1$eig)
library(rgl)
plot3d(dudi1$li)

mutWithStrongestEffect <-
    names(rev(sort(apply(dudi1$li[,2:3]^2, 1, sum)))[1:8])
## mutWithStrongestEffect <-
##     names(rev(sort(apply(dudi1$li[,2:3]^2, 1, sum)))[1:2])
dudi1$li[mutWithStrongestEffect,2:3]

## first axis is healthy -> cancer axis
pdf("mutShiftsPCA.pdf", height=4, width=4)
plot(dudi1$li[,2:3],
     xlim=range(dudi1$li[,2]) * 1.1,
     ylim=range(dudi1$li[,3]) * 1.1)
## text(dudi1$li[mutWithStrongestEffect,2:3],
##      sub("=1", "", mutWithStrongestEffect), pos=3)
s.arrow(dudi1$li[mutWithStrongestEffect,2:3],
        label=sub("=1$", "", mutWithStrongestEffect), add.plot=T)
dev.off();
rm(dudi1)

dudi2 <- dudi.pca(t(E), scale=F, scannf=F, nf=3)

## Xproj <- t(X) %*% as.matrix(dudi2$co) %*%
##     diag(1 / apply(dudi2$co, 2, function(x) { (sum(x^2)) }))

avgE <- apply(E, 1, mean)
## Compute the end points of mutations
## X <- cbind(X, CTRL=avgE)

## save.image("learnMutEffects2.rda")
## load("learnMutEffects2.rda")

Xproj <-
    t(apply(X, 2, function(x) { x - avgE })) %*%
    as.matrix(dudi2$c1)
## Xproj <-
##     t(apply(Xd, 2, function(x) { x - avgE })) %*%
##     as.matrix(dudi2$c1)

## Norms of the PCs:
## apply(dudi2$c1, 2, function(x) { sqrt(sum(x^2)) }) ## 1, good
Xproj[mutWithStrongestEffect,]
## XbackProj <- as.matrix(dudi2$co) %*% t(Xproj)
## tmp <- as.matrix(dudi2$c1) %*% t(Xproj)
## XbackProj <- apply(tmp, 2, function(x) { x + avgE })
XbackProj <- as.matrix(dudi2$c1) %*% t( t(X) %*% as.matrix(dudi2$c1) )
plot(X[,mutWithStrongestEffect[1]], XbackProj[,mutWithStrongestEffect[1]]);
plot(X[,mutWithStrongestEffect[1]] - mean(X[,mutWithStrongestEffect[1]]),
     XbackProj[,mutWithStrongestEffect[1]]);
abline(0,1, col="grey");
plot(X[,1], XbackProj[,1]);
abline(0,1, col="grey");

plot(X[,1] - mean(X[,1]), XbackProj[,1]);
abline(0,1, col="grey");

sum(((X[,1] - mean(X[,1])) - XbackProj[,1]) * XbackProj[,1])
sum(((X[,mutWithStrongestEffect[1]] - mean(X[,mutWithStrongestEffect[1]])) - XbackProj[,mutWithStrongestEffect[1]]) * XbackProj[,mutWithStrongestEffect[1]])

res <- sapply(colnames(X), function(i) {
    ##residuals: SS off Pareto front
    sum(((X[,i] - mean(X[,i])) - XbackProj[,i])^2)
})
onf <- sapply(colnames(X), function(i) {
    #SS on Parto front
    sum(XbackProj[,i]^2)
})
## tot <- sapply(colnames(X), function(i) {
##     ##residuals: SS off Pareto front
##     sum((X[,i] - mean(X[,i]))^2)
## })
## summary(res + onf - tot)
rev(sort(onf/res))[1:50]
ratios <- log10(onf/res);
sum(dudi2$eig[1:3])/sum(dudi2$eig)
ratios <- onf/(onf+res)
hist(ratios)
abline(v=sum(dudi2$eig[1:3])/sum(dudi2$eig), lty=2)
rm(onf); rm(res);

## Randomize mutations and see what the onf/res ratio is for random
## mutations
library(parallel)
SSlist <-
    mclapply(1:10, function(i) {
        cat(sprintf("Iteration %d\n", i))

        ## Mr <- t(apply(Mp, 1, sample))
        ## Xr <-
        ##     t(
        ##     (solve( Mr %*% t(Mr) + 50 * diag(1, nrow=dim(Mr)[1]) ) ) %*%
        ##     Mr %*% t(E)
        ##     )

        ## m <- commonMuts[1]
        ## Xr <- sapply(commonMuts, function(m) {
        ##     cat(sprintf("Iteration %d: %d\n", i, m))
        ##     clab <- sub("=.*_CNA", "=0_CNA", rownames(Mr)[m])
        ##     ## Mp[m,]
        ##     E %*% matrix(Mr[m,], ncol=1) / sum(Mr[m,]) -
        ##         E %*% matrix(Mr[clab,], ncol=1) / sum(Mr[clab,])
        ## })
        ## colnames(Xr) <- rownames(Mp)[commonMuts]

        Mr <- t(apply(Mp[commonMuts,], 1, sample))
        Xr <-
            E %*% t(Mr) %*% diag(1 / apply(Mr, 1, sum)) -
            E %*% t(1-Mr) %*% diag(1 / apply(1-Mr, 1, sum))
        colnames(Xr) <- rownames(Mr)
        
        XrbackProj <- as.matrix(dudi2$c1) %*% t( t(Xr) %*% as.matrix(dudi2$c1) )

        resr <- sapply(colnames(X), function(i) {
            ##residuals: SS off Pareto front
            sum(((Xr[,i] - mean(Xr[,i])) - XrbackProj[,i])^2)
        })
        onfr <- sapply(colnames(X), function(i) {
            ##SS on Parto front
            sum(XrbackProj[,i]^2)
        })
        return(list(resr=resr, onfr=onfr))
    })
## onfr <- unlist(SSlist["onfr",])
onfr <- (sapply(SSlist, function(s) { s[["onfr"]] }))
## resr <- unlist(SSlist["resr",])
resr <- (sapply(SSlist, function(s) { s[["resr"]] }))
## rratios <- log10(onfr/resr);
rratios <- onfr/(onfr+resr);
rratios[1:5,1:5]
 
## pVals <- sort(sapply(ratios, function(x) {
##     sum(x < log10(onfr/resr)) / length(onfr)
## }))

## Let's determine p-values by comparing each mutation to its own
## randomized versions
sigmas <-
    sapply(names(ratios), function(m) {
        ( ratios[m] - mean(rratios[m,]) ) / sd(rratios[m,])
    })
library(fdrtool)
names(sigmas) <- names(ratios)
pVals <- pnorm(- abs(sigmas)) * 2;
hist(pVals)
qVals <- fdrtool(pVals, statistic="pvalue")$q
qVals <- pVals * length(pVals)
qVals[qVals>1] <- 1
hist(qVals)
sort(qVals)[1:10]
sum(qVals < .01)
sort(sigmas[qVals < .01])
plot(sigmas, log10(qVals))

mybreaks <- seq(min(c(ratios, rratios)), max(c(ratios, rratios)), len=20);
pdf("fracSSonVsOff.pdf", height=4, width=4)
hist(ratios,
     breaks=mybreaks, main="",
     col="lightblue", freq=F, ## ylim=c(0,5),
     ## xlab=expression(paste(log[10], " (SSon / SSoff)")),
     xlab="SSon / (SSon + SSoff)")
hist(onfr/(onfr+resr), breaks=mybreaks, border="red",
     col=NULL, add=T, freq=F)
abline(v=sum(dudi2$eig[1:3])/sum(dudi2$eig), lty=2)
legend("topleft", c("measured", "shuffled"),
       fill=c("lightblue", "white"), col=c("black", "red"),
       bg="white")
## ## sapply(names(pVals[pVals<.01 | pVals>.99]), function(x) {
## sapply(names(sigmas[qVals<.01 & sigmas > 0]), function(x) {
## ## sapply(c("PVT1=1_CNA", "KARS=-1_CNA"), function(x) {
##     text(ratios[x], .3, sub("=1", "", x), srt=90)
##     arrows(ratios[x], .13, ratios[x], .02, len=.07)
## })
dev.off();

mutWithStrongestEffect <- names(sigmas[qVals<.01 & sigmas > 0])
## mutWithStrongestEffect <- names(sigmas[qVals<.1 & sigmas > 0])
## mutWithStrongestEffect <- names(sigmas[qVals<.05 & sigmas > 0])

## These genes mutate outside the Pareto front:
pVals[pVals>.999]
## Protocadherin Fat 4, aka FAT4 tumor suppressor; inhibits YAP1
sort(pVals[pVals<.001])

## TP53 and GATA3 mutate inside the cancer space significantly more
## than random mutations.
## TP53 = Basal [2] (IntCl 2) DNA replication, Chromatin, patients 6 yrs
## older patients, many tumor cells and necrosis, Infiltrating Ductal Carcinoma
## GATA3 = LumB [3] (IntCl 3,4) ER+, Peroxisome (Mucinous Carcinoma?)

## CDH1 = unchar archetype [1]: Infiltrating Lobular Carcinoma, Complete
## Response, Locoregional Recurrence, postoperative radiotherapy, new
## tumor event after initial treatment: NO, ER+/PR+, people not
## included in 2012 study; translation / respiration, proteasome

arcsOrig <- read.table("arcsOrig_genes_tabSep.tsv", sep="\t", h=F, as.is=T)[,-1]
plot(as.numeric(arcsOrig[4,]), as.numeric(healthyProfile)); abline(0,1,col="grey")

## genesFilt <- read.table("geneNamesAfterExprFiltering.list", as.is=T)[,1]
## colnames(arcsOrig) <- genesFilt;

apply(arcsOrig, 1, length)
archProj <- t(sapply(1:4, function(i) {
    m0 <- arcsOrig[i,] - healthyProfile - avgE
    as.matrix(m0) %*% as.matrix(dudi2$c1)
}))

hPt <-
    (- matrix(apply(E, 1, mean), nrow=1) ) %*% as.matrix(dudi2$co) %*%
    diag(1 / apply(dudi2$co, 2, function(x) { sqrt(sum(x^2)) }))

s <- 10;
s <- 3;
s <- 5;

## pdf("mutVectors.pdf", height=6*2, width=6);
## par(mfrow=c(2,1))
## plot(dudi2$li[,1:2], col="grey", pch=20,
##      xlim=range(archProj[,1]), ylim=range(archProj[,2]))
## points(archProj[,1:2], col="red", cex=3)
## text(archProj[,1], archProj[,2], 1:4)
## points(dudi2$li[isHealthy,1:2], col="green")
## ## pH <- apply(dudi2$li[isHealthy,], 2, median)
## ## points(pH[1], pH[2], pch=20, col="red", cex=2)

## points(hPt[1], hPt[2], col="green", pch=20, cex=2)
## sapply(mutWithStrongestEffect, function(x) {
##     v <- c(Xproj[x,1] - hPt[1], Xproj[x,2] - hPt[2]);
##     arrows(hPt[1], hPt[2],
##            hPt[1] + s * v[1], hPt[2] + s * v[2],
##            len=.1)
##     text(hPt[1] + s * v[1], hPt[2] + s * v[2],
##          sub("=1", "", x), pos=4)
## })

## legend("bottomleft", "healthy", pch=1, col="green");

## plot(dudi2$li[,c(1,3)], col="grey", pch=20,
##      xlim=range(archProj[,1]), ylim=range(archProj[,3]))
## points(archProj[,c(1,3)], col="red", cex=3)
## text(archProj[,1], archProj[,3], 1:4)
## points(dudi2$li[isHealthy,c(1,3)], col="green")
## ## pH <- apply(dudi2$li[isHealthy,], 2, median)
## ## points(pH[1], pH[2], pch=20, col="red", cex=2)

## points(hPt[1], hPt[3], col="green", pch=20, cex=2)
## x <- mutWithStrongestEffect[1]
## sapply(mutWithStrongestEffect, function(x) {
##     v <- c(Xproj[x,1] - hPt[1], Xproj[x,2] - hPt[2], Xproj[x,3] - hPt[3]);
##     arrows(hPt[1], hPt[3],
##            hPt[1] + s * v[1], hPt[3] + s * v[3],
##            len=.1)
##     text(hPt[1] + s * v[1], hPt[3] + s * v[3],
##          sub("=1", "", x), pos=4)
## })
## dev.off();

hPt

## save.image("plot3d.rda")
## load("plot3d.rda")
plot3d(dudi2$li, alpha=.35, col="grey", box=F)
spheres3d(hPt, col="green", radius=5)
spheres3d(archProj, col="blue", radius=5)
text3d(archProj[,1], archProj[,2], archProj[,3], 1:4, adj=2.5)

sapply(1:3, function(i) {
    sapply(seq(i+1, 4), function(j) {
        segments3d(archProj[c(i,j),1], archProj[c(i,j),2],
                   archProj[c(i,j),3], col="blue", alpha=.5)
    })
})

## Add mutations
sapply(mutWithStrongestEffect, function(m) {
    x <- Xproj[m,1:3]
    v <- x - hPt;
    
    lines3d(c(hPt[1], s * v[1] + hPt[1]),
            c(hPt[2], s * v[2] + hPt[2]),
            c(hPt[3], s * v[3] + hPt[3]), col="black", lwd=2)
    ## text3d(s * v[1] + hPt[1], s * v[2] + hPt[2],
    ##        s * v[3] + hPt[3], sub("_CNA$", "", sub("=1$", "", m)), adj=c(.5, 0))
})

sapply(c("PVT1=1_CNA", "KARS=-1_CNA"), function(m) {
    x <- Xproj[m,1:3]
    v <- x - hPt;
    
    lines3d(c(hPt[1], s * v[1] + hPt[1]),
            c(hPt[2], s * v[2] + hPt[2]),
            c(hPt[3], s * v[3] + hPt[3]), col="black", lwd=2)
    text3d(s * v[1] + hPt[1], s * v[2] + hPt[2],
           s * v[3] + hPt[3], sub("_CNA$", "", sub("=1$", "", m)), adj=c(.5, 0))
})


rev(sort(apply(Xproj^2, 1, sum)))[1:10]
sapply(rownames(Xproj), function(m) {
    x <- Xproj[m,1:3]
    v <- x - hPt;
    
    lines3d(c(hPt[1], s * v[1] + hPt[1]),
            c(hPt[2], s * v[2] + hPt[2]),
            c(hPt[3], s * v[3] + hPt[3]), col="red", lwd=2,
            alpha=.5)
    ## text3d(s * v[1] + hPt[1], s * v[2] + hPt[2],
    ##        s * v[3] + hPt[3], sub("_CNA$", "", sub("=1$", "", m)), adj=c(.5, 0))
})

##################################################

## Show location of Integrative Clusters -> could they be species?

## discClin <-
##     read.table("discreteClinicalData_reOrdered_withTreatment.tsv",
##                sep="\t", h=T)
## summary(discClin[,7:9])

IntClustAnnot <- read.table("TCGAIntClust1100.txt", sep="\t", h=T, as.is=T)
rownames(IntClustAnnot) <- IntClustAnnot[,1]
patientIDs <- read.table("patientIDs.list", h=F, as.is=T)[,1]
ICs <- sapply(patientIDs, function(x) { IntClustAnnot[x,"IntClust"] })
ICsL <- sort(unique(ICs)); ICsL <- ICsL[!is.na(ICsL)]



plot(dudi2$li[,1:2], col="grey")

pdf("mutVectorsICs.pdf", height=6*2, width=6);
par(mfrow=c(2,1))
plot(dudi2$li[,1:2], col="grey", pch=20,
     xlim=range(archProj[,1]), ylim=range(archProj[,2]))
points(archProj[,1:2], col="red", cex=3)
text(archProj[,1], archProj[,2], 1:4)
## points(dudi2$li[isHealthy,1:2], col="green")
## pH <- apply(dudi2$li[isHealthy,], 2, median)
## points(pH[1], pH[2], pch=20, col="red", cex=2)
sapply(ICsL, function(IC) {
    points(dudi2$li[ICs == IC,1:2], col=IC, pch=20)
})

plot(dudi2$li[,c(1,3)], col="grey", pch=20,
     xlim=range(archProj[,1]), ylim=range(archProj[,3]))
points(archProj[,c(1,3)], col="red", cex=3)
text(archProj[,1], archProj[,3], 1:4)
sapply(ICsL, function(IC) {
    points(dudi2$li[ICs == IC,c(1,3)], col=IC, pch=20)
})

legend("topright", legend=ICsL, col=ICsL, pch=20)
dev.off();

##

plot3d(dudi2$li, alpha=.35, col=ICs)
plot3d(dudi2$li, alpha=.35, col="grey", box=F)
## spheres3d(hPt, col="green", radius=5)
spheres3d(archProj, col="blue", radius=5)
text3d(archProj[,1], archProj[,2], archProj[,3], 1:4, adj=2.5)

sapply(1:3, function(i) {
    sapply(seq(i+1, 4), function(j) {
        segments3d(archProj[c(i,j),1], archProj[c(i,j),2],
                   archProj[c(i,j),3], col="blue", alpha=.5)
    })
})

IC <- ICsL[1]
sapply(ICsL, function(IC) {
    ## plot3d(ellipse3d(cov(dudi2$li[IC == ICs,], use="pairwise.complete"),
    ##                  centre=c(0,0,0), level=0.95), add=T)
    midPoint <- apply(dudi2$li[IC == ICs,], 2, function(x) { mean(x, na.rm=T) })
    plot3d(
        ellipse3d(cov(dudi2$li[IC == ICs,], use="pairwise.complete"),
                  centre=midPoint, level=0.66, col=IC, alpha=.1),
        add=T)
    text3d(midPoint[1], midPoint[2], midPoint[3], IC)
})
