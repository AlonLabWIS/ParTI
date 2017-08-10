rm(list=ls()); gc();
source("~/work/jeanLib.R")

altType <- "SNVs"; #SNVs, CNAs
altType <- "CNAs"; #SNVs, CNAs
isMetabric <- length(grep("metabric", getwd())) == 1

library(tidyverse)
geneNames <- read.delim("geneListExp.list", sep="\t", h=F, as.is=T)[,1]
geneNamesFilt <- read.table("geneNamesAfterExprFiltering.list",
                            as.is=T, h=F, sep="\t")[,1]
geneExpression <- read.csv("expMatrix.csv", as.is=T, h=F)
## geneExpression <- read_csv("expMatrix.csv", col_names=F)
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

discClin <- read.table("discreteClinicalData_reOrdered.tsv",
                       as.is=T, sep="\t", h=T)
if ( isMetabric ) {
    ## FIXME --- also look for samples with 0 mutations, not just NaN?
    samplesMutProfiled <- read.table("case_lists_sequenced.list", h=F, as.is=T)[,1]
    patientIDs <- read.table("patientIDs.list", h=F, as.is=T)[,1]
    isMutProfiled <- sapply(patientIDs, function(p) { p %in% samplesMutProfiled })
    isHealthy <-
        apply(is.nan(mutations), 1, sum) == ncol(mutations) & isMutProfiled
    healthyArch <- 2;
} else {
    ## NaNcol <- grep("=NaN", colnames(mutations));
    ## isHealthy <- as.logical(apply(mutations[,NaNcol], 1, mean))
    ## mutations <- mutations[,setdiff(1:ncol(mutations), NaNcol)]
    ## rm(NaNcol)
    isHealthy <- discClin[,"sample_type"] == "Solid Tissue Normal"
    healthyArch <- 4;
}
rm(discClin)

## CNAs <- read.csv("copMatrix_reOrdered_justData.tsv",
##                  as.is=T, h=F)
## copGeneNames <-
##     read.table("copMatrix_reOrdered_booleanized_geneNames.list",
##                h=F, as.is=T)[,1]
## colnames(CNAs) <- copGeneNames
CNAs <- read_delim("copMatrix_reOrdered_booleanized.tsv", delim="\t", col_names=T)
CNAs <- as.data.frame(CNAs)
rownames(CNAs) <- CNAs[,1]
CNAs <- CNAs[,-1]

## CNAs <- read_delim("copMatrix_reOrdered_booleanized_all.tsv",
##                    delim="\t", col_names=T)
## write_rds(CNAs, "CNAs_all.rds")
## CNAs <- read_rds("CNAs_all.rds")
## CNAs <- as.data.frame(CNAs);
## rownames(CNAs) <- CNAs[,1]
## CNAs <- CNAs[,-1]

copGeneNames <- colnames(CNAs)

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

healthyProfile <- apply(geneExpression[isHealthy,], 2, median)

gE0 <- apply(geneExpression, 1, function(x) { x - mean(x) })
rm(geneExpression)
## gE00 <- apply(gE0, 1, function(x) { x - mean(x) })
gE00 <- apply(gE0, 2, function(x) { x - healthyProfile } )
rm(gE0); gc();

library(MASS)
## M <- t(mutations)
## Mp <- rbind(M, t(CNAs))
## Mp <- M
if ( altType == "SNVs" ) {
    Mp <- t(mutations);
} else if ( altType == "CNAs" ) {
    Mp <- t(CNAs);
} else {
    stop("Unknown altType!");
}

rm(mutations)
rm(CNAs)
E <- gE00
rm(gE00, tmp)

## save.image(sprintf("learnMutEffects_%s.rda", altType))
rm(list=setdiff(ls(), "altType")); gc();
load(sprintf("learnMutEffects_%s.rda", altType))

## if ( isMetabric ) {
##     E <- E[,apply(is.na(Mp), 2, sum) == 0]
##     Mp <- Mp[,apply(is.na(Mp), 2, sum) == 0]
## }
## Mp <- M

## 50 is the tikhonov factor (lambda)
## X <-
##     t(
##     (solve( M %*% t(M) + 50 * diag(1, nrow=dim(M)[1]) ) ) %*%
##     M %*% t(E)
##     )
## hist(X, seq(-max(abs(X)), max(abs(X)), len=40))
## hist(apply(X, 2, sd), 20)

gc()
dim(E)
dim(Mp)

sel <- 1:nrow(Mp);

## At least 10 examples so vectors aren't too noisy

commonMuts <-
    intersect(sel,
              setdiff(1:nrow(Mp),
                      c(grep("=NaN", rownames(Mp)),
                        grep("=NA", rownames(Mp)),
                        grep("=0", rownames(Mp)))))
rownames(Mp)[commonMuts][1:5]
CNAfreq <- apply(Mp[commonMuts,], 1, sum);
hist((CNAfreq), 20)

CNAfreq["PTEN=-1"]; CNAfreq["PTEN=1"]; CNAfreq["TP53=-1"]; CNAfreq["TP53=1"]
## minFrequency <- .1 * ncol(Mp); #about 5% of samples
## minFrequency <- 50;
minFrequency <- -1
abline(v=(minFrequency), lty=2)

sum(apply(Mp[commonMuts,], 1, sum) >= minFrequency) / nrow(Mp)
commonMuts <-
    intersect(commonMuts,
              which(apply(Mp, 1, sum) >= minFrequency))
length(commonMuts)
rownames(Mp)[commonMuts][1:5]

MpP <- apply(Mp[commonMuts,], 1, function(x) { x / sum(x) })

if ( altType == "CNAs" ) {
    ctrlMuts <- gsub('=.*$', '=0', rownames(Mp)[commonMuts])
    MpN <- apply(Mp[ctrlMuts,], 1, function(x) { x / sum(x) })
} else {
    MpN <- apply(1 - Mp[commonMuts,], 1, function(x) { x / sum(x) })
}

Xd <- E %*% ( MpP - MpN )

## Estimating mutation effects from linear regression
## Xd <-
##     t(
##     (solve( Mp %*% t(Mp) + 50 * diag(1, nrow=dim(Mp)[1]) ) ) %*%
##     Mp %*% t(E)
##     )[,commonMuts]
## colnames(Xd) <- rownames(Mp)[commonMuts]

dim(Xd)
Xd[1:5,1:5]

dim(Mp[commonMuts,])
X <- Xd;
rm(Xd, MpP, MpN); gc();

## X <- X[,commonMuts]

## These are the mutations with strongest impact on gene expression
rev(sort(apply(X, 2, sd)))[1:10]

library(ade4)
## dudi1 <- dudi.pca(t(X), scale=F, center=F, scannf=F, nf=3)
## cumsum(dudi1$eig)[1:5] / sum(dudi1$eig)
## library(rgl)
## plot3d(dudi1$li)

mutWithStrongestEffect <- names(rev(sort(apply(X, 2, sd)))[1:10])
## mutWithStrongestEffect <-
##     names(rev(sort(apply(dudi1$li[,2:3]^2, 1, sum)))[1:8])
## mutWithStrongestEffect <-
##     names(rev(sort(apply(dudi1$li[,2:3]^2, 1, sum)))[1:2])
## dudi1$li[mutWithStrongestEffect,2:3]

## ## first axis is healthy -> cancer axis
## pdf("mutShiftsPCA.pdf", height=4, width=4)
## plot(dudi1$li[,2:3],
##      xlim=range(dudi1$li[,2]) * 1.1,
##      ylim=range(dudi1$li[,3]) * 1.1)
## ## text(dudi1$li[mutWithStrongestEffect,2:3],
## ##      sub("=1", "", mutWithStrongestEffect), pos=3)
## s.arrow(dudi1$li[mutWithStrongestEffect,2:3],
##         label=sub("=1$", "", mutWithStrongestEffect), add.plot=T)
## dev.off();
## rm(dudi1)

dudi2 <- dudi.pca(t(E), scale=F, scannf=F, nf=3)

## Xproj <- t(X) %*% as.matrix(dudi2$co) %*%
##     diag(1 / apply(dudi2$co, 2, function(x) { (sum(x^2)) }))

avgE <- apply(E, 1, mean)
## Compute the end points of mutations
## X <- cbind(X, CTRL=avgE)

save.image(sprintf("learnMutEffects2_%s.rda", altType))
## load(sprintf("learnMutEffects2_%s.rda", altType))

projMat <- as.matrix(dudi2$c1);
EPCAeig <- dudi2$eig;
rm(dudi2); gc()
## save(X, avgE, projMat, file="pre.rda");
## load("pre.rda")
X0 <- apply(X, 2, function(x) { x - avgE })
Xproj <- t(X0) %*% projMat;
save(Xproj, file="Xproj.rda")
## Xproj <-
##     t(apply(Xd, 2, function(x) { x - avgE })) %*%
##     as.matrix(dudi2$c1)

## Norms of the PCs:
## apply(dudi2$c1, 2, function(x) { sqrt(sum(x^2)) }) ## 1, good
Xproj[mutWithStrongestEffect,]
## XbackProj <- as.matrix(dudi2$co) %*% t(Xproj)
## tmp <- as.matrix(dudi2$c1) %*% t(Xproj)
## XbackProj <- apply(tmp, 2, function(x) { x + avgE })
XbackProj <- projMat %*% t( t(X) %*% projMat )
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

rm(XbackProj)

## XbackProj <- projMat %*% t( t(X) %*% projMat )
## res <- sapply(colnames(X), function(i) {
##     ##residuals: SS off Pareto front
##     sum(((X[,i] - mean(X[,i])) - XbackProj[,i])^2)
## })
## onf <- sapply(colnames(X), function(i) {
##     #SS on Parto front
##     sum(XbackProj[,i]^2)
## })
## ## tot <- sapply(colnames(X), function(i) {
## ##     ##residuals: SS off Pareto front
## ##     sum((X[,i] - mean(X[,i]))^2)
## ## })
## ## summary(res + onf - tot)
## rev(sort(onf/res))[1:50]
## ratios <- log10(onf/res);
## sum(dudi2$eig[1:3])/sum(dudi2$eig)
## ratios <- onf/(onf+res)
## hist(ratios)
## abline(v=sum(dudi2$eig[1:3])/sum(dudi2$eig), lty=2)
## rm(onf); rm(res);

arcsOrig <-
    read.table("arcsOrig_genes_tabSep.tsv", sep="\t", h=F, as.is=T)[,-1]

## Randomize mutations and see what the onf/res ratio is for random
## mutations

save.image(sprintf("learnMutEffects3_%s.rda", altType))
## load(sprintf("learnMutEffects3_%s.rda", altType))

## 3D "Pareto front" subspace
subspaces <- list()
subspaces[["Front"]] <- projMat;

## healthy -> cancer axis
tmp <- as.numeric(apply(arcsOrig[-healthyArch,], 2, mean) - arcsOrig[healthyArch,])
tmp <- tmp / sqrt(sum(tmp^2))
sum(tmp^2)
summary(tmp)
summary(projMat)
subspace <- matrix(tmp, ncol=1)
subspaces[["HealthyCancer"]] <- subspace;

## cancer archetypes subspace
is <- setdiff(1:nrow(arcsOrig), healthyArch)
tmp <- cbind(as.numeric(arcsOrig[is[2],] - arcsOrig[is[1],]),
             as.numeric(arcsOrig[is[3],] - arcsOrig[is[1],]))
tmp[,1] <- tmp[,1] / sqrt(sum(tmp[,1]^2))
tmp[,2] <- tmp[,2] - tmp[,1] * sum(tmp[,1] * tmp[,2])
tmp[,2] <- tmp[,2] / sqrt(sum(tmp[,2]^2))
tmp[1:5,]
summary(tmp)
sum(tmp[,1] * tmp[,2])
subspace <- tmp
subspaces[["CancerPlane"]] <- subspace;

subspace <- subspaces[["Front"]]

nPoss <-
    round(exp(seq(log(10), log(max(apply(Mp[commonMuts,], 1, sum))), len=12)))
i <- 1

save(E, nPoss, subspaces, file=sprintf("SSonDistInput_%s.rda", altType))

## Run SSonDist.R (in a different R session for memory efficiency to compute SSonDistRevCumApprox.rda
load(sprintf("SSonDistRevCumApprox_%s.rda", altType))
## load(sprintf("SSonDistRevCumApprox_%s.rda", "CNAs"))

rm(subspace)
names(subspaces)
## sIdx <- 1
## for (sIdx in 1:length(subspaces)) {

drivers <- list()
drivers[["Santarius2010"]] <-
    read.csv("/home/jhausser/work/cancerTaskAtlas/drivers/Santarius2010.csv",
             as.is=T, h=F)[,1]
drivers[["Pereira2016"]] <-
    read.csv("/home/jhausser/work/cancerTaskAtlas/drivers/Pereira2016.list",
             as.is=T, h=F)[,1]
drivers[["NikZainal2016"]] <-
    read.csv("/home/jhausser/work/cancerTaskAtlas/drivers/NikZainal2016-S14.Type.csv",
             as.is=T, h=T, skip=1)[,-1]
drivers[["NikZainal2016"]] <- unique(drivers[["NikZainal2016"]][,"Gene"])

head(drivers[["NikZainal2016"]])
sapply(drivers, length)
overlap <- function(A,B) {
    list(both=intersect(A,B),
         onlyA=setdiff(A,B),
         onlyB=setdiff(B,A))
}

library(tidyverse)
geneChr <- read_tsv("../hs_gene_coordinates_biomart.tsv", col_names=T)[,-1]
geneChr <- geneChr[,-ncol(geneChr)]
geneChr <- geneChr[,ncol(geneChr):1] %>% filter(grepl("^[0-9]*$", chr))
## geneChr[,"chr"] <- as.numeric(geneChr[,"chr"])
chrLens <- group_by(geneChr, chr) %>% summarize(chrLen=max(end)) %>% arrange(chr)
chrOffsets <- c(0, cumsum(as.numeric(unlist(chrLens[,"chrLen"]))))
chrOffsets <- chrOffsets[-length(chrOffsets)]
names(chrOffsets) <- unlist(chrLens[,"chr"])
genomeLen <- sum(as.numeric(unlist(chrLens[,"chrLen"])))

universeGeneNames <- unique(sub("=.*$", "", rownames(Mp)))
qValThreshold <- .01;

clusterGenomCoords <- function(coords, expand=1) {
    ordIdcs <- order(coords)
    clustIdx <- rep(-1, length(coords))
    clustIdx[1] <- 0;
    i <- 2
    for (i in 2:length(ordIdcs)) {
        if ( abs(coords[ordIdcs[i]] - coords[ordIdcs[i-1]]) <
             2 * expand ) {
            clustIdx[i] <- clustIdx[i-1]
        } else {
            clustIdx[i] <- clustIdx[i-1] + 1
        }
    }
    oriClustIdx <- numeric(length(clustIdx))
    oriClustIdx[ordIdcs] <- clustIdx
    return(oriClustIdx)
}


library(parallel)
s <- names(subspaces)[1]
## s <- names(subspaces)[2]
mutsAlongFront <- list()
for (s in names(subspaces)) {
    cat(sprintf("%s\n", s))
    subspace <- subspaces[[s]]
    
    XbackProj <- as.matrix(subspace) %*% t( t(X) %*% as.matrix(subspace) )
    res <- sapply(colnames(X), function(i) {
        ##residuals: SS off Pareto front
        sum((X[,i] - XbackProj[,i])^2)
    })
    onf <- sapply(colnames(X), function(i) {
                                        #SS on Parto front
        sum(XbackProj[,i]^2)
    })
    ratios <- onf/(onf+res)
    norms <- sqrt(onf+res)
    nPos <- apply(Mp[commonMuts,], 1, sum)
    
    SSlist <-
        mclapply(1:5, function(i) {
            cat(sprintf("Iteration %d\n", i))

            ## summary(apply(Mp[commonMuts,], 1, sum))
            Mr <- t(apply(Mp[commonMuts,], 1, sample))
            ## Xr <-
            ##     E %*% t(Mr) %*% diag(1 / apply(Mr, 1, sum)) -
            ##     E %*% t(1-Mr) %*% diag(1 / apply(1-Mr, 1, sum))
            MrP <- apply(Mr, 1, function(x) { x / sum(x) })
            MrN <- apply(1-Mr, 1, function(x) { x / sum(x) })
            Xr <- E %*% ( MrP - MrN )
            colnames(Xr) <- rownames(Mr)
            
            XrbackProj <- as.matrix(subspace) %*% t( t(Xr) %*% as.matrix(subspace) )
            
            resr <- sapply(colnames(X), function(i) {
                ##residuals: SS off Pareto front
                sum(((Xr[,i] - mean(Xr[,i])) - XrbackProj[,i])^2)
            })
            onfr <- sapply(colnames(X), function(i) {
                ##SS on Parto front
                sum(XrbackProj[,i]^2)
            })
            return(list(resr=resr, onfr=onfr,
                        nPos=apply(Mp[commonMuts,], 1, sum)))
        })

    onfr <- (sapply(SSlist, function(s) { s[["onfr"]] }))
    resr <- (sapply(SSlist, function(s) { s[["resr"]] }))
    nPosr <- (sapply(SSlist, function(s) { s[["nPos"]] }))
    rratios <- onfr/(onfr+resr);
    rnorms <- sqrt(onfr + resr)
    rratios[1:5,1:5]
    rnorms[1:5,1:5]

    x <- Mp[commonMuts,][1,]
    ## x <- Mp[commonMuts,][1082,]
    nPosIdx <- apply(Mp[commonMuts,], 1, function(x) {
        which( sum(x) <= nPoss )[1] - 1
    })

    ## plot(apply(Mp[commonMuts,], 1, sum), nPoss[nPosIdx]); abline(0,1)
    i <- 1
    ## i <- 1082
    log10pSS <- sapply(1:length(norms), function(i) {
        ## sum(pSSonParms[[s]][,nPosIdx[i]] * c(1, onf[i]))
        sum(pSSonParms[["SS"]][,nPosIdx[i]] * c(1, norms[i]^2))
    })
    log10pRatios <- sapply(1:length(ratios), function(i) {
        x <- ratios[i]
        n <- length(pSSonParms[[s]])
        sum(pSSonParms[[s]] * x^seq(0,n-1))
    })
    summary(log10pSS)
    summary(log10pRatios)
    ## hist(10^log10pRatios)
    ## Bonferroni correction
    ## log10q <- log10p + log10(length(onf))
    ## log10q[log10q>0] <- 0
    ## hist(log10q)
    sum(10^log10pSS < qValThreshold/length(norms))
    sum(10^log10pRatios < qValThreshold/length(ratios))
    
    ## x <- Mp[commonMuts,][1,]
    ## nPosrIdx <- sapply(nPosr, function(x) {
    ##     which( x < nPoss )[1] - 1
    ## })
    ## plot(nPosr, nPoss[nPosrIdx]); abline(0,1)
    ## i <- 1
    ## log10pr <- sapply(1:length(onfr), function(i) {
    ##     sum(pSSonParms[[s]][,nPosrIdx[i]] * c(1, onfr[i]))
    ## })
    ## log10pr[log10pr > 0] <- 0
    ## summary(log10pr)
    ## hist(10^log10pr)
    ## histogram of p-values looks flat!
    ## It's just that many of the mutations are well aligned with the
    ## Pareto front
    
    ## Let's determine p-values by comparing each mutation to its own
    ## randomized versions
    ## sigmas <-
    ##     sapply(names(ratios), function(m) {
    ##         ( ratios[m] - mean(rratios[m,]) ) / sd(rratios[m,])
    ##     })
    
    ## names(sigmas) <- names(ratios)
    ## pVals <- pnorm(- abs(sigmas)) * 2;
    pVals <- 10^log10pRatios;
    names(pVals) <- rownames(Mp[commonMuts,])
    
    ##     hist(pVals)
    ##     library(fdrtool)
    ## qVals <- fdrtool(pVals, statistic="pvalue")$q
    ## plot(sigmas, log10(qVals))
    qVals <- pVals * length(pVals)
    qVals[qVals>1] <- 1
    ## hist(qVals)
    sort(qVals)[1:10]

    mybreaks <- seq(min(c(ratios, rratios)), max(c(ratios, rratios)), len=20);
    pdf(sprintf("fracSSonVsOff_%s_%s.pdf", altType, s), height=4, width=4)
    hist(ratios,
         breaks=mybreaks, main="",
         col="lightblue", freq=F, ## ylim=c(0,5),
         xlim=c(0,1),
         ylim=c(0, max(table(cut(onfr/(onfr+resr), breaks=mybreaks)) /
                       (length(onfr)) / diff(mybreaks)[1])),
         ## xlab=expression(paste(log[10], " (SSon / SSoff)")),
         xlab=sprintf("frac of SS on %s", s))
    hr <- hist(onfr/(onfr+resr), breaks=mybreaks, border="red",
               col=NULL, add=T, freq=F)
    if ( s == "Front" ) {
        abline(v=sum(EPCAeig[1:3])/sum(EPCAeig), lty=2)
    }
    
    legend("topright", c("measured", "shuffled"),
           fill=c("lightblue", "white"), col=c("black", "red"),
           bg="white")
    ## ## sapply(names(pVals[pVals<.01 | pVals>.99]), function(x) {
    ## sapply(names(sigmas[qVals<.01 & sigmas > 0]), function(x) {
    ## ## sapply(c("PVT1=1_CNA", "KARS=-1_CNA"), function(x) {
    ##     text(ratios[x], .3, sub("=1", "", x), srt=90)
    ##     arrows(ratios[x], .13, ratios[x], .02, len=.07)
    ## })
    dev.off();

    allMuts <- data.frame(
        ratio=ratios,
        norm=norms,
        q=qVals,
        frequency=apply(Mp[commonMuts,], 1, sum),
        log10pNorm=log10pSS,
        log10pRatio=log10pRatios)
    mutsAlongFront[[s]] <- allMuts[qVals<qValThreshold,] 
    
    ## test for match between drivers and found mutations
    tDrivers <- drivers[[1]]
    driversOverlap <-
        sapply(drivers, function(tDrivers) {
            ## Only look at drivers that are in the universe
            ## tDrivers <- tDrivers[tDrivers %in% universeGeneNames]
            
            ##total white balls
            m <- length(unique(sub("=.*$", "",
                                   rownames(mutsAlongFront[[s]]))))
            
            ##total black balls
            n <- length(universeGeneNames) - m
            n <- 25000 - m
            
            ##number of drawn balls
            k <- length(tDrivers)

            ##number of while balls in the draw
            q <- sum(tDrivers %in%
                     unique(sub("=.*$", "",
                                rownames(mutsAlongFront[[s]]))))
            c(sapply(overlap(tDrivers,
                             unique(sub("=.*$", "",
                                        rownames(mutsAlongFront[[s]])))),
                     length),
              p=phyper(q=q, m=m, n=n, k=k, lower.tail=F))
        })

    ## Add 'driver' indicator to table
    toAdd <- names(which(driversOverlap["p",] < .01))
    if ( length(toAdd) > 0 ) {
        x <- toAdd[1]
        for ( x in toAdd ) {
            tmp <- cbind(mutsAlongFront[[s]],
                         sub("=.*$", "", rownames(mutsAlongFront[[s]])) %in%
                         drivers[[x]])
            colnames(tmp)[ncol(tmp)] <- sprintf("%s (driver)", x)
            mutsAlongFront[[s]] <- tmp
        }
    }
    myDrivers <- unique(unlist(sapply(toAdd, function(x) { drivers[[x]] })))

    library(tidyverse)
    allMuts[,"CNA"] <- rownames(allMuts)
    library(stringr)
    allMuts <- as_tibble(allMuts) %>%
        mutate(name=str_replace(CNA, "=.*$", ""))
    allMutsChr <- left_join(allMuts, geneChr)
    allMutsChr <-
        left_join(allMutsChr, tibble(chr=names(chrOffsets),
                                     offset=chrOffsets)) %>%
        filter(!is.na(chr)) %>%
        mutate(pos=start+offset)
    driversTibble <- tibble(name=myDrivers, driver=T)
    allMutsChrD <-
        left_join(allMutsChr, driversTibble) %>%
        mutate(driver=!is.na(driver))

    mychr <- "1"
    for (mychr in unique(unlist(allMutsChrD[,"chr"]))) {
        gIdcs <- which(allMutsChrD[,"chr"] == mychr)
        myClusts <- clusterGenomCoords(unlist(allMutsChrD[gIdcs,"start"]), expand=10e6)
        allMutsChrD[gIdcs,"chrCluster"] <- paste(myClusts, mychr, sep="-")
    }

    allMutsChrD[allMutsChrD[,"chrCluster"] == "16-0",
                c("chr", "start")] %>% arrange(start) %>% print(n=100)
    isParallel <- allMutsChrD[,"log10pRatio"] <= -5;
    isDriver <- unlist(allMutsChrD[,"driver"])
    with(allMutsChrD,
         plot(pos[isParallel], frequency[isParallel],
              col=as.numeric(factor(chrCluster))[isParallel], ylim=c(0,800)))
    with(allMutsChrD,
         points(pos[isDriver], frequency[isDriver],
                col=as.numeric(factor(chrCluster))[isDriver], pch=3))

    library(ggrepel)
    ggplot(filter(allMutsChrD, driver | log10pRatio<pValThreshold)) +
        geom_point(mapping=aes(x=pos, y=-log10pRatio, shape=driver,
                               col=factor(chrCluster)), show.legend=F) +
        geom_vline(tibble(boundaries=chrOffsets[-1]),
                   mapping=aes(xintercept=boundaries), linetype=2) +
        geom_hline(mapping=aes(yintercept=-pValThreshold)) +
        geom_label_repel(mapping=aes(x=pos, y=-log10pRatio, label=CNA),
                         data=filter(allMutsChrD, driver & log10pRatio<pValThreshold)) +
        geom_label(mapping=aes(x=offset, y=1/2, label=chr),
                   data=tibble(chr=names(chrOffsets),
                               offset=chrOffsets),
                   nudge_x=2.5e7) +
        coord_cartesian(xlim=c(0,2.8e9)) +
        labs(x="Genomic position", y=quote(paste(log[10], " ", p))) +
        scale_x_continuous(breaks=scales::pretty_breaks(n=10))
    ggsave(sprintf("chrAlnFrontDrivers_%s_%s.pdf", altType, s), height=3, width=12)

    pValThreshold <- log10(qValThreshold / length(commonMuts))
    ggplot(filter(allMutsChrD, driver | log10pRatio<pValThreshold)) +
        geom_point(mapping=aes(x=pos, y=-log10pRatio, color=driver)) +
        geom_vline(tibble(boundaries=chrOffsets[-1]),
                   mapping=aes(xintercept=boundaries), linetype=2) +
        geom_hline(mapping=aes(yintercept=-pValThreshold)) +
        geom_label_repel(mapping=aes(x=pos, y=-log10pRatio, label=CNA),
                         data=filter(allMutsChrD, driver & log10pRatio<pValThreshold)) +
        geom_label(mapping=aes(x=offset, y=1/2, label=chr),
                   data=tibble(chr=names(chrOffsets),
                               offset=chrOffsets),
                   nudge_x=2.5e7) +
        coord_cartesian(xlim=c(0,2.8e9)) +
        labs(x="Genomic position", y=quote(paste(log[10], " ", p))) +
        scale_x_continuous(breaks=scales::pretty_breaks(n=10))
    ggsave(sprintf("chrAlnFrontDriversAlt_%s_%s.pdf", altType, s), height=3, width=12)
    
    ## distribution of distances to nearest driver
    g <- "NUP93"
    g <- "FBXO31"

    getDistToNextDriver <- function(g) {
        ## cat(sprintf("%s\n", g))
        gPos <-
            mean(as.numeric(unlist(unique(
                allMutsChrD[allMutsChrD[,"name"] == g,
                            "start"]))))
        gChr <-
            as.character(unique(
                allMutsChrD[allMutsChrD[,"name"] == g,
                            "chr"]))
        driversOnChr <- allMutsChrD %>%
            filter(chr == gChr, driver)
        if ( nrow(driversOnChr) == 0 ) {
            return(as.numeric(chrLens[chrLens[,"chr"] == gChr,"chrLen"]))
        } else {
            return(min(abs(gPos -
                           unique(driversOnChr[,"start"]))))
        }
    }

    rndPosDists <-
        sapply(1:1000, function(i) {
            gChr <- sample(unlist(chrLens[,"chr"]), size=1)
            gPos <- runif(1, min=0,
                          max=as.numeric(chrLens[chrLens[,"chr"] == gChr, "chrLen"]))
            driversOnChr <- allMutsChrD %>%
                filter(chr == gChr, driver)
            if ( nrow(driversOnChr) == 0 ) {
                return(as.numeric(chrLens[chrLens[,"chr"] == gChr,"chrLen"]))
            } else {
                return(min(abs(gPos -
                               unique(driversOnChr[,"start"]))))
            }
        })
    
    alnGenes <- unique(sub("=.*$", "", rownames(mutsAlongFront[[s]])))
    alnGenesDists <- sapply(alnGenes, getDistToNextDriver)
    
    allGenes <- unique(unlist(allMuts[,"name"]))
    allGenesDists <- sapply(allGenes, getDistToNextDriver)

    log10(c(median(alnGenesDists), median(allGenesDists, na.rm=T), median(rndPosDists)))

    pdf(sprintf("distToNearestDriver_%s_%s.pdf", altType, s),
        height=4.5, width=4.5)
    plot(ecdf(allGenesDists), col="black", main="",
         xlab="distance to closest driver",
        ylab="empirical cumulative distribution")
    lines(ecdf(rndPosDists), col="grey", lwd=2)
    lines(ecdf(alnGenesDists), col="blue")
    text(1e8, .45,
         sprintf("p = %.1e",
                 wilcox.test(log10(alnGenesDists), log10(rndPosDists))$p.value))
    legend("bottomright",
           sprintf("%d %s",
                   c(length(alnGenes), length(allGenes), length(rndPosDists)),
                   c("genes parallel to front", "cancer genes", "random locations")),
           col=c("blue", "black", "grey"), lwd=2, bg="white")
    dev.off();
    
    ## hist(log10(alnGenesDists), breaks=seq(3.5, 9, by=.5), freq=F,
    ##      col="lightblue", main="",
    ##      xlab=expression(paste(log[10], " distance to closest driver [nt]")))
    ## hist(log10(allGenesDists), breaks=seq(3.5, 9, by=.5), freq=F,
    ##      add=T)
    ## abline(v=log10(median(alnGenesDists, na.rm=T)), lty=2, col="grey")
    ## ## hist(log10(rndPosDists), breaks=seq(3.5, 9, by=.5), freq=F, add=T,
    ## ##      border="green")
    ## abline(v=log10(median(rndPosDists, na.rm=T)), lty=2, col="green")
    ## legend("topleft",
    ##        sprintf("%d %s",
    ##                c(length(alnGenes), length(allGenes)),
    ##                c("genes parallel to front", "cancer genes")),
    ##        fill=c("lightblue", "white"))
    ## text(8.5, .9,
    ##      sprintf("p = %.4f",
    ##              wilcox.test(log10(alnGenesDists), log10(allGenesDists))$p.value))
    
    write.csv(mutsAlongFront[[s]][order(mutsAlongFront[[s]][,"q"]),],
              file=sprintf("onFront_%s_%s.csv", altType, s))

    pdf(sprintf("fracSSonFront_SS_freq_%s_%s.pdf", altType, s), height=4, width=4*3)
    par(mfrow=c(1,3))
    plot(rratios, rnorms^2, col="grey",
         xlab=sprintf("fraction of ||v||^2 on %s", s),
         ylab="effect size ||v||^2", xlim=c(0,1), log="y")
    with(allMutsChrD,
         points(ratio, norm^2## , col=factor(chrCluster)
                ))
    myp <- 1e-6
    sapply(10^c(-6, -5, -4, -3, -2), function(myp) {
        res <- optimize(function(r) {
            abs(sum(pSSonParms[[s]] * r^seq(0, length(pSSonParms[[s]])-1)) - log10(myp))
        }, c(0,1))
        abline(v=res$minimum, lty=2, col="grey")
        text(res$minimum, max(rnorms^2), sprintf("%.1e", myp), srt=90, pos=2)
    })
    sapply(myDrivers, function(g) {
        sel <- sub("=.*$", "", names(ratios)) == g
        points(ratios[sel], norms[sel]^2, pch=20, col="red")
    })
    
    plot(rratios, nPosr, col="grey",
         xlab=sprintf("fraction of ||v||^2 on %s", s), ylab="frequency",
         xlim=c(0,1))
    with(allMutsChrD,
         points(ratio, frequency## , col=as.factor(unlist(allMutsChrD[,"chrCluster"]))
           ))
    myp <- 1e-6
    sapply(10^c(-6, -5, -4, -3, -2), function(myp) {
        res <- optimize(function(r) {
            abs(sum(pSSonParms[[s]] * r^seq(0, length(pSSonParms[[s]])-1)) - log10(myp))
        }, c(0,1))
        abline(v=res$minimum, lty=2, col="grey")
        text(res$minimum, max(nPosr), sprintf("%.1e", myp), srt=90, pos=2)
    })
    sapply(myDrivers, function(g) {
        sel <- sub("=.*$", "", names(ratios)) == g
        points(ratios[sel], nPos[sel], pch=20, col="red")
    })

    plot(rnorms^2, nPosr, col="grey",
         xlab="effect size ||v||^2", ylab="frequency", log="x")
    points(norms^2, nPos## , col=as.factor(unlist(allMutsChrD[,"chrCluster"]))
           )
    myp <- 1e-6
    sapply(10^c(-6, -5, -4, -3, -2), function(myp) {
        i <- 1
        xs <-
            sapply(1:length(nPoss), function(i) {
                (log10(myp) - pSSonParms[["SS"]]["a",i]) /
                    pSSonParms[["SS"]]["b",i]
            })
        lines(xs, nPoss, lty=2, col="grey")
        if ( myp < 5e-6 || myp > 10^(-2.5) ) {
            text(min(xs) * 0.8 + as.numeric(myp < 5e-6) * min(xs) * 0.5,
                 max(nPosr) * 0.95, sprintf("%.1e", myp),
                 srt=90)
        }
    })
    sapply(myDrivers, function(g) {
        sel <- sub("=.*$", "", names(norms)) == g
        points(norms[sel]^2, nPos[sel], pch=20, col="red")
    })
    legend("topright",
           c("CNAs", "random controls",
             "drivers"),
           ## sprintf("drivers (%s)", paste(toAdd, collapse=", "))),
           pch=c(1, 1, 20), col=c("black", "grey", "red"))
    dev.off();

    library(rgl)
    plot3d(rratios, log10(rnorms^2), nPosr, col="grey", alpha=0,
           xlab="frac of ||v||^2 on front", ylab="log10(||v||^2)", zlab="frequency")
    with(allMutsChrD,
         points3d(ratio, log10(norm^2), frequency,
                  col=as.numeric(factor(chrCluster))+1,
                  alpha=.8))

    ## driversCloud <- sapply(myDrivers, function(g) {
    ##     sel <- sub("=.*$", "", names(norms)) == g
    ##     c(ratio=ratios[sel], log10norm=log10(norms[sel]^2), freq=nPos[sel])
    ## })
    ## points3d(driversCloud[,1], driversCloud[,2], driversCloud[,3], col="red")
    
    ## if ( s == "Front" ) {
    ##     save(nPosr, rratios, file="nPosr_rratios.rda")
    ##     write.csv(allMuts[qVals==1,][rev(order(allMuts[qVals==1,"frequency"])),],
    ##               file=sprintf("offFront_%s.csv", altType))
    ## }

    if ( nrow(mutsAlongFront[[s]]) > 1 ) {
        mutsAlongFront[[s]] <-
            mutsAlongFront[[s]][order(mutsAlongFront[[s]][,"q"]),]
    }
    ## mutWithStrongestEffect <- names(sigmas[qVals<.1 & sigmas > 0])
    ## mutWithStrongestEffect <- names(sigmas[qVals<.05 & sigmas > 0])
}

sapply(mutsAlongFront, nrow)

sapply(overlap(drivers[["Santarius2010"]],
               drivers[["Pereira2016"]]), length)
sapply(overlap(drivers[["Santarius2010"]],
               drivers[["NikZainal2016"]]), length)
sapply(overlap(drivers[["Pereira2016"]],
               drivers[["NikZainal2016"]]), length)

if ( altType == "CNAs" ) {
    universeGeneNames <- unique(sub("=.*$", "", copGeneNames))
} else {
    universeGeneNames <- unique(sub("=.*$", "", rownames(Mp)))
}


library(topGO)

GOtab <- lapply(names(mutsAlongFront), function(s) {
    geneList <- factor(as.integer(
        universeGeneNames %in%
        unique(gsub("=.*$", "", rownames(mutsAlongFront[[s]])))
    ))
    names(geneList) <- universeGeneNames;
    
    GOdataBP <- new("topGOdata",
                    description = s, ontology = "BP",
                    allGenes = geneList, nodeSize = 10,
                    annot = annFUN.org,
                    mapping="org.Hs.eg.db",
                    ID="symbol")
    resultDef <- 
        runTest(GOdataBP, algorithm = "weight01", statistic = "fisher")
    allResBP <- GenTable(GOdataBP, Weighted=resultDef, topNodes=100)
    
    GOdataMF <- new("topGOdata",
                    description = s, ontology = "MF",
                    allGenes = geneList, nodeSize = 10,
                    annot = annFUN.org,
                    mapping="org.Hs.eg.db",
                    ID="symbol")
    resultDef <- 
        runTest(GOdataMF, algorithm = "weight01", statistic = "fisher")
    allResMF <- GenTable(GOdataMF, Weighted=resultDef, topNodes=100)

    allRes <- rbind(allResBP, allResMF)
    allRes <- allRes[order(as.numeric(allRes[,"Weighted"])),]
    head(allRes)
    ## hist(as.numeric(allRes[,"Weighted"]))
    return(allRes)
})

lapply(GOtab, head)

## For SNVs, too little mutations to make a difference.

## No big differences between enriched GO categories in the three sets
## of CNAs:
## poly(A) RNA binding
## mitochondrial translation elongation & termination
## mitosis

## Poly(A) mRNA binding
intersect(annFUN.org(whichOnto="MF", mapping="org.Hs.eg.db",
                     ID="symbol")[["GO:0044822"]],
          unique(gsub("=.*$", "", rownames(mutsAlongFront[["Front"]]))))
## Translation regulation, splicing / NMD factors, DNA repair




head(mutsAlongFront[["Front"]])
plot(mutsAlongFront[["Front"]][,"ratio"],
     sqrt(mutsAlongFront[["Front"]][,"SSon"]))
plot(mutsAlongFront[["Front"]][,"ratio"],
     1 / sqrt(mutsAlongFront[["Front"]][,"ratio"] / mutsAlongFront[["Front"]][,"SSon"]))
sapply(seq(10, 50, by=10), function(offset) {
    lines(seq(0, 1, len=100),
          offset / seq(0, 1, len=100), lty=2, col="grey")
    })
rownames(mutsAlongFront[["Front"]])[
    identify(mutsAlongFront[["Front"]][,"ratio"],
             1 / sqrt(mutsAlongFront[["Front"]][,"ratio"] /
                      mutsAlongFront[["Front"]][,"SSon"]),
             rownames(mutsAlongFront[["Front"]]))
]
## All the genes in the top right quadrant have something to do with
## transcription control and the cell cycle!


    

## mutations with large effects on the front are usually aligned with
## the front
library(rgl)
plot3d(mutsAlongFront[["Front"]][,"ratio"],
       1 / sqrt(mutsAlongFront[["Front"]][,"ratio"] / mutsAlongFront[["Front"]][,"SSon"]),
       mutsAlongFront[["Front"]][,"frequency"],
       xlab="on/(on+off)",
       ylab="mutation effect ||v||",
       zlab="frequency")
for (i in grep("(driver)", colnames(mutsAlongFront[["Front"]]))) {
    sel <- mutsAlongFront[["Front"]][,i]
    spheres3d(
        mutsAlongFront[["Front"]][sel,"ratio"],
        1 / sqrt(mutsAlongFront[["Front"]][sel,"ratio"] /
                 mutsAlongFront[["Front"]][sel,"SSon"]),
        mutsAlongFront[["Front"]][sel,"frequency"],
        col="red", radius=5)
}


plot(sqrt(mutsAlongFront[["Front"]][,"SSon"]),
     mutsAlongFront[["Front"]][,"frequency"],
     xlab="||v||", ylab="frequency")
identify(sqrt(mutsAlongFront[["Front"]][,"SSon"]),
         mutsAlongFront[["Front"]][,"frequency"],
         rownames(mutsAlongFront[["Front"]]))

pSSonParms[["Front"]]
nPoss

load("nPosr_rratios.rda")
range(nPosr)
rratios

pm70 <- (abs(apply(nPosr, 1, unique) - 70) < 7)
pm250 <- (abs(apply(nPosr, 1, unique) - 250) < 25)
pm500 <- (abs(apply(nPosr, 1, unique) - 500) < 50)
prod(dim(rratios[pm500,]))
prod(dim(rratios[pm250,]))
prod(dim(rratios[pm70,]))

library(Hmisc)

myecdf <- Ecdf(rratios[pm70,], what="1-F", log="",
               xlab="ratios", ylab="1-F", lwd=2, col="blue");
myecdf <- Ecdf(rratios[pm250,], what="1-F", log="",
               xlab="ratios", ylab="1-F", lwd=2, col="grey", add=T);
myecdf <- Ecdf(rratios[pm500,], what="1-F", log="",
               xlab="ratios", ylab="1-F", lwd=2, col="black", add=T);

pdf(sprintf("frequency_fracOnFront_%s.pdf", altType),
    height=5, width=5)
plot(as.numeric(nPosr), as.numeric(rratios), col="grey",
     ylab="fraction on front", xlab="frequency",
     ylim=c(0,1), log="x")
points(mutsAlongFront[["Front"]][,"frequency"],
       mutsAlongFront[["Front"]][,"ratio"])
for (i in grep("(driver)", colnames(mutsAlongFront[["Front"]]))) {
    points(mutsAlongFront[["Front"]][mutsAlongFront[["Front"]][,i],"frequency"],
           mutsAlongFront[["Front"]][mutsAlongFront[["Front"]][,i],"ratio"],
           col="red", pch=20)

}
if ( nrow(mutsAlongFront[["Front"]]) < 10 ) {
    text(mutsAlongFront[["Front"]][,"frequency"],
         mutsAlongFront[["Front"]][,"ratio"],
         rownames(mutsAlongFront[["Front"]]), pos=3)
}
legend("topright",
       colnames(mutsAlongFront[["Front"]])[grep("(driver)",
                                                colnames(mutsAlongFront[["Front"]]))],
       pch=20, col="red")
dev.off();

identify(sqrt(mutsAlongFront[["Front"]][,"SSon"]),
         mutsAlongFront[["Front"]][,"frequency"],
         rownames(mutsAlongFront[["Front"]]))


intersect(rownames(mutsAlongFront[["HealthyCancer"]]),
          rownames(mutsAlongFront[["CancerPlane"]]))
## There are some mutations common to progression axis and cancer
## plane! These are the mutations generic to cancer:
mutsProg <- setdiff(rownames(mutsAlongFront[["HealthyCancer"]]),
                    rownames(mutsAlongFront[["CancerPlane"]]))
## These the mutations that distinguish cancers
mutsTypes <- setdiff(
    rownames(mutsAlongFront[["CancerPlane"]]),
    rownames(mutsAlongFront[["HealthyCancer"]]))

plot(as.numeric(arcsOrig[4,]), as.numeric(healthyProfile)); abline(0,1,col="grey")

## genesFilt <- read.table("geneNamesAfterExprFiltering.list", as.is=T)[,1]
## colnames(arcsOrig) <- genesFilt;

##################################################
## 3D plot of samples and mutations

apply(arcsOrig, 1, length)
archProj <- t(sapply(1:nrow(arcsOrig), function(i) {
    m0 <- arcsOrig[i,] - healthyProfile - avgE
    as.matrix(m0) %*% projMat
}))

hPt <-
    (- matrix(apply(E, 1, mean), nrow=1) ) %*% projMat 

s <- 10;
s <- 3;

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

## save.image(sprintf("plot3d_%s.rda", altType))
## load(sprintf("plot3d_%s.rda", altType))
library(rgl)
open3d()
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

if ( nrow(mutsAlongFront[["Front"]]) < 10 ) {
    ## sapply(rownames(mutsAlongFront[["HealthyCancer"]]), function(m) {
    sapply(rownames(mutsAlongFront[["Front"]]), function(m) {
        x <- Xproj[m,1:3]
        v <- x - hPt;
        
        lines3d(c(hPt[1], s * v[1] + hPt[1]),
                c(hPt[2], s * v[2] + hPt[2]),
                c(hPt[3], s * v[3] + hPt[3]), col="blue", lwd=2)
        text3d(s * v[1] + hPt[1], s * v[2] + hPt[2],
               s * v[3] + hPt[3], sub("_CNA$", "", sub("=1$", "", m)),
               adj=c(.5, 0))
    })
} else {
    ## Add mutations healthy -> cancer
    sapply(mutsProg, function(m) {
        x <- Xproj[m,1:3]
        v <- x - hPt;
        
        lines3d(c(hPt[1], s * v[1] + hPt[1]),
                c(hPt[2], s * v[2] + hPt[2]),
                c(hPt[3], s * v[3] + hPt[3]), col="black", lwd=2)
        ## text3d(s * v[1] + hPt[1], s * v[2] + hPt[2],
        ##        s * v[3] + hPt[3], sub("_CNA$", "", sub("=1$", "", m)), adj=c(.5, 0))
    })

    ## Add mutations parallel to the cancer plane
    sapply(mutsTypes, function(m) {
        x <- Xproj[m,1:3]
        v <- x - hPt;
        
        lines3d(c(hPt[1], s * v[1] + hPt[1]),
                c(hPt[2], s * v[2] + hPt[2]),
                c(hPt[3], s * v[3] + hPt[3]), col="red", lwd=2)
        ## text3d(s * v[1] + hPt[1], s * v[2] + hPt[2],
        ##        s * v[3] + hPt[3], sub("_CNA$", "", sub("=1$", "", m)), adj=c(.5, 0))
    })
}

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
