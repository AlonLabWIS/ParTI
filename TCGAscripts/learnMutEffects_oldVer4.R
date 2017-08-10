rm(list=ls()); gc();
source("~/work/jeanLib.R")
library(MASS)
library(rgl)
library(ade4)
library(parallel)
library(fdrtool)
library(topGO)
library(Hmisc)
library(cowplot)

library(tidyverse)
library(ggrepel)
library(stringr)

altType <- "CNAs"; #SNVs, CNAs
altType <- "SNVs"; #SNVs, CNAs
isMetabric <- length(grep("metabric", getwd())) == 1

##################################################

geneNames <- read.delim("geneListExp.list", sep="\t", h=F, as.is=T)[,1]
geneNamesFilt <- read.table("geneNamesAfterExprFiltering.list",
                            as.is=T, h=F, sep="\t")[,1]
geneExpression <- as.matrix(read_csv("expMatrix.csv", col_names=F))
## geneExpression <- read_csv("expMatrix.csv", col_names=F)
colnames(geneExpression) <- geneNames;
rm(geneNames);
tmp <- sapply(geneNamesFilt, function(g) {
    geneExpression[,g]
})
geneExpression <- tmp
geneExprAvg <- apply(geneExpression, 2, mean)

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
    isHealthy <- discClin[,"sample_type"] == "Solid Tissue Normal"
    healthyArch <-
        as.numeric(
            read_csv("clinicalEnrichment_discrete_significant.csv") %>%
            filter(`Feature Name` == 'sample type: Solid Tissue Normal') %>%
            select(`archetype #`)
        )
}
rm(discClin)

CNAs <- read_delim("copMatrix_reOrdered_booleanized.tsv", delim="\t", col_names=T)
CNAs <- as.data.frame(CNAs)
rownames(CNAs) <- CNAs[,1]
CNAs <- CNAs[,-1]

copGeneNames <- colnames(CNAs)
x <- copGeneNames[1]
copTab <- t(sapply(copGeneNames, function(x) {
    strsplit(x, "=")[[1]]
}))

healthyProfile <- apply(geneExpression[isHealthy,], 2, median)

gE0 <- apply(geneExpression, 1, function(x) { x - mean(x) })
rm(geneExpression)
gE00 <- apply(gE0, 2, function(x) { x - healthyProfile } )
## gE00 is centered on healthyProfile
rm(gE0); gc();

E <- gE00
rm(gE00, tmp)

## Mp <- rbind(M, t(CNAs))
if ( altType == "SNVs" ) {
    Mp <- t(mutations);
} else if ( altType == "CNAs" ) {
    Mp <- t(CNAs);
} else {
    stop("Unknown altType!");
}

rm(mutations, CNAs); gc();
save.image(sprintf("learnMutEffects_%s.rda", altType))

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
altFreq <- apply(Mp[commonMuts,], 1, function(x) { sum(x, na.rm=T) })
hist((altFreq), 20)

altFreq["PTEN=-1"]; altFreq["PTEN=1"]; altFreq["TP53=-1"]; altFreq["TP53=1"]
minFrequency <- -1
abline(v=(minFrequency), lty=2)

sum(apply(Mp[commonMuts,], 1, function(x) {sum(x, na.rm=T)}) >= minFrequency) / nrow(Mp)
commonMuts <-
    intersect(commonMuts,
              which(apply(Mp, 1, function(x) { sum(x, na.rm=T) }) >= minFrequency))
length(commonMuts)
rownames(Mp)[commonMuts][1:5]

## MpP <- apply(Mp[commonMuts,], 1, function(x) { x / sum(x, na.rm=T) })
MpP <- apply(Mp[commonMuts,], 1, function(x) { x / sum(x) })

## FIXME
## Controls are healthy samples! If we use all samples that don't have
## the mutations, we will partly assign the effect of a given mutation to
## other mutations. We use the simpler 'not mutated' control for now
## because any mutation may align with the front if healthy samples
## are the reference.
## MpN <-
##     sapply(1:length(commonMuts), function(i) {
##         isHealthy / sum(isHealthy)})
if ( altType == "CNAs" ) {
    ctrlMuts <- gsub('=.*$', '=0', rownames(Mp)[commonMuts])
    MpN <- apply(Mp[ctrlMuts,], 1, function(x) { x / sum(x) })
} else {
    MpN <- apply(1 - Mp[commonMuts,], 1, function(x) { x / sum(x, na.rm=T) })
}

Xd <- E %*% ( MpP - MpN )
Xd[1:5,1:5]
X <- Xd;
rm(Xd, MpP, MpN); gc();

## These are the mutations with strongest impact on gene expression
rev(sort(apply(X, 2, sd)))[1:10]

mutWithStrongestEffect <- names(rev(sort(apply(X, 2, sd)))[1:10])

dudi2 <- dudi.pca(t(E), scale=F, scannf=F, nf=3)

avgE <- apply(E, 1, mean)

save.image(sprintf("learnMutEffects2_%s.rda", altType))
## load(sprintf("learnMutEffects2_%s.rda", altType))

projMat <- as.matrix(dudi2$c1);
posInPCspace <- as.matrix(dudi2$li);
EPCAeig <- dudi2$eig;
rm(dudi2); gc()

## X is a vector. No need to rewrite it in order to project it!
## X0 <- apply(X, 2, function(x) { x - avgE })
X0 <- X;
Xproj <- t(X0) %*% projMat;
save(Xproj, file="Xproj.rda")

Xproj[mutWithStrongestEffect,]
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

arcsOrig <-
    as.matrix(
        read_csv("arcsOrig_genes.csv",
                 col_names=unlist(read_csv("geneNamesAfterExprFiltering.list", col_names=F))))

save.image(sprintf("learnMutEffects3_%s.rda", altType))
## load(sprintf("learnMutEffects3_%s.rda", altType))

## 3D "Pareto front" subspace
subspaces <- list()
subspaces[["Front"]] <- projMat;

## healthy -> cancer axis
tmp <-
    as.numeric(
        apply(as.matrix(arcsOrig)[-healthyArch,], 2, mean) -
        as.matrix(arcsOrig)[healthyArch,]
    )
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
## load(sprintf("SSonDistRevCumApprox_%s.rda", altType))
## Results are the same for CNAs and SNVs:
load(sprintf("SSonDistRevCumApprox_%s.rda", "CNAs"))

rm(subspace)
names(subspaces)
## sIdx <- 1
## for (sIdx in 1:length(subspaces)) {

drivers <- list()
overlap <- function(A,B) {
    list(both=intersect(A,B),
         onlyA=setdiff(A,B),
         onlyB=setdiff(B,A))
}
library(gplots)
if ( grepl("brca", tolower(getwd())) ) {
    drivers[["Santarius2010"]] <-
        read.csv("/home/jhausser/work/cancerTaskAtlas/drivers/Santarius2010.csv",
                 as.is=T, h=F)
    colnames(drivers[["Santarius2010"]]) <- c("gene", "type")
    drivers[["Santarius2010"]][,"type"] <- "ONC";
    drivers[["Santarius2010"]] <-
        as_tibble(drivers[["Santarius2010"]]) %>%
        mutate(ONC=T, TSG=F)

    drivers[["Pereira2016"]] <-
        read.table("/home/jhausser/work/cancerTaskAtlas/drivers/Pereira2016.tsv",
                   as.is=T, h=F, sep="\t")
    colnames(drivers[["Pereira2016"]]) <- c("gene", "type")
    drivers[["Pereira2016"]] <-
        as_tibble(drivers[["Pereira2016"]]) %>% mutate(value=T) %>%
        complete(gene, type, fill=list(value=F)) %>%
        spread(type, value)
    ## drivers[["NikZainal2016"]] <-
    ##     read.csv("/home/jhausser/work/cancerTaskAtlas/drivers/NikZainal2016-S14.Type.csv",
    ##              as.is=T, h=T, skip=1)[,-1]
    drivers[["NikZainal2016"]] <-
        read_csv("/home/jhausser/work/cancerTaskAtlas/drivers/NikZainal2016-summarized.csv") %>%
        mutate(effects=str_replace_all(effects, "\\|", "")) %>%
        separate_rows(effects, sep=", ") %>%
        mutate(status=T) %>% complete(gene, effects, fill=list(status=F)) %>%
        spread(effects, status);
    drivers[["NikZainal2016"]] <- 
        drivers[["NikZainal2016"]] %>%
        select(everything(), -matches("^[1-4]")) %>%
        mutate(ONC=amp | missense | inframe,
               TSG=frameshift | HD | nonsense | `ess splice`)
    drivers[["NikZainal2016"]] %>% filter(ONC) %>% count()
    drivers[["NikZainal2016"]] %>% filter(TSG) %>% count()
    drivers[["NikZainal2016"]] %>% count()
    drivers[["NikZainal2016"]] %>% filter(ONC & TSG)
    drivers[["NikZainal2016"]] %>% filter(!ONC & !TSG)

    sapply(drivers, function(x) { x %>% count() })
    pdf("driversVenn.pdf", height=5, width=5)
    venn(sapply(drivers, function(x) { as.data.frame(x)[,1] }))
    dev.off()
} else {
    drivers[["intogen"]] <- 
        read_delim(
            sprintf("intogen-%s-drivers-data.tsv",
                    sub("^.*TCGA.([A-Z]*)_UCSC", "\\1", getwd())),
            delim="\t") %>% filter(KNOWN_DRIVER=="True") %>%
        select(SYMBOL) %>% distinct() %>%
        mutate(ONC=T, TSG=T) %>% rename(gene=SYMBOL)
}

geneChr <-
    read_tsv("~/work/cancerTaskAtlas/hs_gene_coordinates_biomart.tsv", col_names=T) %>%
    dplyr::select(-one_of("ID", "description")) %>%
    filter(grepl("^[0-9]*$", chr) | chr == "X" | chr == "Y" ) %>%
    mutate(chrNum=parse_integer(chr)) %>%
    rename(gene=name)
distinct(geneChr[,"chr"]) %>% print(n=100)
nLoci <-
    select(geneChr, gene) %>% group_by(gene) %>%
    summarize(n=n()) %>% mutate(uniqueLocus=(n<=10))
inner_join(geneChr, nLoci) %>% arrange(desc(n), gene) %>%
    select(gene,n) %>% distinct()  %>% print(n=100)

geneChrU <- inner_join(geneChr, nLoci) %>% filter(uniqueLocus)
## geneChr[,"chr"] <- as.numeric(geneChr[,"chr"])
chrLens <-
    group_by(geneChr, chr) %>%
    summarize(chrLen=max(end), chrNum=mean(chrNum)) %>%
    arrange(chrNum) %>%
    mutate(offsets=lag(cumsum(as.numeric(chrLen)))) %>%
    replace_na(replace=list(offsets=0))
print(chrLens, n=100)

genomeLen <- as.numeric(
    summarize(chrLens, sum(as.numeric(chrLen))))

universeGeneNames <- unique(sub("=.*$", "", rownames(Mp)))
qValThreshold <- .01;

clusterGenomCoords <- function(coords, expand=1) {
    if ( length(coords) == 1 ) { return(0); }
    ordIdcs <- order(coords)
    clustIdx <- rep(-1, length(coords))
    clustIdx[1] <- 0;
    i <- 2
    for (i in 2:length(ordIdcs)) {
        if ( abs(coords[ordIdcs[i]] - coords[ordIdcs[i-1]]) < expand ) {
            clustIdx[i] <- clustIdx[i-1]
        } else {
            clustIdx[i] <- clustIdx[i-1] + 1
        }
    }
    oriClustIdx <- numeric(length(clustIdx))
    oriClustIdx[ordIdcs] <- clustIdx
    return(oriClustIdx)
}


s <- names(subspaces)[1]
## s <- names(subspaces)[2]
mutsAlongFront <- list()
## for (s in names(subspaces)) {
cat(sprintf("%s\n", s))
subspace <- subspaces[[s]]

## Load randomized controls
load(sprintf("SSonDistSamples_%s_%s.rda", altType, s))
rratios <- ratios; rm(ratios);

##################################################

XbackProj <- as.matrix(subspace) %*% t( t(X) %*% as.matrix(subspace) )
res <- sapply(colnames(X), function(i) {
    ##residuals: SS off Pareto front
    sum((X[,i] - XbackProj[,i])^2)
})
onf <- sapply(colnames(X), function(i) {
    ##SS on Parto front
    sum(XbackProj[,i]^2)
})
ratios <- onf / (onf+res)
norms <- sqrt(onf + res)
normsOnFront <- sqrt(onf)

nPos <- apply(Mp[commonMuts,], 1, sum)

x <- Mp[commonMuts,][1,]
## x <- Mp[commonMuts,][1082,]
nPosIdx <- apply(Mp[commonMuts,], 1, function(x) {
    which( sum(x) <= nPoss )[1] - 1
})
if ( sum(is.na(nPosIdx)) ) {
    cat("These alterations occur at frequencies than we didn't sample:\n")
    cat(paste(names(nPosIdx[is.na(nPosIdx)]), collapse="\n"))
    nPosIdx[is.na(nPosIdx)] <- length(nPoss) - 1
}

## plot(apply(Mp[commonMuts,], 1, sum), nPoss[nPosIdx]); abline(0,1)
i <- 1
## i <- 1082
log10pSS <- sapply(1:length(norms), function(i) {
    ## sum(pSSonParms[[s]][,nPosIdx[i]] * c(1, onf[i]))
    log10p <- sum(pSSonParms[["SS"]][,nPosIdx[i]]
                  * c(1, norms[i]^2))
    if ( log10p > 0 ) { return (0) }
    return(log10p)
})
log10pSSonFront <-
    sapply(1:length(normsOnFront), function(i) {
        ## sum(pSSonParms[[s]][,nPosIdx[i]] * c(1, onf[i]))
        log10p <- sum(pSSonParms[["SSonFront"]][,nPosIdx[i]]
                      * c(1, normsOnFront[i]^2))
        if ( log10p > 0 ) { return (0) }
        return(log10p)
    })
i <- which(names(ratios) == "TP53=1")
i <- which(names(ratios) == "GATA3=1")
log10pRatios <- sapply(1:length(ratios), function(i) {
    x <- ratios[i]
    n <- length(pSSonParms[[s]])
    log10p <- sum(pSSonParms[[s]] * x^seq(0,n-1))
    if ( log10p > 0 ) { return (0) }
    return(log10p)
})
summary(log10pSS)
summary(log10pSSonFront)
summary(log10pRatios)
## hist(10^log10pRatios)
## Bonferroni correction
## log10q <- log10p + log10(length(onf))
## log10q[log10q>0] <- 0
## hist(log10q)
sum(10^log10pSS < qValThreshold/length(norms))
sum(10^log10pSSonFront < qValThreshold/length(norms))
sum(10^log10pRatios < qValThreshold/length(ratios))

log10qSS <- log10(fdrtool(10^log10pSS, statistic="pvalue")$lfdr)
## hist(10^log10qSS, seq(0,1,by=.1), col="lightblue")
## hist(10^log10pSS, seq(0,1,by=.1), add=T)
## plot(log10pSS, log10qSS); abline(0,1)

log10qSSonFront <- log10(fdrtool(10^log10pSSonFront,
                                 statistic="pvalue")$lfdr)
## hist(10^log10qSSonFront, seq(0,1,by=.1), col="lightblue")
## hist(10^log10pSSonFront, seq(0,1,by=.1), add=T)
## plot(log10pSSonFront, log10qSSonFront); abline(0,1)

log10qRatios <- log10(fdrtool(10^log10pRatios, statistic="pvalue")$lfdr)
if ( diff(range(log10qRatios)) < 1e-6 ) {
    cat("FDRs don't vary; fdrtool probably failed. Falling back on bonferroni.\n")
    log10qRatios <- log10pRatios + log10(length(log10pRatios))
    log10qRatios[log10qRatios>0] <- 0
}
## hist(10^log10qRatios, seq(0,1,by=.1), col="lightblue")
## hist(10^log10pRatios, seq(0,1,by=.1), add=T)
## plot(log10pRatios, log10qRatios); abline(0,1)

sort(log10qRatios)[1:10]

## SSlist <-
##     mclapply(1:5, function(i) {
##         cat(sprintf("Iteration %d\n", i))

##         ## summary(apply(Mp[commonMuts,], 1, sum))
##         Mr <- t(apply(Mp[commonMuts,], 1, sample))
##         ## Xr <-
##         ##     E %*% t(Mr) %*% diag(1 / apply(Mr, 1, sum)) -
##         ##     E %*% t(1-Mr) %*% diag(1 / apply(1-Mr, 1, sum))
##         MrP <- apply(Mr, 1, function(x) { x / sum(x) })
##         MrN <- apply(1-Mr, 1, function(x) { x / sum(x) })
##         Xr <- E %*% ( MrP - MrN )
##         colnames(Xr) <- rownames(Mr)

##         XrbackProj <- as.matrix(subspace) %*% t( t(Xr) %*% as.matrix(subspace) )

##         resr <- sapply(colnames(X), function(i) {
##             ##residuals: SS off Pareto front
##             sum(((Xr[,i] - mean(Xr[,i])) - XrbackProj[,i])^2)
##         })
##         onfr <- sapply(colnames(X), function(i) {
##             ##SS on Parto front
##             sum(XrbackProj[,i]^2)
##         })
##         return(list(resr=resr, onfr=onfr,
##                     nPos=apply(Mp[commonMuts,], 1, sum)))
##     })
## onfr <- (sapply(SSlist, function(s) { s[["onfr"]] }))
## resr <- (sapply(SSlist, function(s) { s[["resr"]] }))
## nPosr <- (sapply(SSlist, function(s) { s[["nPos"]] }))
## rm(SSlist);
## rratios <- onfr/(onfr+resr);

mybreaks <- seq(min(c(ratios, rratios)), max(c(ratios, rratios)), len=20);
pdf(sprintf("fracSSonVsOff_%s_%s.pdf", altType, s), height=4, width=4)
hist(ratios,
     breaks=mybreaks, main="",
     col="lightblue", freq=F, ## ylim=c(0,5),
     xlim=c(0,1),
     ylim=c(0, max(table(cut(rratios, breaks=mybreaks)) /
                   (length(rratios)) / diff(mybreaks)[1])),
     xlab=sprintf("fraction of mutation effect aligned with %s", s))
hr <- hist(rratios, breaks=mybreaks, border="red",
           col=NULL, add=T, freq=F)
if ( s == "Front" ) {
    abline(v=sum(EPCAeig[1:3])/sum(EPCAeig), lty=2)
}
legend("topright", c(altType, "shuffled"),
       fill=c("lightblue", "white"), col=c("black", "red"),
       bg="white")
dev.off();

allMuts <- tibble(
    alteration=names(ratios),
    fracOnFront=ratios,
    log10qFracOnFront=log10qRatios,
    normOnFront=normsOnFront,
    log10qNormOnFront=log10qSSonFront,
    norm=norms,
    log10qNorm=log10qSS,
    frequency=apply(Mp[commonMuts,], 1, sum)) %>%
    mutate(gene=str_replace(alteration, "=.*$",""))
allMuts <-
    mutate(allMuts,
           isAln=log10qNormOnFront<log10(qValThreshold)
           ## isAln=log10qFracOnFront<log10(qValThreshold)
           ## isAln=log10qFracOnFront<log10(qValThreshold) &
           ##     log10qNorm<log10(qValThreshold)
           )

mutsAlongFront[[s]] <-
    allMuts %>% filter(isAln) %>%
    arrange(log10qFracOnFront, desc(fracOnFront))

allMuts %>% filter(grepl("ERBB2$", gene))
allMuts %>% filter(grepl("PTEN$", gene))
allMuts %>% filter(grepl("TP53$", gene))
allMuts %>% filter(grepl("MYC$", gene))

allMuts %>% filter(grepl("AKT1$", gene))
allMuts %>% filter(grepl("BRAF$", gene))
allMuts %>% filter(grepl("DDR2$", gene))
allMuts %>% filter(grepl("EGFR$", gene))
allMuts %>% filter(grepl("KRAS$", gene))

if ( isMetabric ) {
    rndGenes <- read_tsv("random_genes.list", col_names=F) %>%
        rename(gene=X1) %>% mutate(isRndGene=T)
    allMuts <- left_join(allMuts, rndGenes) %>% replace_na(list(isRndGene=F))
    filter(allMuts, isRndGene) %>% select(gene) %>% distinct() %>% count()
    topGenes <- read_tsv("top_cancer_genes.list", col_names=F) %>%
        rename(gene=X1) %>% mutate(isTopGene=T)
    allMuts <- left_join(allMuts, topGenes) %>% replace_na(list(isTopGene=F))
    filter(allMuts, isTopGene) %>% select(gene) %>% distinct() %>% count()
}

## Add all driving alterations
x <- names(drivers)[1]
allDriversList <-
    lapply(names(drivers), function(x) {
        c(sprintf("%s=1", as.data.frame(filter(drivers[[x]], ONC))[,"gene"]),
          sprintf("%s=2", as.data.frame(filter(drivers[[x]], ONC))[,"gene"]),
          sprintf("%s=-1", as.data.frame(filter(drivers[[x]], TSG))[,"gene"]),
          sprintf("%s=-2", as.data.frame(filter(drivers[[x]], TSG))[,"gene"]))
    })
names(allDriversList) <- names(drivers)
allDrivers <- unique(unlist(allDriversList))

allMuts <- allMuts %>%
    mutate(isDrivingAlteration=alteration %in% allDrivers)
if ( ! isMetabric ) {
    allMuts <- allMuts %>%
        mutate(isTopGene=isDrivingAlteration, isRndGene=!isDrivingAlteration)
}
## filter(allMuts, isDrivingAlteration) %>% print(n=100)

## Switch driver status to false if sibling alteration (same direction,
## different copy number) is parallel to front
## g <- "AKT2"
## g <- "ZNF703"
g <- "BRCA1"
for (g in unlist(allMuts %>%
                 filter(isDrivingAlteration) %>% select(gene) %>%
                 distinct())) {
    ## filter(allMuts, gene==g) %>%
    ##     select(alteration, log10qFracOnFront, isAln, isDrivingAlteration)
    cat(sprintf("Processing %s\n", g))

    ## If at least one alteration of a gene is aligned with the front,
    ## remove driver status from alterations not aligned with the front
    if ( filter(allMuts, gene==g) %>%
         select(alteration, log10qFracOnFront, isAln, isDrivingAlteration) %>%
         filter(isAln) %>% count() > 0 ) {
        allMuts[which(allMuts[,"gene"] == g & !allMuts[,"isAln"]),
                "isDrivingAlteration"] <- F
    }
}

minFreq <- .05 * ncol(Mp)

pdf("scatterDriverPredictors.pdf", height=8, width=8)
plot(
    as.data.frame(
        filter(allMuts, frequency>minFreq))[,c("fracOnFront", "log10qNorm",
                                               "log10qNormOnFront",
                                               "frequency")])
dev.off();

lm1 <- glm(isDrivingAlteration ~ log10qNormOnFront * log10qFracOnFront * frequency,
           family=binomial, data=allMuts)
summary(lm1)
allMutsGlm <- add_column(allMuts, glm=predict(lm1, newdata=allMuts))

p1 <- ggplot(data=allMutsGlm) +
    geom_violin(aes(x=isDrivingAlteration, y=log10qNormOnFront)) +
    theme_gray()
p2 <- ggplot(data=allMutsGlm) +
    geom_violin(aes(x=isDrivingAlteration, y=log10qNorm)) +
    theme_gray()
p3 <- ggplot(data=allMutsGlm) +
    geom_violin(aes(x=isDrivingAlteration, y=log10qFracOnFront)) +
    theme_gray()
p4 <- ggplot(data=allMutsGlm) +
    geom_violin(aes(x=isDrivingAlteration, y=glm)) +
    theme_gray()
plot_grid(p1, p2, p3, p4, nrow=2, align="v")
ggsave("predictDriver_fourPredictors.pdf", height=8, width=8)

filter(allMuts, isDrivingAlteration) %>% arrange(desc(fracOnFront))
ggplot(allMuts %>% select(fracOnFront, isDrivingAlteration) %>%
       filter(isDrivingAlteration) %>%
       bind_rows(tibble(fracOnFront=rratios,
                        isDrivingAlteration=F))) +
    geom_violin(aes(x=isDrivingAlteration, y=100*fracOnFront)) +
    geom_hline(aes(yintercept=100*sum(EPCAeig[1:3])/sum(EPCAeig)),
               linetype=2) +
    labs(x="Is driving alteration?", y="% alignment with Pareto front") +
    theme_gray()
ggsave(sprintf("%s_alnWithFront.pdf", altType), height=4, width=4);
wilcox.test(rratios,
            unlist(filter(allMuts, isDrivingAlteration) %>% select(fracOnFront)))

## toCmp <- rep(NA, nrow(allMutsGlm))
## toCmp[allMutsGlm$isRndGene & !allMutsGlm$isDrivingAlteration] <- "random genes";
## toCmp[allMutsGlm$isDrivingAlteration] <- "driver genes";
## ggplot(allMutsGlm %>% mutate(toCmp=factor(toCmp)) %>% drop_na(toCmp)) +
##     geom_violin(aes(x=toCmp, y=100*fracOnFront)) +
##     geom_hline(aes(yintercept=100*sum(EPCAeig[1:3])/sum(EPCAeig)),
##                linetype=2) +
##     labs(x="Is driving alteration?", y="% alignment with Pareto front") +
##     theme_gray()
## ggsave(sprintf("%s_alnWithFront.pdf", altType), height=4, width=4);
## wilcox.test(
##     unlist(filter(allMutsGlm, isDrivingAlteration)[,"fracOnFront"]),
##     unlist(filter(allMutsGlm, isRndGene & !isDrivingAlteration)[,"fracOnFront"]))

qCutOffs <- seq(quantile(as.data.frame(allMutsGlm)[,"glm"], .01),
                quantile(as.data.frame(allMutsGlm)[,"glm"], .99),
                len=100)
qCutOff <- median(qCutOffs)
ROC <-
    t(sapply(qCutOffs,
             function(qCutOff) {
                 c(specificity=as.numeric(allMutsGlm %>%
                                          filter(frequency>=minFrequency) %>% 
                                          filter(glm>=qCutOff) %>%
                                          select(isDrivingAlteration, log10qNormOnFront) %>% 
                                          summarize(specificity=mean(isDrivingAlteration))),
                   sensitivity=as.numeric(allMutsGlm %>%
                                          filter(frequency>=minFrequency) %>% 
                                          filter(isDrivingAlteration) %>%
                                          transmute(caught=glm>=qCutOff) %>% 
                                          summarize(sensitivity=mean(caught))))
             }))

pdf("ROCanalysis.pdf", height=3, width=3*3)
par(mfrow=c(1,3))
plot(qCutOffs, ROC[,"specificity"], type="l",
     xlab="cut-off", ylab="specificity");
## abline(v=log10(qValThreshold), lty=2, col="grey")

plot(qCutOffs, ROC[,"sensitivity"], type="l",
     xlab="cut-off", ylab="sensitivity");
## abline(v=log10(qValThreshold), lty=2, col="grey")

plot(ROC[,"sensitivity"], ROC[,"specificity"],
     xlab="sensitivity", ylab="specificity", type="l",
     xlim=c(0, max(ROC[,"sensitivity"])), ylim=c(0, max(ROC[,"specificity"])))
dev.off()

## test for match between drivers and found mutations    
driversIdx <- names(drivers)[1]
pdf(sprintf("vennAlnDrivers_%s_%s.pdf", altType, s),
    height=4*length(allDriversList), width=8)
par(mfrow=c(length(allDriversList),2))

driversOverlap <-
    sapply(names(drivers), function(driversIdx) {
        ## Only look at drivers that are in the universe
        tDrivers <- allDriversList[[driversIdx]]

        ## Overlap in top cancer genes set
        allMuts <- mutate(allMuts, isDrivingAlteration=alteration %in% tDrivers)
        
        ## total white balls
        m <-
            as.numeric(
                filter(allMuts, isTopGene & isAln) %>% 
                select(gene) %>% distinct() %>% count()
            )
        
        ##total black balls
        n <-
            as.numeric(
                filter(allMuts, isTopGene & !isAln) %>% 
                select(gene) %>% distinct() %>% count()
            )
        
        ##number of drawn balls
        k <-
            as.numeric(
                filter(allMuts, isTopGene & isDrivingAlteration) %>% 
                select(gene) %>% distinct() %>% count()
            )

        ##number of while balls in the draw
        q <-
            as.numeric(
                filter(allMuts, isTopGene & isAln & isDrivingAlteration) %>% 
                select(gene) %>% distinct() %>% count()
            )
        
        retVal <-
            c(sapply(
                overlap(
                    A=unlist(filter(allMuts, isTopGene & isAln) %>% 
                             select(gene) %>% distinct()),
                    B=unlist(filter(allMuts, isTopGene & isDrivingAlteration) %>% 
                             select(gene) %>% distinct())), length),
              p=phyper(q=q-1, m=m, n=n, k=k, lower.tail=F))
        
        vennList <-
            list(unlist(filter(allMuts, isTopGene & isAln) %>%
                        select(gene) %>% distinct()),
                 unlist(filter(allMuts, isTopGene & isDrivingAlteration) %>%
                        select(gene) %>% distinct()),
                 unlist(filter(allMuts, isTopGene) %>%
                        select(gene) %>% distinct()))
        names(vennList) <- c("aligned genes", driversIdx, "all genes");
        venn(vennList)
        title("top cancer genes",
              sub=sprintf("p = %.3f",
                          phyper(q=q-1, m=m, n=n, k=k, lower.tail=F)))
        
        ## Overlap in random genes set

        ## total white balls
        m <-
            as.numeric(
                filter(allMuts, isRndGene & isAln) %>% 
                select(gene) %>% distinct() %>% count()
            )
        
        ##total black balls
        n <-
            as.numeric(
                filter(allMuts, isRndGene & !isAln) %>% 
                select(gene) %>% distinct() %>% count()
            )
        
        ##number of drawn balls
        k <-
            as.numeric(
                filter(allMuts, isRndGene & isDrivingAlteration) %>% 
                select(gene) %>% distinct() %>% count()
            )

        ##number of while balls in the draw
        q <-
            as.numeric(
                filter(allMuts, isRndGene & isAln & isDrivingAlteration) %>% 
                select(gene) %>% distinct() %>% count()
            )
        
        c(sapply(
            overlap(
                A=unlist(filter(allMuts, isRndGene & isAln) %>%
                         select(gene) %>% distinct()),
                B=unlist(filter(allMuts, isRndGene & isDrivingAlteration) %>%
                         select(gene) %>% distinct())),
            length),
          p=phyper(q=q, m=m, n=n, k=k, lower.tail=F))

        vennList <-
            list(unlist(filter(allMuts, isRndGene & isAln) %>%
                        select(gene) %>% distinct()),
                 unlist(filter(allMuts, isRndGene & isDrivingAlteration) %>%
                        select(gene) %>% distinct()),
                 unlist(filter(allMuts, isRndGene) %>%
                        select(gene) %>% distinct()))
        names(vennList) <- c("aligned genes", driversIdx, "all genes");
        venn(vennList)
        title("1000 random genes",
              sub=sprintf("p = %.3f",
                          phyper(q=q, m=m, n=n, k=k, lower.tail=F)))
        return(retVal)
    })
dev.off()

## Add 'driver' indicator to table
toAdd <- names(which(driversOverlap["p",] < .01))
if ( length(toAdd) > 0 ) {
    x <- toAdd[1]
    for ( x in toAdd ) {
        drvVarName <- sprintf("%s (driver)", x)
        allMuts[[drvVarName]] <-
            with(allMuts, alteration %in% allDriversList[[x]])
    }
}
## myDrivers <- unique(unlist(sapply(toAdd, function(x) { drivers[[x]] })))

## allMuts %>% group_by(isTopGene) %>%
##     summarize(nAln=sum(isAln), total=n())
## allMuts %>% group_by(isRndGene) %>%
##     summarize(nAln=sum(isAln), total=n())

allMutsChr <- inner_join(allMuts, geneChrU) %>% filter(isTopGene)
## %>% select(-isTopGene, -isRndGene)

allMutsChrD <-
    left_join(allMutsChr, chrLens, by="chr") %>%
    select(-chrNum.x) %>% rename(chrNum=chrNum.y) %>% 
    ## filter(!is.na(chr)) %>%
    mutate(pos=start+offsets)

expandClusterBy <- 5e6;
mychr <- "1"
mychr <- "16"
for (mychr in unlist(distinct(allMutsChrD, chr))) {
    gIdcs <- which(allMutsChrD[,"chr"] == mychr &
                   ( allMutsChrD[,"isAln"] | allMutsChrD[,"isDrivingAlteration"] ))
    if ( length(gIdcs) == 0) { next; }
    myClusts <-
        clusterGenomCoords(
        (unlist(allMutsChrD[gIdcs,"start"]) +
         unlist(allMutsChrD[gIdcs,"end"]))/2,
        expand=expandClusterBy)
    allMutsChrD[gIdcs, "chrCluster"] <- paste(myClusts, mychr, sep="-")
}

## Let's find out if random genes are in clusters of possible
## cancer genes
clusters <-
    group_by(allMutsChrD, chrCluster) %>% drop_na(chrCluster) %>% 
    summarize(start=min(start), end=max(end), chr=first(chr))
clusters %>% arrange(chr,start) %>% print(n=120)
ggplot(clusters %>% mutate(len=end-start)) +
    geom_histogram(mapping=aes(x=len)) +
    scale_x_log10()
clusters %>% arrange(chr, start) %>%
    mutate(len=end-start) %>% summarize(totLen=sum(len)) / 3e9
## with 5e6 expansions, cancer clusters make up 17% of the genome

alignedRndGenes <-
    inner_join(allMuts, geneChrU) %>% filter(isRndGene) %>%
    select(gene, start, end, chr) %>% distinct()
gv <- alignedRndGenes[2,]
isInKnownCluster <-
    apply(alignedRndGenes, 1,
          function(gv) {
              filter(clusters, start - expandClusterBy < gv["start"] &
                               gv["end"] < end + expandClusterBy &
                               chr == gv["chr"]) %>% count() > 0
          })
mean(isInKnownCluster)
## 63% of random genes are located in the 17% of the genome in which
## cancer genes (drivers + top 1000 eQTL) are located.

## allMutsChrD %>% filter(frequency>700) %>% summarize(chr=unique(chrNum))
## All the high-frequency alterations are on chr 1

mychr <- as.data.frame(chrLens)[1,"chr"]
mychr <- "1"
mychr <- "17"
mychr <- "11"
myalpha <- .65
## for (mychr in as.data.frame(chrLens)[,"chr"]) {
##     if ( filter(allMutsChrD, chr==mychr) %>% count() == 0 ) {
##         next;
##     }
##     xlab <- sprintf("Position on chromosome %s [MB]", mychr)
##     tMaxFreq <-
##         as.numeric(
##             filter(allMutsChrD, chr==mychr) %>%
##             ## select(frequency) %>% 
##             summarize(maxFreq=max(frequency)))
##     p1 <-
##         ggplot(filter(allMutsChrD, chr==mychr)) +
##         geom_point(mapping=aes(x=start/1e6, y=frequency, color=isAln)) +
##         coord_cartesian(xlim=c(0,
##                                as.numeric(filter(chrLens, chr==mychr) %>%
##                                           select(chrLen))/1e6)) +
##         geom_line(data=clusters %>% filter(chr==mychr) %>%
##                       gather(`start`, `end`, key=type, value=pos),
##                   mapping=aes(x=pos/1e6, y=tMaxFreq*1.1,
##                               group=chrCluster)) +
##         geom_text(data=clusters %>% filter(chr==mychr) %>%
##                       mutate(pos=(start+end)/2),
##                   mapping=aes(x=pos/1e6, y=tMaxFreq*1.2,
##                               label=chrCluster)) +
##         geom_label_repel(data=filter(allMutsChrD,
##                                      chr==mychr & isDrivingAlteration),
##                          mapping=aes(x=start/1e6, y=frequency,
##                                      label=alteration, color=isAln),
##                          alpha=myalpha, show.legend=F) +
##         labs(x=xlab) +
##         theme_gray()

##     p2 <-
##         ggplot(filter(allMutsChrD, chr==mychr)) +
##         geom_point(mapping=aes(x=start/1e6, y=-log10qFracOnFront, color=isAln)) +
##         coord_cartesian(xlim=c(0,
##                                as.numeric(filter(chrLens, chr==mychr) %>%
##                                           select(chrLen))/1e6)) +
##         geom_label_repel(data=filter(allMutsChrD, chr==mychr & isDrivingAlteration),
##                          mapping=aes(x=start/1e6, y=-log10qFracOnFront,
##                                      label=alteration, color=isAln),
##                          alpha=myalpha, show.legend=F) +
##         labs(x=xlab) +
##         theme_gray()

##     p3 <-
##         ggplot(filter(allMutsChrD, chr==mychr)) +
##         geom_point(mapping=aes(x=start/1e6, y=-log10qNorm, color=isAln)) +
##         coord_cartesian(xlim=c(0,
##                                as.numeric(filter(chrLens, chr==mychr) %>%
##                                           select(chrLen))/1e6)) +
##         geom_label_repel(data=filter(allMutsChrD, chr==mychr & isDrivingAlteration),
##                          mapping=aes(x=start/1e6, y=-log10qNorm,
##                                      label=alteration, color=isAln),
##                          alpha=myalpha, show.legend=F) +
##         labs(x=xlab) +
##         theme_gray()

##     p4 <-
##         ggplot(filter(allMutsChrD, chr==mychr)) +
##         geom_point(mapping=aes(x=start/1e6,
##                                y=-log10qNormOnFront,
##                                color=isAln)) +
##         coord_cartesian(xlim=c(0,
##                                as.numeric(filter(chrLens, chr==mychr) %>%
##                                           select(chrLen))/1e6)) +
##         geom_label_repel(data=filter(allMutsChrD, chr==mychr & isDrivingAlteration),
##                          mapping=aes(x=start/1e6,
##                                      y=-log10qNormOnFront,
##                                      label=alteration, color=isAln),
##                          alpha=myalpha, show.legend=F) +
##         labs(x=xlab) +
##         theme_gray()

##     p5 <-
##         ggplot(filter(allMutsChrD, chr==mychr)) +
##         geom_point(mapping=aes(x=start/1e6, y=100*fracOnFront## , color=isAln
##                                )) +
##         coord_cartesian(xlim=c(0,
##                                as.numeric(filter(chrLens, chr==mychr) %>%
##                                           select(chrLen))/1e6)) +
##         geom_label_repel(data=filter(allMutsChrD, chr==mychr & isDrivingAlteration),
##                          mapping=aes(x=start/1e6, y=100*fracOnFront,
##                                      label=alteration, color=isAln),
##                          alpha=myalpha, show.legend=F) +
##         labs(x=xlab, y="% alignment with Pareto front") +
##         theme_gray()
    
##     plot_grid(p1, p2, p3, p4, p5, labels=toupper(letters[1:5]),
##               nrow=5, ncol=1, align="v")
##     ggsave(sprintf("chrPos_frequency_fracOnFront_log10qNorm_%s_%s_%s.pdf",
##                    altType, s, mychr),
##            height=3*5, width=12*1)
## }

allMutsChrD[allMutsChrD[,"chrCluster"] == "0-14",
            c("chr", "start")] %>% arrange(start)

## Chr clusters which have at least one gene parallel to the front
allMutsChrD %>% filter(isAln) %>% select(chrCluster) %>% distinct() %>% count()
## All regions
select(allMutsChrD, chrCluster) %>% distinct() %>% count()

write_csv(allMutsChrD %>% filter(isAln | isDrivingAlteration) %>%
          arrange(chrNum, chrCluster, ## gene,
                  log10qNormOnFront) %>%
          mutate(`start (Mb)`=round(start/1e6, 2)) %>%
          mutate(`end (Mb)`=round(end/1e6, 2)) %>%
          select(c(-pos,-offsets,-chrNum, -chrLen, -gene,
                   -uniqueLocus, -n, -isTopGene, -isRndGene,
                   -start, -end)),
          sprintf("onFront_%s_%s.csv", altType, s))

## What are genes whose amp and del are parallel to the front and
## are known drivers?
ampDelAlnGenes <-
    allMutsChrD %>%
    filter((isAln | isDrivingAlteration) &
           ## log10qNorm<log10(qValThreshold) &
           frequency>minFreq) %>%
    group_by(gene) %>%
    summarize(hasAmp=sum(grepl("=[12]", alteration) & isAln)>0,
              hasDel=sum(grepl("=-[12]", alteration) & isAln)>0,
              isDrivingAlteration=sum(isDrivingAlteration)>0) %>%
    filter(hasAmp & hasDel & isDrivingAlteration) %>%
    select(gene) %>% mutate(mm=T)
ampDelAlnTab <-
    inner_join(allMutsChrD, ampDelAlnGenes) %>%
    ## filter((isAln|isDrivingAlteration) &
    ##        log10qNorm<log10(qValThreshold) &
    ##        frequency>minFreq) %>% 
    arrange(gene) %>%
    select(gene, alteration, fracOnFront, log10qFracOnFront,
           normOnFront, log10qNormOnFront,
           norm, log10qNorm,
           frequency, isAln, isDrivingAlteration, chrCluster)
ampDelAlnTab
write_csv(ampDelAlnTab, sprintf("ampDelAlnTab_%s_%s.csv", altType, s))

ampDelAlnTabPlot <- 
    ampDelAlnTab %>%
    separate(alteration, into=c("tmp", "dosage"), sep="=") %>%
    mutate(dosage=factor(as.numeric(dosage), order=T)) %>%
    arrange(gene, dosage)
## Fix cluster names
gene2clust <-
    ampDelAlnTabPlot %>% select(gene, chrCluster) %>% drop_na() %>% distinct()
ampDelAlnTabPlot <-
    inner_join(select(ampDelAlnTabPlot, -chrCluster, -tmp),
               gene2clust) %>%
    mutate(gene=sprintf("%s (%s)", gene, chrCluster))

plot_grid(
    ggplot(data=ampDelAlnTabPlot) +
    geom_point(aes(x=dosage, y=-log10qNormOnFront)) +
    facet_wrap(~gene) + theme_gray(),
    
    ggplot(data=ampDelAlnTabPlot) +
    geom_point(aes(x=dosage, y=fracOnFront)) +
    facet_wrap(~gene) + theme_gray(),

    ggplot(data=ampDelAlnTabPlot) +
    geom_point(aes(x=dosage, y=-log10qNorm)) +
    facet_wrap(~gene) + theme_gray(),
    
    ggplot(data=ampDelAlnTabPlot) +
    geom_point(aes(x=dosage, y=frequency)) +
    facet_wrap(~gene) + theme_gray(),
    
    align="v", nrow=2)
ggsave("geneDoseResponses.pdf", height=6, width=8);

topDriverInEachCluster <-
    filter(allMutsChrD, isAln & isDrivingAlteration) %>%
    select(alteration, chrCluster, log10qNormOnFront) %>% 
    group_by(chrCluster) %>% top_n(1, -log10qNormOnFront) %>% 
    summarize(alteration=first(alteration)) %>%
    mutate(isTopDriverInCluster=T)
allMutsChrD <-
    left_join(allMutsChrD, topDriverInEachCluster) %>%
    replace_na(list(isTopDriverInCluster=F)) %>%
    mutate(chrClusterFact=factor(chrCluster)) ## for colors

allMutsChrD  %>% group_by(chrCluster)  %>%
    summarize(isThereDriver=max(isTopDriverInCluster)) %>%
    summarize(nWithDriver=sum(isThereDriver),
              nWithoutDriver=sum(!isThereDriver))

ggplot(filter(allMutsChrD, isDrivingAlteration | isAln)) +
    geom_point(mapping=aes(x=pos/1e9, y=-log10qNormOnFront, shape=isDrivingAlteration,
                           col=chrClusterFact), show.legend=F) +
    geom_vline(data=chrLens,
               mapping=aes(xintercept=offsets/1e9), linetype=2) +
    geom_hline(mapping=aes(yintercept=-log10(qValThreshold))) +
    geom_label_repel(mapping=aes(x=pos/1e9, y=-log10qNormOnFront,
                                 label=alteration, col=chrClusterFact),
                     data=filter(allMutsChrD, isTopDriverInCluster),
                     show.legend=F, alpha=.8) +
    geom_label_repel(mapping=aes(x=offsets/1e9, y=.3, label=chr),
                     data=chrLens,
                     nudge_x=2.5e7/1e9, alpha=.5) +
    ## coord_cartesian(xlim=c(0,2.8e9)) +
    labs(x="Genomic position [GB]", y=quote(paste(log[10], " ", q))) +
    scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
    theme_gray()
ggsave(sprintf("chrAlnFrontDrivers_%s_%s.pdf", altType, s), height=3, width=12)

## distribution of distances to nearest driver
if ( altType == "CNAs" ) {
    g <- "NUP93"
    g <- "FBXO31"

    getDistToNextDriver <- function(g, allMutsChr) {
        ## cat(sprintf("%s\n", g))
        ## browser();
        gPos <-
            as.numeric(filter(allMutsChr, gene == g) %>%
                       select(start) %>% summarize(start=mean(start)))

        tChr <-
            filter(allMutsChr, gene == g) %>%
            select(chr) %>% distinct()
        if ( count(tChr) != 1 ) { return(NA) }
        gChr <- as.character(tChr)

        driversOnChr <-
            allMutsChr %>% filter(chr == gChr & isDrivingAlteration) %>%
            select(start) %>% distinct()
        
        if ( count(driversOnChr) == 0 ) {
            return(as.numeric(filter(chrLens, chr == gChr) %>% select(chrLen)))
        } else {
            return(min(abs(gPos - as.numeric(unlist(driversOnChr)))))
        }
    }

    rndPosDists <-
        sapply(1:1000, function(i) {
            gChr <- sample(unlist(chrLens[,"chr"]), size=1)
            gPos <- runif(1, min=0,
                          max=as.numeric(chrLens[chrLens[,"chr"] == gChr, "chrLen"]))
            driversOnChr <- allMutsChrD %>%
                filter(chr == gChr & isDrivingAlteration)
            if ( nrow(driversOnChr) == 0 ) {
                return(as.numeric(chrLens[chrLens[,"chr"] == gChr,"chrLen"]))
            } else {
                return(min(abs(gPos -
                               unique(driversOnChr[,"start"]))))
            }
        })
    rndPosDists <- log10(rndPosDists)

    topGeneAlnDriver <- 
        filter(allMutsChr, isTopGene) %>% group_by(gene) %>%
        select(gene, isAln, isDrivingAlteration) %>%
        summarize(isAln=sum(isAln)>0, isDrivingAlteration=sum(isDrivingAlteration)>0)
    alnGenes <- unique(unlist(
        filter(topGeneAlnDriver, isAln | isDrivingAlteration) %>% select(gene)
    ))
    alnGenesDists <-
        log10(
            sapply(alnGenes, function(g) {
                getDistToNextDriver(g, allMutsChr) })
        )

    topAlnGenes <- unique(unlist(filter(
        allMutsChr,
        isTopGene & log10qFracOnFront<log10(qValThreshold/1e3))[,"gene"]))
    topAlnGenesDists <-
        log10(
            sapply(topAlnGenes, function(g) {
                getDistToNextDriver(g, allMutsChr) })
        )

    notAlnGenes <- unique(unlist(
        filter(topGeneAlnDriver, !isAln) %>% select(gene)
    ))
    notAlnGenesDists <-
        log10(
            sapply(notAlnGenes, function(g) {
                getDeistToNextDriver(g, allMutsChr) })
        )
    ## allGenes <- unique(unlist(allMutsChrD[,"gene"]))
    ## allGenesDists <- sapply(allGenes, function(g) {
    ##     getDistToNextDriver(g, allMutsChr) })

    rndMutsChr <-
        inner_join(filter(allMuts, isRndGene | isDrivingAlteration), geneChrU)
    rndGenes <-
        unique(unlist(
            filter(rndMutsChr, isRndGene) %>% select(gene, isAln, isDrivingAlteration) %>%
            group_by(gene) %>%
            summarize(isAln=sum(isAln)>0,
                      isDrivingAlteration=sum(isDrivingAlteration)>0) %>%
            filter(!isAln) %>% select(gene)
        ))
    rndGenesDists <-
        log10(
            sapply(rndGenes, function(g) {
                getDistToNextDriver(g, rndMutsChr) })
        )

    getDistToNextDriver("SNORA4", rndMutsChr)

    log10(c(median(topAlnGenesDists, na.rm=T),
            median(alnGenesDists, na.rm=T), median(allGenesDists, na.rm=T),
            median(rndGenesDists, na.rm=T), median(rndPosDists)))

    pdf(sprintf("distToNearestDriver_%s_%s.pdf", altType, s),
        height=5.5, width=5.5)
    plot(ecdf(rndPosDists), col="grey", main="",
         xlab=expression(paste(log[10], " distance to closest driver [nt]")),
         ylab="empirical cumulative distribution",
         xlim=c(3.5, 8.5))

    ## lines(ecdf(topAlnGenesDists), col="red", lwd=2)
    lines(ecdf(alnGenesDists), col="orange", lwd=2)
    ## lines(ecdf(notAlnGenesDists), col="black", lwd=2)
    ## lines(ecdf(alnGenesDists), col="blue")
    lines(ecdf(rndGenesDists), col="blue")

    text(4.5, .4,
         sprintf("aligned vs random\np = %.2e (KS-test)",
                 ks.test(alnGenesDists, rndGenesDists)$p.value))
    
    
    ## text(6e7, .75, "Compared to random genes", col="blue")
    ## text(6e7, .7,
    ##      sprintf("p = %.3f",
    ##              wilcox.test(log10(topAlnGenesDists),
    ##                          log10(rndGenesDists))$p.value),
    ##      col="red")
    ## text(6e7, .65,
    ##      sprintf("p = %.3f",
    ##              wilcox.test(log10(alnGenesDists),
    ##                          log10(rndGenesDists))$p.value),
    ##      col="orange")

    ## text(6e7, .45, "Compared to random locations", col="grey")
    ## text(6e7, .4,
    ##      sprintf("p = %.3f",
    ##              wilcox.test(log10(topAlnGenesDists),
    ##                          log10(rndPosDists))$p.value),
    ##      col="red")
    ## text(6e7, .35,
    ##      sprintf("p = %.3f",
    ##              wilcox.test(log10(alnGenesDists),
    ##                          log10(rndPosDists))$p.value),
    ##      col="orange")
    
    legend("topleft",
           sprintf("%d %s",
                   c(## length(topAlnGenesDists),
                       length(alnGenesDists),
                       length(rndGenesDists), length(rndPosDists)),
                   c(## "top genes aligned to front",
                       "genes aligned to front",
                       "not aligned random genes", "random loci")),
           col=c("orange", "blue", "grey"), lwd=2, bg="white")
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
    
    ## write.csv(mutsAlongFront[[s]][order(mutsAlongFront[[s]][,"q"]),],
    ##           file=sprintf("onFront_%s_%s.csv", altType, s))

    ## pdf(sprintf("fracSSonFront_SS_freq_%s_%s.pdf", altType, s), height=4, width=4)
}
allMutsPlot <- 
    allMutsChrD %>% filter(frequency>=minFreq) %>%
    unite(alnDrv, isAln, isDrivingAlteration, sep=",", remove=F) %>%
    mutate(alnDrv=str_replace(alnDrv,
                              "FALSE,FALSE", "not driver, not aligned")) %>%
    mutate(alnDrv=str_replace(alnDrv,
                              "FALSE,TRUE", "driver but not aligned")) %>%
    mutate(alnDrv=str_replace(alnDrv,
                              "TRUE,FALSE", "not driver but yet aligned")) %>%
    mutate(alnDrv=str_replace(alnDrv,
                              "TRUE,TRUE", "driver and aligned"))

pValContours <- 10^c(-6, -5, -4, -3, -2);
tpc <-
    tibble(contour=pValContours,
           fracOnFront=
               sapply(pValContours, function(myp) {
                   res <- optimize(function(r) {
                       abs(sum(pSSonParms[[s]] *
                               r^seq(0, length(pSSonParms[[s]])-1)) -
                           log10(myp))
                   }, c(0,1))$minimum
               }))

## By cluster
ggplot(allMutsPlot %>% mutate(isAmplified=grepl("=[12]$", alteration))) +
    geom_vline(data=tpc, mapping=aes(xintercept=fracOnFront),
               linetype=2) +
    geom_point(aes(x=fracOnFront, y=frequency, col=chrClusterFact,
                   shape=isDrivingAlteration),
               show.legend=F) +
    geom_text(data=tpc, mapping=aes(x=fracOnFront,
                                    y=summarize(allMutsChrD,
                                                max=.8*max(frequency)),
                                    label=pValContours), angle=90,
              nudge_x=-.02) +
    theme_gray() +
    labs(x=sprintf("fraction of ||v||^2 on %s", s), y="frequency")
ggsave(sprintf("fracSSonFront_freq_byCluster_%s_%s.pdf", altType, s),
       height=4, width=4)

## By driver and amp/del
ggplot(allMutsPlot %>% mutate(isAmplified=grepl("=[12]$", alteration))) +
    geom_vline(data=tpc, mapping=aes(xintercept=fracOnFront),
               linetype=2) +
    geom_point(aes(x=fracOnFront, y=frequency, shape=isDrivingAlteration,
                   col=isAmplified),
               show.legend=T) +
    geom_text(data=tpc, mapping=aes(x=fracOnFront,
                                    y=summarize(allMutsChrD,
                                                max=.8*max(frequency)),
                                    label=pValContours), angle=90,
              nudge_x=-.02) +
    theme_gray() +
    labs(x=sprintf("fraction of ||v||^2 on %s", s), y="frequency")
ggsave(sprintf("fracSSonFront_freq_byDriverAmp_%s_%s.pdf", altType, s),
       height=4, width=7)

p1 <-
    ggplot(allMutsPlot) +
    geom_vline(data=tpc, mapping=aes(xintercept=fracOnFront), linetype=2) +
    geom_point(aes(x=fracOnFront, y=norm^2, col=alnDrv), show.legend=T) +
    geom_text(data=tpc, mapping=aes(x=fracOnFront, y=50,
                                    label=pValContours), angle=90,
              nudge_x=-.02) +
    scale_y_log10() +
    theme_gray() +
    labs(x=sprintf("fraction of ||v||^2 on %s", s),
         y="effect size ||v||^2")

p2 <-
    ggplot(allMutsPlot) +
    geom_vline(data=tpc, mapping=aes(xintercept=fracOnFront), linetype=2) +
    geom_point(aes(x=fracOnFront, y=frequency, col=alnDrv), show.legend=T) +
    geom_text(data=tpc, mapping=aes(x=fracOnFront,
                                    y=summarize(allMutsChrD, max=.8*max(frequency)),
                                    label=pValContours), angle=90,
              nudge_x=-.02) +
    theme_gray() +
    labs(x=sprintf("fraction of ||v||^2 on %s", s), y="frequency")

tcl <- tibble(frequency=0, es=0, p=0)
for (myp in pValContours) {
    xs <-
        sapply(1:length(nPoss), function(i) {
            (log10(myp) - pSSonParms[["SS"]]["a",i]) /
                pSSonParms[["SS"]]["b",i]
        })
    tcl <-
        bind_rows(tcl, tibble(frequency=nPoss, es=xs, p=myp))
}
tcl <- filter(tcl, frequency>0) %>%
    mutate(isLabel=(frequency==max(frequency)) &
               (p==min(p) | p==max(p)))

p3 <-
    ggplot(allMutsPlot) +
    geom_line(data=tcl, mapping=aes(x=es, y=frequency, group=p), linetype=2) +
    geom_point(aes(x=norm^2, y=frequency, col=alnDrv)) +
    geom_text(data=filter(tcl, isLabel),
              mapping=aes(x=es, y=frequency, label=p),
              angle=90, nudge_x=.05, nudge_y=-25) +
    scale_x_log10() +
    theme_gray() +
    labs(x=sprintf("effect size ||v||^2", s), y="frequency")

tcl <- tibble(frequency=0, es=0, p=0)
for (myp in pValContours) {
    xs <-
        sapply(1:length(nPoss), function(i) {
            (log10(myp) - pSSonParms[["SSonFront"]]["a",i]) /
                pSSonParms[["SS"]]["b",i]
        })
    tcl <-
        bind_rows(tcl, tibble(frequency=nPoss, es=xs, p=myp))
}
tcl <- filter(tcl, frequency>0) %>%
    mutate(isLabel=(frequency==max(frequency)) &
               (p==min(p) | p==max(p)))

p4 <- 
    ggplot(allMutsPlot) +
    geom_line(data=tcl, mapping=aes(x=es, y=frequency, group=p), linetype=2) +
    geom_point(aes(x=normOnFront^2, y=frequency, col=alnDrv)) +
    geom_text(data=filter(tcl, isLabel),
              mapping=aes(x=es, y=frequency, label=p),
              angle=90, nudge_x=.05, nudge_y=-25) +
    scale_x_log10() +
    theme_gray() +
    labs(x=sprintf("effect size on front ||v_on||^2", s), y="frequency")

plot_grid(p1, p2, p3, p4, labels=toupper(letters[1:4]),
          nrow=2, ncol=2, align="v")
ggsave(sprintf("fracSSonFront_SS_freq_%s_%s.pdf", altType, s),
       height=4*2, width=7*2)

## Missed drivers
write_csv(
    filter(allMutsPlot,
           isDrivingAlteration &
           (log10qFracOnFront > log10(qValThreshold) &
            log10qNormOnFront > log10(qValThreshold) ) &
           frequency>minFreq),
    sprintf("missedDrivers_%s_%s.csv", altType, s))

## Where are the known driver genes which aren't aligned?
p1 <-
    ggplot(allMutsPlot) +
    geom_vline(data=tpc, mapping=aes(xintercept=fracOnFront), linetype=2) +
    geom_point(aes(x=fracOnFront, y=norm^2, col=alnDrv), show.legend=T) +
    ## geom_text(data=tpc, mapping=aes(x=fracOnFront, y=50,
    ##                                 label=pValContours), angle=90,
    ##           nudge_x=-.02) +
    geom_label_repel(data=filter(allMutsPlot,
                                 !isAln & isDrivingAlteration &
                                 log10qNorm < log10(qValThreshold) &
                                 frequency>minFreq),
                     mapping=aes(x=fracOnFront, y=norm^2,
                                 col=alnDrv, label=alteration),
                     show.legend=F, alpha=.7) +
    scale_y_log10() +
    theme_gray() +
    labs(x=sprintf("fraction of ||v||^2 on %s", s),
         y="effect size ||v||^2")

p2 <-
    ggplot(allMutsPlot) +
    geom_vline(data=tpc, mapping=aes(xintercept=fracOnFront), linetype=2) +
    geom_point(aes(x=fracOnFront, y=frequency, col=alnDrv), show.legend=T) +
    geom_label_repel(data=filter(allMutsPlot,
                                 !isAln & isDrivingAlteration &
                                 log10qNorm < log10(qValThreshold) &
                                 frequency>minFreq),
                     mapping=aes(x=fracOnFront, y=frequency,
                                 col=alnDrv, label=alteration),
                     show.legend=F, alpha=.7) +
    theme_gray() +
    labs(x=sprintf("fraction of ||v||^2 on %s", s), y="frequency")

p3 <- 
    ggplot(allMutsPlot) +
    geom_line(data=tcl, mapping=aes(x=es, y=frequency, group=p), linetype=2) +
    geom_point(aes(x=normOnFront^2, y=frequency, col=alnDrv)) +
    geom_label_repel(data=filter(allMutsPlot,
                                 !isAln & isDrivingAlteration &
                                 log10qNorm < log10(qValThreshold) &
                                 frequency>minFreq),
                     mapping=aes(x=normOnFront^2, y=frequency,
                                 col=alnDrv, label=alteration),
                     show.legend=F, alpha=.7) +
    scale_x_log10() +
    theme_gray() +
    labs(x=sprintf("effect size on front ||v_on||^2", s), y="frequency")


plot_grid(p1, p2, p3, labels=toupper(letters[1:2]), ncol=3, align="v")
ggsave(sprintf("fracSSonFront_SS_freq_missedDrivers_%s_%s.pdf", altType, s),
       height=4, width=7*3)
## }

sapply(mutsAlongFront, nrow)

## sapply(overlap(drivers[["Santarius2010"]],
##                drivers[["Pereira2016"]]), length)
## sapply(overlap(drivers[["Santarius2010"]],
##                drivers[["NikZainal2016"]]), length)
## sapply(overlap(drivers[["Pereira2016"]],
##                drivers[["NikZainal2016"]]), length)

## if ( altType == "CNAs" ) {
##     universeGeneNames <- unique(sub("=.*$", "", copGeneNames))
## } else {
##     universeGeneNames <- unique(sub("=.*$", "", rownames(Mp)))
## }

## GOtab <- lapply(names(mutsAlongFront), function(s) {
##     geneList <- factor(as.integer(
##         universeGeneNames %in%
##         unique(gsub("=.*$", "", rownames(mutsAlongFront[[s]])))
##     ))
##     names(geneList) <- universeGeneNames;

##     GOdataBP <- new("topGOdata",
##                     description = s, ontology = "BP",
##                     allGenes = geneList, nodeSize = 10,
##                     annot = annFUN.org,
##                     mapping="org.Hs.eg.db",
##                     ID="symbol")
##     resultDef <- 
##         runTest(GOdataBP, algorithm = "weight01", statistic = "fisher")
##     allResBP <- GenTable(GOdataBP, Weighted=resultDef, topNodes=100)

##     GOdataMF <- new("topGOdata",
##                     description = s, ontology = "MF",
##                     allGenes = geneList, nodeSize = 10,
##                     annot = annFUN.org,
##                     mapping="org.Hs.eg.db",
##                     ID="symbol")
##     resultDef <- 
##         runTest(GOdataMF, algorithm = "weight01", statistic = "fisher")
##     allResMF <- GenTable(GOdataMF, Weighted=resultDef, topNodes=100)

##     allRes <- rbind(allResBP, allResMF)
##     allRes <- allRes[order(as.numeric(allRes[,"Weighted"])),]
##     head(allRes)
##     ## hist(as.numeric(allRes[,"Weighted"]))
##     return(allRes)
## })

## lapply(GOtab, head)

## ## For SNVs, too little mutations to make a difference.

## ## No big differences between enriched GO categories in the three sets
## ## of CNAs:
## ## poly(A) RNA binding
## ## mitochondrial translation elongation & termination
## ## mitosis

## ## Poly(A) mRNA binding
## intersect(annFUN.org(whichOnto="MF", mapping="org.Hs.eg.db",
##                      ID="symbol")[["GO:0044822"]],
##           unique(gsub("=.*$", "", rownames(mutsAlongFront[["Front"]]))))
## ## Translation regulation, splicing / NMD factors, DNA repair

##################################################
## 3D plot of samples and mutations

## Compute scalar product of muts with vectors pointing to
## different archetypes, to later select best aligned mutation
## with each archetype
i <- 1
Xproj[1:5,]
archProj <- t(sapply(1:nrow(arcsOrig), function(i) {
    m0 <- arcsOrig[i,] - geneExprAvg
                                        #- healthyProfile - avgE
    matrix(m0, nrow=1) %*% projMat
}))

## plot3d(t(apply(E, 2, function(x) { x + healthyProfile - geneExprAvg })) %*% projMat)
## spheres3d(archProj, radius=5)

hPt <-
    (- matrix(apply(E, 1, mean), nrow=1) ) %*% projMat 

s <- 10;
s <- 3;

## Set archetype colors
arcCols <- c("#4078fa", "#ffff00", "#e379e3", "#ca5e5e", "#28c928")
superArcs <-
    read_csv("~/work/cancerTaskAtlas/TCGA/ALL_UCSC/arcsOrig_genes.csv",
             col_names=unlist(read_csv("~/work/cancerTaskAtlas/TCGA/ALL_UCSC/geneNamesAfterExprFiltering.list",
                                       col_names=F)))

commonGenes <- intersect(colnames(arcsOrig), colnames(superArcs))
##################################################
## Compute distances to super archetypes, pick smallest
matchSAtab <-
    sapply(1:nrow(arcsOrig), function(i) {
        apply(superArcs[,commonGenes], 1, function(x) {
            sqrt(sum(( (arcsOrig[i,commonGenes] - geneExprAvg[commonGenes]) - x )^2))
        })
    })
## apply(matchSAtab, 2, function(x) { arcCols[x^2 / sum(x^2) < .1] })
tCols <- apply(matchSAtab, 2, function(x) { arcCols[which.min(x^2 / sum(x^2))] })

##################################################
## Pick least square permutation
enumPerms <- function(m=4, S=1:5) {
    pLists <- list();
    if ( m > 0 ) {
        for (i in 1:length(S)) {
            subLists <- enumPerms(m-1, S[-i])
            ## browser();
            if ( length(subLists) > 0 ) {
                for (j in 1:length(subLists)) {
                    pLists[[length(pLists)+1]] <- c(S[i], subLists[[j]])
                }
            } else {
                pLists[[length(pLists)+1]] <- S[i]
            }
        }
    }
    return(pLists);
}
allPerms <- enumPerms(m=nrow(archProj), S=1:nrow(superArcs))
p <- allPerms[[1]]
mAO <- arcsOrig[,commonGenes] - geneExprAvg[commonGenes];
mSA <- superArcs[,commonGenes]
if ( file.exists("superArcs_to_arcs_allDists.rda") ) {
    load("superArcs_to_arcs_allDists.rda");
    load("topPerm.rda")
} else {
    allDists <- sapply(allPerms, function(p) {
        sqrt(sum(
        ( mAO - mSA[p,] )^2
        ))
    })
    save(allDists, file="superArcs_to_arcs_allDists.rda")
    p53arc <- which.max(Xproj["TP53=1",] %*% vRefs)
    doesP53match <- sapply(allPerms, function(x) { x[p53arc] == 3 || x[p53arc] == 4 })
    doesHealthyMatch <- sapply(allPerms, function(x) { x[healthyArch] == 5})
    
    plot(allDists, pch=20, col=(doesHealthyMatch & doesP53match) + 2)
    plot(sort(allDists), pch=20, col=(doesHealthyMatch&doesP53match)[order(allDists)] + 2)
    
    sapply(order(allDists)[(doesHealthyMatch&doesP53match)[order(allDists)]][1:5],
           function(i) { allPerms[[i]] })
    topPerm <- order(allDists)[(doesHealthyMatch&doesP53match)[order(allDists)]][1]
    save(topPerm, file="topPerm.rda")
}

tCols <- arcCols[allPerms[[topPerm]]]
setdiff(1:nrow(archProj), healthyArch)
vRefs <-
    sapply(1:nrow(archProj), function(i) {
        ## vRef <- archProj[i,] - archProj[healthyArch,]
        vRef <- archProj[i,]
        vRef <- vRef / sqrt(sum(vRef^2))
        return(vRef)
    })
## sum(vRefs[,2]^2)

mutOrientToArchs <-
    Xproj[unlist(allMuts %>% 
                 select(alteration)),] %*%
    vRefs
mutOrientToArchs <-
    Xproj[unlist(allMuts %>% filter(isAln & isDrivingAlteration) %>% 
                 select(alteration)),] %*%
    vRefs
mutOrientToArchs <-
    t(apply(mutOrientToArchs, 1, function(x) {
        x[x < max(x)] <- 0
        return(x)
    }))
sapply(1:nrow(archProj), function(i) {
    list(rev(sort(mutOrientToArchs[,i]))[1:5])
})
m <- "GATA3=1"
m <- "BRAF=1"
m <- "SMAD4=-1";
m <- "ERBB2=1";
m <- "BRCA1=1";
m <- "CCND3=2"
i <- 1
i <- 2
i <- 3

save.image(sprintf("alt3Dplots_%s.rda", altType))
## load(sprintf("alt3Dplots_%s.rda", altType))

sapply(1:nrow(archProj), function(i) {
    m <- names(rev(sort(mutOrientToArchs[,i]))[1])
    if ( mutOrientToArchs[m,i] == 0 ) { return(); }
    ## x <- Xproj[m,1:3]
    ## v <- x - hPt;
    v <- Xproj[m,1:3]

    ptCol <- rep("grey", ncol(Mp))
    ptCol[Mp[m,]==1] <- "red"

    open3d()
    plot3d(posInPCspace,
           ## alpha=.35, col="grey",
           col=ptCol, type="p", size=4, shininess=100, alpha=.5,
           box=F, axes=F,
           ## xlab="PC1", ylab="PC2", zlab="PC3",
           xlab="", ylab="", zlab="",
           xlim=range(archProj[,1]),
           ylim=range(archProj[,2]),
           zlim=range(archProj[,3]))
    ## spheres3d(hPt, col="green", radius=5)
    spheres3d(archProj, radius=10, col=tCols)
    ## text3d(archProj[,1], archProj[,2], archProj[,3],
    ##        allPerms[[topPerm]],
    ##        ## 1:nrow(archProj),
    ##        adj=2.5)
    sapply(1:(nrow(archProj)-1), function(i) {
        sapply(seq(i+1, nrow(archProj)), function(j) {
            segments3d(archProj[c(i,j),1], archProj[c(i,j),2],
                       archProj[c(i,j),3], col="grey")
        })
    })
    
    arrow3d(c(0, 0, 0), v*s, width=.1)
    cat(paste(m, "\n"))
    ## text3d(s * v[1], s * v[2], s * v[3],
    ##        sub("_CNA$", "", sub("=1$", "", m)),
    ##        adj=c(.5, 0))
})

## Stain MSI-H samples
discClin <- read.table("discreteClinicalData_reOrdered.tsv",
                       as.is=T, sep="\t", h=T)

unique(discClin[,"CDE_ID_3226963"])
ptCol <- (discClin[,"CDE_ID_3226963"] == "MSI-H") + 1
ptCol <- (discClin[,"CDE_ID_3226963"] == "MSS") + 1
ptCol <- (discClin[,"CDE_ID_3226963"] == "MSI-L") + 1
## ptCol <- as.numeric(as.factor(discClin[,"CDE_ID_3226963"])) - 2

plot3d(posInPCspace,
       ## alpha=.35, col="grey",
       col=ptCol, type="p", size=4, shininess=100, alpha=1,
       box=F, axes=F,
       ## xlab="PC1", ylab="PC2", zlab="PC3",
       xlab="", ylab="", zlab="",
       xlim=range(archProj[,1]),
       ylim=range(archProj[,2]),
       zlim=range(archProj[,3]))
## spheres3d(hPt, col="green", radius=5)
spheres3d(archProj, radius=10, col=tCols)
## text3d(archProj[,1], archProj[,2], archProj[,3],
##        allPerms[[topPerm]],
##        ## 1:nrow(archProj),
##        adj=2.5)
sapply(1:(nrow(archProj)-1), function(i) {
    sapply(seq(i+1, nrow(archProj)), function(j) {
        segments3d(archProj[c(i,j),1], archProj[c(i,j),2],
                   archProj[c(i,j),3], col="grey")
    })
})
