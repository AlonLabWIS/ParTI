# Init --------------------------------------------------------------------

rm(list=ls()); gc();
source("~/work/jeanLib.R")
library(MASS)
library(rgl)
library(ade4)
library(parallel)
library(fdrtool)
library(Hmisc)
library(cowplot)

library(tidyverse)
library(ggrepel)
library(stringr)

altType <- "SNVs"; #SNVs, CNAs
## altType <- "CNAs"; #SNVs, CNAs
isMetabric <- length(grep("metabric", getwd())) == 1
cancerType <- getwd() %>%
    str_replace("^.*/", "") %>% str_replace("_UCSC$", "")

# load data ---------------------------------------------------------------

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

## ## Centering is commutative:
## dim(geneExpression)
## length(apply(geneExpression, 1, function(x) { mean(x) })) #for each sample
## length(apply(geneExpression, 2, function(x) { mean(x) })) #for each gene
## dblCentered <-
##     apply(
##         t(apply(geneExpression, 1, function(x) { x - mean(x) })),
##         2, function(x) { x - mean(x) })
## dblCentered2 <-
##     t(apply(
##         apply(geneExpression, 2, function(x) { x - mean(x) }),
##         1, function(x) { x - mean(x) }))
## dim(dblCentered)
## dim(dblCentered2)
## sum(abs(dblCentered - dblCentered2))
## ## End

## Diagnostics
## cols <- rep(1, nrow(geneExpression))
## cols[270] <- 2

## par(mfrow=c(2,3))
## plot(apply(geneExpression, 1, function(x) { sum(exp(x)) }),
##      pch=20, col=cols, cex=cols)
## plot(apply(geneExpression, 1, function(x) { sum(x) }),
##      pch=20, col=cols, cex=cols)
## ## Renormlize gene expression
## ## geneExpression <- t(apply(geneExpression, 1, function(x) { log(exp(x) * 1e9 / sum(exp(x))) }))

## plot(apply(geneExpression, 1, mean), apply(geneExpression, 1, sd),
##      pch=20, col=cols, cex=cols)
## plot(apply(geneExpression, 1, median),
##      apply(geneExpression, 1, function(x) { quantile(x, .75) - quantile(x, .25) }),
##      pch=20, col=cols, cex=cols)
## plot(apply(geneExpression, 2, mean), apply(geneExpression, 2, sd))

## par(mfrow=c(3,3))
## sapply(names(which(apply(geneExpression, 2, sd)>4)), function(g) {
##     hist(geneExpression[,g], main=g)
## })
## hist(geneExpression[,"TP53"], main="TP53")

## End diagnostics

which(apply(is.na(geneExpression), 1, sum)>0)
geneExprAvg <- apply(geneExpression, 2, mean)

discClin <- read.table("discreteClinicalData_reOrdered.tsv",
                       as.is=T, sep="\t", h=T)
if ( isMetabric ) {
    ## ## FIXME --- also look for samples with 0 mutations, not just NaN?
    ## samplesMutProfiled <- read.table("case_lists_sequenced.list", h=F, as.is=T)[,1]
    ## patientIDs <- read.table("patientIDs.list", h=F, as.is=T)[,1]
    ## isMutProfiled <- sapply(patientIDs, function(p) { p %in% samplesMutProfiled })
    ## isHealthy <-
    ##     apply(is.nan(mutations), 1, sum) == ncol(mutations) & isMutProfiled
    ## healthyArch <- 2;
    isPT <- rep(T, nrow(discClin))
} else {
    ## isHealthy <- discClin[,"sample_type"] == "Solid Tissue Normal"
    ## healthyArch <-
    ##     as.numeric(
    ##         read_csv("clinicalEnrichment_discrete_significant.csv") %>%
    ##         filter(`Feature Name` == 'sample type: Solid Tissue Normal') %>%
    ##         select(`archetype #`)
    ##     )
    isPT <- discClin[,"sample_type"] == "Primary Tumor"
}
## rm(discClin)

mutations <- read.csv("mutMatrix_reOrdered_booleanized_justData.csv",
                      as.is=T, h=F)
mutGeneNames <-
    read.table("mutMatrix_reOrdered_booleanized_geneNames.list",
               h=F, as.is=T)[,1]
colnames(mutations) <- mutGeneNames;
rm(mutGeneNames)
dim(mutations)
mutations <- as.matrix(mutations)[isPT,]

if ( altType == "CNAs" ) {
    CNAs <- read_delim("copMatrix_reOrdered_booleanized.tsv", delim="\t", col_names=T)
    CNAs <- as.data.frame(CNAs)
    rownames(CNAs) <- CNAs[,1]
    CNAs <- CNAs[isPT,-1]

    copGeneNames <- colnames(CNAs)
    x <- copGeneNames[1]
    copTab <- t(sapply(copGeneNames, function(x) {
        strsplit(x, "=")[[1]]
    }))
}
length(apply(geneExpression, 2, function(x) { 1 }))
## gE0 <- apply(geneExpression, 1, function(x) { x - mean(x) })
## rm(geneExpression)
## gE00 <- apply(gE0, 2, function(x) { x - geneExprAvg } )
## gE00 is centered on healthyProfile

gE00 <- t(apply(geneExpression[isPT,], 2, function(x) { x - mean(x) }))
## gE00 <-
##     t(apply(
##         t(apply(geneExpression[isPT,], 1, function(x) { x - mean(x) })),
##         2, function(x) { x - mean(x) }))

E <- gE00
rm(gE00, tmp, geneExpression); gc();

## Mp <- rbind(M, t(CNAs))
if ( altType == "SNVs" ) {
    Mp <- t(mutations);
} else if ( altType == "CNAs" ) {
    Mp <- t(CNAs);
} else {
    stop("Unknown altType!");
}

rm(mutations, CNAs); gc();

arcsOrig <-
  as.matrix(
    read_csv("arcsOrig_genes.csv",
             col_names=unlist(read_csv("geneNamesAfterExprFiltering.list", col_names=F))))
nArchs <- nrow(arcsOrig);
nPCs <- nArchs - 1;
if ( nPCs == 2 ) {
    nPCs <- 3; #so we can visualize with the same 3D tools
}

save.image(sprintf("learnMutEffects_%s.rda", altType))

# prepare mutations data -------------------------------------------------------

# load(sprintf("learnMutEffects_%s.rda", altType))

dim(E)
dim(Mp)

sel <- 1:nrow(Mp);

## At least 10 examples so vectors aren't too noisy

## commonMuts <-
##     intersect(sel,
##               setdiff(1:nrow(Mp),
##                       c(grep("=NaN", rownames(Mp)),
##                         grep("=NA", rownames(Mp)),
##                         grep("=0", rownames(Mp)))))
commonMuts <-
    intersect(sel,
              setdiff(1:nrow(Mp),
                      c(grep("=NaN", rownames(Mp)),
                        grep("=NA", rownames(Mp)) )))
rownames(Mp)[commonMuts][1:5]
altFreq <- apply(Mp[commonMuts,], 1, function(x) { sum(x, na.rm=T) })
hist((altFreq), 20)

## number of unprofiled mutations per sample
if ( isMetabric ) {
    mutUnprofiled <- apply(is.nan(Mp), 2, sum) == nrow(Mp)
} else {
    mutUnprofiled <-
        apply(Mp[grep("=NaN", rownames(Mp)),], 2, sum) == length(grep("=NaN", rownames(Mp)))
}
if ( length(grep("ALL_UCSC$", getwd())) == 1 ) {
    mutUnprofiled <- rep(F, ncol(Mp))
}

## number of mutations found per sample
## plot(ecdf(apply(Mp[commonMuts,], 2, sum))); abline(v=400)
plot(ecdf(apply(Mp[commonMuts,], 2, function(x) { sum(x, na.rm=T) }))); abline(v=400)

if ( cancerType == "COAD" && altType == "SNVs" ) {
    ggplot(tibble(nMuts=apply(Mp[commonMuts,], 2, sum),
                  MS=discClin$CDE_ID_3226963[isPT])) +
        geom_boxplot(aes(MS, nMuts))
    hyperMutated <- apply(Mp[commonMuts,], 2, sum) >= 400
    table(discClin$CDE_ID_3226963[isPT], hyperMutated)
} else {
    hyperMutated <- rep(F, ncol(Mp))
}

# library(gplots)
# MpShow <- Mp[commonMuts,!mutUnprofiled];
# rownames(MpShow) <-
#     sapply(rownames(MpShow), function(x) {
#         ifelse(x %in% allDrivers, x, "")
#     })
# ## Number of samples in which mutations are found
# mutFreq <- apply(Mp[commonMuts,!mutUnprofiled & !hyperMutated], 1, sum);
# rev(sort(mutFreq))[1:5]
# 
# sort(mutFreq[allDrivers])
# 
# plot(ecdf(mutFreq));
# abline(v=5);
# sum(mutFreq>=5)
# dim(MpShow)
# png("mutClust.png", height=1600, width=1200, pointsize=32)
# heatmap.2(MpShow[mutFreq>=5,], scale="none", trace="none",
#           distfun=function(x) { dist(x, method="manhattan") },
#           hclustfun=function(x) { hclust(x, method="ward.D2") })
# dev.off();



## Do some mutations happen in exactly the same samples?
distMat <- dist(t(Mp[commonMuts,!mutUnprofiled]))
sum(distMat == 0) #almost 3000 identical pairs! This may explain why
                  #non-driver non cancer genes align so well with
                  #front: they behave like drivers!
## There are no frequency 0 mutations -> identical profiles came from
## the 25% of samples with no profiled mutations
## mutInSamples <-
##     apply(Mp[commonMuts,], 1, function(x) {
##         paste(which(as.logical(x)), collapse=",")
##     })
## mutInSamples[1:5]
## rev(sort(table(mutInSamples)))[1:10]
## enframe(table(mutInSamples)) %>% arrange(desc(value))
## which(mutInSamples == "62") ## these mutations occur only in sample '62'
## which(Mp[commonMuts,]["CTSG=1",] == 1)
## which(Mp[commonMuts,]["EMB=1",] == 1)
## sum(Mp[commonMuts,62]) ## This sample has 2715 mutations!

## which(mutInSamples == "142") ## these mutations occur only in sample '142'
## which(Mp[commonMuts,]["ARL2=1",] == 1)
## which(Mp[commonMuts,]["CMAS=1",] == 1)
## sum(Mp[commonMuts,142]) ## This sample has 1426 mutations!

## summary(apply(Mp[commonMuts,], 2, sum)) 
## sum(apply(Mp[commonMuts,], 2, sum) == 0) / ncol(Mp)
## ## 25% of samples have no mutations... because they haven't been
## ## profiled?

plot(ecdf(apply(Mp[commonMuts,!mutUnprofiled], 2, sum))); abline(v=400)

altFreq["PTEN=-1"]; altFreq["PTEN=1"]; altFreq["TP53=-1"]; altFreq["TP53=1"]
minFrequency <- 10; ## setting 1 here makes it difficult to pass FDR
abline(v=(minFrequency), lty=2)

sum(apply(Mp[commonMuts,!hyperMutated], 1,
          function(x) {sum(x, na.rm=T)}) >= minFrequency) / nrow(Mp)
commonMuts <-
    intersect(commonMuts,
              which(apply(Mp[,!hyperMutated], 1, function(x) { sum(x, na.rm=T) }) >= minFrequency))
length(commonMuts)
rownames(Mp)[commonMuts][1:5]
## Minimal number of control (copy-unaltered) genes:
## min(apply(Mp[grep("=0$", rownames(Mp)),], 1, sum))
if ( altType == "CNAs" ) {
    ## Remove copy-unaltered alterations
    commonMuts <- setdiff(commonMuts, grep("=0$", rownames(Mp)))
}

## MpP <- apply(Mp[commonMuts,], 1, function(x) { x / sum(x, na.rm=T) })
MpP <- apply(Mp[commonMuts,!mutUnprofiled], 1, function(x) { x / sum(x) })

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
    ## ctrlMuts <- gsub('=.*$', '=0', rownames(Mp)[commonMuts])
    ## MpN <- apply(Mp[ctrlMuts,], 1, function(x) { x / sum(x) })
    ctrlMuts <- gsub('=.*$', '=0', rownames(Mp)[commonMuts])
    MpN <- apply(Mp[ctrlMuts, !mutUnprofiled], 1, function(x) { x / sum(x) })
} else {
    MpN <- apply(1 - Mp[commonMuts, !mutUnprofiled], 1, function(x) { x / sum(x, na.rm=T) })
}

MpP[is.nan(MpP) | is.na(MpP)] <- 0
MpN[is.nan(MpN) | is.na(MpN)] <- 0

Xd <- E[,!mutUnprofiled] %*% ( MpP - MpN )
Xd[1:5,1:5]
X <- Xd;
rm(Xd, MpP, MpN); gc();

## These are the mutations with strongest impact on gene expression
rev(sort(apply(X, 2, sd)))[1:10]

mutWithStrongestEffect <- names(rev(sort(apply(X, 2, sd)))[1:10])

dudi2 <- dudi.pca(t(E), scale=F, scannf=F, nf=nPCs)

## avgE <- apply(E, 1, mean)

## save.image(sprintf("learnMutEffects2_%s.rda", altType))
## load(sprintf("learnMutEffects2_%s.rda", altType))

projMat <- as.matrix(dudi2$c1);
posInPCspace <- as.matrix(dudi2$li);
EPCAeig <- dudi2$eig;
rm(dudi2); gc()

## X is a vector. No need to rewrite it in order to project it!
## X0 <- apply(X, 2, function(x) { x - avgE })
## X0 <- X; rm(X);
Xproj <- t(X) %*% projMat;
save(Xproj, file=sprintf("Xproj_%s.rda", altType))

Xproj[mutWithStrongestEffect,]
XbackProj <- projMat %*% t( t(X) %*% projMat )

plot(X[,mutWithStrongestEffect[10]], XbackProj[,mutWithStrongestEffect[10]]);
plot(X[,mutWithStrongestEffect[10]] - mean(X[,mutWithStrongestEffect[10]]),
     XbackProj[,mutWithStrongestEffect[10]]);
abline(0,1, col="grey");
plot(X[,50], XbackProj[,50]);
abline(0,1, col="grey");

plot(X[,50] - mean(X[,50]), XbackProj[,50]);
abline(0,1, col="grey");

rm(XbackProj)

save.image(sprintf("learnMutEffects3_%s.rda", altType))

# compare alignment of mutations to front(s) ------------------------------

# load(sprintf("learnMutEffects3_%s.rda", altType))

## 3D "Pareto front" subspace
subspaces <- list()
subspaces[["Front"]] <- projMat[,1:(nArchs-1)];

## ## healthy -> cancer axis
## tmp <-
##     as.numeric(
##         apply(as.matrix(arcsOrig)[-healthyArch,], 2, mean) -
##         as.matrix(arcsOrig)[healthyArch,]
##     )
## tmp <- tmp / sqrt(sum(tmp^2))
## sum(tmp^2)
## summary(tmp)
## summary(projMat)
## subspace <- matrix(tmp, ncol=1)
## subspaces[["HealthyCancer"]] <- subspace;

## ## cancer archetypes subspace
## is <- setdiff(1:nrow(arcsOrig), healthyArch)
## tmp <- cbind(as.numeric(arcsOrig[is[2],] - arcsOrig[is[1],]),
##              as.numeric(arcsOrig[is[3],] - arcsOrig[is[1],]))
## tmp[,1] <- tmp[,1] / sqrt(sum(tmp[,1]^2))
## tmp[,2] <- tmp[,2] - tmp[,1] * sum(tmp[,1] * tmp[,2])
## tmp[,2] <- tmp[,2] / sqrt(sum(tmp[,2]^2))
## tmp[1:5,]
## summary(tmp)
## sum(tmp[,1] * tmp[,2])
## subspace <- tmp
## subspaces[["CancerPlane"]] <- subspace;

subspace <- subspaces[["Front"]]

## nPoss <-
##     round(exp(seq(log(10), log(max(apply(Mp[commonMuts,], 1, sum))), len=12)))
nPoss <-
    round(exp(seq(log(10),
                  log(floor(.99 * ncol(Mp[commonMuts,]))),
                  len=12)))
i <- 1

save(E, nPoss, subspaces, file=sprintf("SSonDistInput_%s.rda", altType))

## Run SSonDist.R (in a different R session for memory efficiency to compute SSonDistRevCumApprox.rda
## load(sprintf("SSonDistRevCumApprox_%s.rda", altType))
## Results are the same for CNAs and SNVs:
load(sprintf("SSonDistRevCumApprox_%s.rda", "SNVs"))

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
    dplyr::summarize(n=n()) %>% mutate(uniqueLocus=(n<=10))
inner_join(geneChr, nLoci) %>% arrange(desc(n), gene) %>%
    select(gene,n) %>% distinct()  %>% print(n=100)

geneChrU <- inner_join(geneChr, nLoci) %>% filter(uniqueLocus)
## geneChr[,"chr"] <- as.numeric(geneChr[,"chr"])
chrLens <-
    group_by(geneChr, chr) %>%
    dplyr::summarize(chrLen=max(end), chrNum=mean(chrNum)) %>%
    arrange(chrNum) %>%
    mutate(offsets=lag(cumsum(as.numeric(chrLen)))) %>%
    replace_na(replace=list(offsets=0))
print(chrLens, n=100)

genomeLen <- as.numeric(
    summarize(chrLens, sum(as.numeric(chrLen))))

universeGeneNames <- unique(sub("=.*$", "", rownames(Mp)))
## qValThreshold <- .01;
qValThreshold <- .1; ## to get a few aligned CNAs

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
## load(sprintf("SSonDistSamples_%s_%s.rda", altType, s))
load(sprintf("SSonDistSamples_%s_%s.rda", "SNVs", s))
rratios <- ratios; rm(ratios);


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
    ## browser();
    which( sum(x, na.rm=T) <= nPoss )[1] - 1
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

getBHcutOff <- function(p, fdr=.1) {
    ps <- sort(p);
    cutOffIdx <- sum(ps < fdr*(1:length(ps))/length(ps))
    if ( cutOffIdx == 0 ) {
        return(0)
    }  else {
        ps[cutOffIdx]
    }
}

## hist(10^log10pRatios)
## Bonferroni correction
## log10q <- log10p + log10(length(onf))
## log10q[log10q>0] <- 0
## hist(log10q)
nBHSS <- sum(10^sort(log10pSS) <
             qValThreshold * (1:length(log10pSS))/length(log10pSS)
             )
nBHSSonFront <- sum(10^sort(log10pSSonFront) <
                    qValThreshold *
                    (1:length(log10pSSonFront))/length(log10pSSonFront)
                    )
nBHratios <- sum(
    10^sort(log10pRatios) <
    qValThreshold * (1:length(log10pRatios)) /
    length(log10pRatios)
)

log10qSS <- log10(fdrtool(10^log10pSS, statistic="pvalue")$qval)
## hist(10^log10qSS, seq(0,1,by=.1), col="lightblue")
## hist(10^log10pSS, seq(0,1,by=.1), add=T)
## plot(log10pSS, log10qSS); abline(0,1)

log10qSSonFront <- log10(fdrtool(10^log10pSSonFront,
                                 statistic="pvalue")$qval)
## hist(10^log10qSSonFront, seq(0,1,by=.1), col="lightblue")
## hist(10^log10pSSonFront, seq(0,1,by=.1), add=T)
## plot(log10pSSonFront, log10qSSonFront); abline(0,1)

## log10qRatios <- log10(fdrtool(10^log10pRatios, statistic="pvalue")$qval)
## if ( diff(range(log10qRatios)) < 1 ) {
##     cat("FDRs don't vary; fdrtool probably failed. Falling back on BH.\n")
    log10qRatios <- log10pRatios;
    log10qRatios[order(log10pRatios)[seq(nBHratios+1, length(log10pRatios))]] <- 0
## }
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

## minFreq <- .05 * ncol(Mp)
minFreq <- min(nPoss);
## minFreq <- 5;

allMuts <- tibble(
    alteration=names(ratios),
    fracOnFront=ratios,
    log10qFracOnFront=log10qRatios,
    normOnFront=normsOnFront,
    log10qNormOnFront=log10qSSonFront,
    norm=norms,
    log10qNorm=log10qSS,
    frequency=apply(Mp[commonMuts,!mutUnprofiled], 1, sum),
    freqNHM=apply(Mp[commonMuts,!hyperMutated&!mutUnprofiled], 1, sum)) %>%
    mutate(gene=str_replace(alteration, "=.*$",""));
allMuts <-
    mutate(allMuts,
           ## isAln=log10qNormOnFront<log10(qValThreshold)
           isAln=log10qFracOnFront<log10(qValThreshold)
           ## isAln=log10qFracOnFront<log10(qValThreshold) &
           ##     log10qNorm<log10(qValThreshold)
           )
if ( isMetabric && altType == "SNVs") {
    allMuts <-
        allMuts %>% mutate(alteration=paste(alteration, "=1", sep=""))
}
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

allMuts %>% filter(grepl("IDH1$", gene))

allMuts %>% filter(grepl("APC$", gene))

if ( isMetabric && altType == "CNAs" ) {
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
if ( altType == "CNAs" ) {
    allDriversList <-
        lapply(names(drivers), function(x) {
            c(sprintf("%s=1", as.data.frame(filter(drivers[[x]], ONC))[,"gene"]),
              sprintf("%s=2", as.data.frame(filter(drivers[[x]], ONC))[,"gene"]),
              sprintf("%s=-1", as.data.frame(filter(drivers[[x]], TSG))[,"gene"]),
              sprintf("%s=-2", as.data.frame(filter(drivers[[x]], TSG))[,"gene"]))
        })
} else {
    allDriversList <-
        lapply(names(drivers), function(x) {
            sprintf("%s=1",
                    as.data.frame(filter(drivers[[x]],
                                         ONC | TSG))[,"gene"])
        })
}
names(allDriversList) <- names(drivers)
allDrivers <- unique(unlist(allDriversList))

allMuts <- allMuts %>% mutate(isDrivingAlteration=alteration %in% allDrivers)
if ( ! isMetabric || ( isMetabric & altType == "SNVs" ) ) {
    allMuts <- allMuts %>% mutate(isTopGene=T, isRndGene=F)
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

# predict drivers by aln to front -----------------------------------------

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

cancerGenes <-
    (read_tsv("~/work/cancerTaskAtlas/SandersNatGenetics2013/cancer_genes.list",
              col_names=F) %>% unlist %>% as.character)

library(forcats)
cmpAlns <-
    allMuts %>%
    select(gene, alteration, fracOnFront, isDrivingAlteration, freqNHM, frequency);
cmpAlns <-
    cmpAlns %>% 
    mutate(isCancerGene=map_lgl(cmpAlns[,"gene"] %>% unlist %>% as.character,
                                ~ . %in% cancerGenes));
newFac <-
    factor(
        pmap_chr(cmpAlns %>%
                 select(isDrivingAlteration, isCancerGene),
                 function(isDrivingAlteration, isCancerGene) {
                     if ( isDrivingAlteration ) {
                         return("driver gene");
                     } else {
                         if ( isCancerGene ) {
                             return("cancer gene");
                         } else {
                             return("other gene");
                         }
                     }
                 }),
        levels=c("driver gene", "cancer gene", "other gene", "shuffled"))
cmpAlns <- cmpAlns %>% mutate(altType=newFac);
if ( sum(hyperMutated) > 10 ) {
    ggplot(cmpAlns) +
        geom_point(aes(frequency-freqNHM, freqNHM, color=altType)) +
        theme_gray() + scale_x_log10() + scale_y_log10();
}

## How many non hyper mutated sample for each hyper mutated sample?
nNHM <- sum(!mutUnprofiled & !hyperMutated)
nHM <- sum(!mutUnprofiled & hyperMutated)
cmpAlns <- 
    cmpAlns %>%
    ## filter(frequency>=minFreq) %>% 
    filter((altType=="driver gene") |
           (altType=="cancer gene") |
           (altType=="other gene" &
            (freqNHM == frequency |
             freqNHM / nNHM >= (frequency-freqNHM) / nHM))
           ) %>%
    ## & frequency-freqNHM>freqNHM*5)) %>%
    ## altType=="other gene") %>%
    filter(frequency>minFreq) %>% 
    select(-frequency,-freqNHM) %>% 
    ## Add shuffle control mutations
    bind_rows(tibble(fracOnFront=rratios,
                     altType=factor(rep("shuffled", length(rratios)),
                                    levels=levels(newFac))
                     ));
## cmpAlns <-
##     cmpAlns %>%
##     left_join(
##         cmpAlns %>% group_by(altType) %>% summarize(n=n()) %>%
##         mutate(altLabel=sprintf("%s\n(n=%d)", altType, n)) %>%
##         select(-n)
##     )
## Compare to shuffled
if ( altType == "SNVs" ) {
    pWilcox <-
        sapply(c("driver gene", "cancer gene",
                 "other gene", "shuffled"), function(lab) {
            wilcox.test(rratios,
                        cmpAlns %>% filter(altType==!!lab) %>%
                        select(fracOnFront) %>% unlist, 
                        alternative="less")$p.value
                 }) %>% enframe(name="altType", value="p") %>%
        mutate(p=sprintf("p = %.1e", p)) %>%
        mutate(altType=as_factor(altType, levels(cmpAlns["altType"])))
} else {
    labs <- c("driver gene", "cancer gene", "other gene", "shuffled")
    pWilcox <- 
        sapply(labs[1:3], function(lab) {
            i <- which(labs == lab)
            X <- cmpAlns %>% filter(altType==!!labs[i]) %>%
                select(fracOnFront) %>% unlist
            if ( labs[i+1] == "shuffled" ) {
                Y <- rratios
            } else {
                Y <- cmpAlns %>% filter(altType==!!labs[i+1]) %>%
                    select(fracOnFront) %>% unlist
            }
            wilcox.test(X, Y, alternative="greater")$p.value
        }) %>% enframe(name="altType", value="p") %>%
        mutate(p=sprintf("p = %.1e", p)) %>%
        mutate(altType=as_factor(altType, levels(cmpAlns["altType"])))
}

## cmpAlns %>% inner_join(pWilcox) %>% unite(altType, p, col="altType", sep="\n")
ggplot(cmpAlns) +
    geom_boxplot(aes(x=altType, y=100*fracOnFront), width=.4) +
    ## geom_crossbar(aes(x=altType, y=y, ymin=ymin, ymax=ymax), width=.4, fill="white",
    ##               data=cmpAlns %>% group_by(altType) %>%
    ##                   summarize(y=median(100*fracOnFront),
    ##                             ymin=quantile(100*fracOnFront, .25),
    ##                             ymax=quantile(100*fracOnFront, .75))
    ##               ) +
    ## geom_jitter(aes(x=altType, y=100*fracOnFront, color="#AA3333"), width=.1, 
    ##             data=cmpAlns %>% filter(altType != "shuffled")) +
    geom_hline(aes(yintercept=100*sum(EPCAeig[1:(nArchs-1)])/sum(EPCAeig)),
               linetype=2) +
    labs(x="Type of alteration", y="% alignment with Pareto front") +
    theme_gray() + theme(legend.position="none")
ggsave(sprintf("%s_alnWithFront.pdf", altType), height=2.5, width=3);

## The same, in degrees
toFlash <- allMuts %>%
    ## filter(isAln & isDrivingAlteration) %>%
    filter(isDrivingAlteration) %>%
    arrange(desc(isAln), desc(frequency)) %>%
    slice(1:(3+3)) %>% 
    select(alteration) %>%
    inner_join(cmpAlns)
if ( altType == "SNVs" ) {
    toFlash <- toFlash %>% mutate(alteration=gene)
} else if ( altType == "CNAs" ) {
    toFlash <- toFlash %>%
        mutate(alteration=str_replace_all(alteration,
                                          c("=-2"="(-2)",
                                            "=-1"="(-1)",
                                            "=1"="(+1)",
                                            "=2"="(+2)")))
}
## if ( altType == "SNVs" ) {
##     toFlash <- toFlash %>%
##         transmute(gene=str_replace(alteration, "=.*$", ""))
## }

## Downsample random genes to as many genes as driver genes for visualization
nDrivers <- cmpAlns %>% filter(altType == "driver gene") %>% nrow
## bind_rows(cmpAlns %>% filter(altType != "shuffled"),
##                  cmpAlns %>%
##                  filter(altType == "shuffled") %>%
##                  .[sample(1:length(rratios),
##                           ## 100),])) +
##                           nDrivers),]
##                  )

ggplot(bind_rows(
    cmpAlns %>% filter(altType != "shuffled"),
    cmpAlns %>%
    filter(altType == "shuffled") %>%
    arrange(fracOnFront) %>% 
    .[seq(1, length(rratios),
          len=ifelse(nDrivers<25, 25, nDrivers)),] )) +
    geom_boxplot(aes(x=altType, y=90*(1-fracOnFront)), width=.4) +
    geom_hline(aes(yintercept=90*(1-sum(EPCAeig[1:(nArchs-1)])/sum(EPCAeig))),
               linetype=2) +
    geom_point(aes(x=altType, y=90*(1-fracOnFront)),
               color="black", fill="white", shape=21, size=2,
               data=toFlash) +
    geom_text_repel(aes(x=altType, y=90*(1-fracOnFront), label=alteration),
                     size=3, data=toFlash) +
    geom_label(aes(x=altType, y=87.5, label=p), data=pWilcox) +
    labs(x="Type of alteration",
         y=expression("Angle to Pareto front (degrees)")) +
    theme_gray() + theme(legend.position="none")
## Main figure
## ggsave(sprintf("%s_alnWithFront_deg.pdf", altType), height=2.5, width=3);
## Suppl figure
ggsave(sprintf("%s_alnWithFront_deg.pdf", altType), height=4, width=5);

## Some other tests
length(unlist(cmpAlns %>% filter(altType=="driver gene") %>% select(fracOnFront)))
length(unlist(cmpAlns %>% filter(altType=="cancer gene") %>% select(fracOnFront)))
wilcox.test(unlist(cmpAlns %>% filter(altType=="driver gene") %>% select(fracOnFront)),
            unlist(cmpAlns %>% filter(altType=="cancer gene") %>% select(fracOnFront)),
            alternative="greater")
length(unlist(cmpAlns %>% filter(altType=="other gene") %>% select(fracOnFront)))
wilcox.test(unlist(cmpAlns %>% filter(altType=="driver gene") %>% select(fracOnFront)),
            unlist(cmpAlns %>% filter(altType=="other gene") %>% select(fracOnFront)),
            alternative="greater")

ggplot(allMuts) +
    geom_point(aes(x=frequency, y=normOnFront,
                   group=isDrivingAlteration, col=isDrivingAlteration)) +
    scale_x_log10()

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
                                          filter(frequency>=minFreq) %>% 
                                          filter(glm>=qCutOff) %>%
                                          select(isDrivingAlteration, log10qNormOnFront) %>% 
                                          summarize(specificity=mean(isDrivingAlteration))),
                   sensitivity=as.numeric(allMutsGlm %>%
                                          filter(frequency>=minFreq) %>% 
                                          filter(isDrivingAlteration) %>%
                                          transmute(caught=glm>=qCutOff) %>% 
                                          summarize(sensitivity=mean(caught))))
             }))

pdf(sprintf("ROCanalysis_%s.pdf", altType), height=3, width=3*3)
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

driversIdx <- names(drivers)[1]
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
        if ( nrow(filter(allMuts, isRndGene)) > 0 ) {
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
        }
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

# examine CNAs along chromosomes ------------------------------------------

if ( altType == "CNAs" ) {
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
    mychr <- "11"
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
    ## ## ## Make per chromosome plot
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

    ##     ggplot(allMutsChrD %>%
    ##            filter(chr==mychr, str_detect(alteration, "=-1")) %>%
    ##            mutate(angle=acos(fracOnFront) * 180 / pi) %>% 
    ##            mutate(isATM=str_detect(alteration, "ATM=-1"))
    ##            ) +
    ##         geom_point(mapping=aes(x=start/1e6, y=angle, color=isATM ##, color=isAln
    ##                                )) +
    ##         coord_cartesian(xlim=c(0,
    ##                                as.numeric(filter(chrLens, chr==mychr) %>%
    ##                                           select(chrLen))/1e6)) +
    ##         ## geom_label_repel(data=filter(allMutsChrD,
    ##         ##                              chr==mychr & isDrivingAlteration) %>%
    ##         ##                      mutate(angle=acos(fracOnFront) * 180 / pi),
    ##         ##                  mapping=aes(x=start/1e6, y=angle,
    ##         ##                              label=alteration## , color=isAln
    ##         ##                              ),
    ##         ##                  alpha=myalpha, col="red", show.legend=F) +
    ##         labs(x=xlab, y="angle to Pareto front") +
    ##         theme_gray()
    ##     ggsave("ATM1_CNA.pdf", height=3, width=5)
        
    ##     plot_grid(p1, p2, p3, p4, p5, labels=toupper(letters[1:5]),
    ##               nrow=5, ncol=1, align="v")
    ##     ggsave(sprintf("chrPos_frequency_fracOnFront_log10qNorm_%s_%s_%s.pdf",
    ##                    altType, s, mychr),
    ##            height=3*5, width=12*1)
    ## }
    ## ## ## End of per chromosome plot
    
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
}

# Characterize archetypes ---------------------------------------

MSigDB <- read_csv("MSigDBenrichment_continuous_All.csv")

toShow <- as_tibble(unlist(
    list(A2=c("REACTOME RESPIRATORY ELECTRON TRANSPORT"),
         A3=c("KEGG DNA REPLICATION", "MITOTIC CELL CYCLE"),
         A4=c("KEGG ALLOGRAFT REJECTION"))
)) %>% rename(`Feature Name`=value)

ggplot(inner_join(MSigDB, toShow) %>% 
    select(`archetype #`, `Feature Name`, `Mean Difference`) %>%
    complete(`archetype #`, `Feature Name`)  %>%
    replace_na(list(`Mean Difference`=0))) +
geom_line(aes(x=`archetype #`, 
              y=`Mean Difference`, col=`Feature Name`,
              group=`Feature Name`), size=1.5) +
    theme_grey()

source("../ParTI/TCGAscripts/hallmarkOmeter.inc.R")
toShow <- list_to_tibble(toShow, key="Group", value="Feature Name")

MSigDBord <- read_rds(ifelse(isMetabric,"../TCGA/MSigDBord.rds","../MSigDBord.rds"))
archProps <-
    inner_join(MSigDB, toShow) %>% 
    select(`archetype #`, `Feature Name`, `Mean Difference`, Group) %>%
mutate(`Feature Name`=fmtMSigDBnames(`Feature Name`))##  %>%
##     mutate(`Feature Name`=factor(`Feature Name`, levels=MSigDBord)) %>%
##     ## complete(`archetype #`, `Feature Name`) %>%
## replace_na(list(`Mean Difference`=0))

ggplot(archProps %>%
       mutate(`Feature Name`=factor(`Feature Name`, levels=MSigDBord))) +
    geom_line(aes(x=`archetype #`, y=`Mean Difference`)) +
ylab(expression(paste(log[2], " fold change"))) +
facet_wrap(Group ~ `Feature Name` , scales="free_y") +
    geom_hline(aes(yintercept=0), linetype=2) +
    theme_grey()
ggsave("hallmarkOmeter.pdf", height=10, width=14);

ggplot(archProps %>%
       mutate(`Feature Name`=factor(`Feature Name`, levels=MSigDBord))) +
    geom_col(aes(x=`Feature Name`, y=`Mean Difference`)) +
    labs(y="mean diff. in log expression", x="MSigDB gene group")+ 
    coord_flip() +
    ## geom_errorbarh(aes(y=MSigDB,
    ##                    xmin=expression - SE,
    ##                    xmax=expression + SE)) +
    ## facet_grid(. ~ `archetype #`, scales="free_x") +
    facet_grid(. ~ `archetype #`) +
    theme_grey()
ggsave("hallmarkOmeter_byArch.pdf", height=6, width=14)

# 3D plots of samples and mutations ---------------------------------------

## Compute scalar product of muts with vectors pointing to
## different archetypes, to later select best aligned mutation
## with each archetype
i <- 1
load(sprintf("Xproj_%s.rda", altType))
if ( isMetabric & altType == "SNVs") {
    rownames(Xproj) <- paste(rownames(Xproj), "=1", sep="")
}
Xproj[1:5,]
archProj <- t(sapply(1:nrow(arcsOrig), function(i) {
    m0 <- arcsOrig[i,] - geneExprAvg
                                        #- healthyProfile - avgE
    matrix(m0, nrow=1) %*% projMat
}))

## plot3d(t(apply(E, 2, function(x) { x + healthyProfile - geneExprAvg })) %*% projMat)
## spheres3d(archProj, radius=5)

## hPt <-
##     (- matrix(apply(E, 1, mean), nrow=1) ) %*% projMat 

s <- ifelse(altType=="CNAs", 10, 3);

## Import archetype color scheme
source("../ParTI/TCGAscripts/hallmarkOmeter.inc.R")

superArcs <-
    read_csv("~/work/cancerTaskAtlas/TCGA/ALL_UCSC/arcsOrig_genes.csv",
             col_names=unlist(read_csv("~/work/cancerTaskAtlas/TCGA/ALL_UCSC/geneNamesAfterExprFiltering.list",
                                       col_names=F)))

commonGenes <- intersect(colnames(arcsOrig), colnames(superArcs))

## Compute distances to super archetypes, pick smallest
matchSAtab <-
    sapply(1:nrow(arcsOrig), function(i) {
        apply(superArcs[,commonGenes], 1, function(x) {
            sqrt(sum(( (arcsOrig[i,commonGenes] - geneExprAvg[commonGenes]) - x )^2))
        })
    })
## Columns are archetypes in current cancer type; rows are
## super-archetypes

## apply(matchSAtab, 2, function(x) { arcCols[x^2 / sum(x^2) < .1] })
## tCols <- apply(matchSAtab, 2, function(x) { arcCols[which.min(x^2 / sum(x^2))] })

##################################################
## 1. Align by selected MSigDBs

arcsMSig <- read_csv("MSigDBenrichment_continuous_significant.csv")
## arcsMSig <- arcsMSig %>% filter(`P value (Mann-Whitney)` < 1e-3)
arcsMSig <- arcsMSig %>% filter(`Mean Difference` > 0.1)
if ( isMetabric ) {
    SAMSig <-
        read_csv("../TCGA/ALL_UCSC/MSigDBenrichment_continuous_significant.csv")
} else {
    SAMSig <-
        read_csv("../ALL_UCSC/MSigDBenrichment_continuous_significant.csv")
}

## SAMSig <- SAMSig %>% filter(`P value (Mann-Whitney)` < 1e-3)
SAMSig <- SAMSig %>% filter(`Mean Difference` > 0.1)
featsUniv <- union(SAMSig %>% select(`Feature Name`),
                   arcsMSig %>% select(`Feature Name`))
## featsUniv <- SAMSig %>% select(`Feature Name`) %>% unique()
## read_rds("../MSigDBord.rds")
arcIdx <- 1;
## Match each archetype to the super archetype with highest overlap in
## enriched MSigDBs
SAmapping <-
    map_int(unlist(arcsMSig %>% select(`archetype #`) %>% unique), function(arcIdx) {
        SAidx <- 1;
        SAidx <- 2;
        arcScores <-
            map_dbl(unlist(SAMSig %>% select(`archetype #`) %>% unique),
                    function(SAidx) {
                        SAfeats <- SAMSig %>% filter(`archetype #` == SAidx) %>% select(`Feature Name`)
                        arcFeats <- arcsMSig %>%
                            filter(`archetype #` == arcIdx) %>%
                            select(`Feature Name`)
                        ## nrow(SAfeats)
                        ## nrow(arcFeats)

                        ## myUniv <- union(SAfeats, arcFeats)
                        ## expIntersect <- nrow(SAfeats) * nrow(arcFeats) / nrow(myUniv)
                        expIntersect <- nrow(SAfeats) * nrow(arcFeats) / nrow(featsUniv)

                        ## phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
                        ## q,x: number of white balls drawn without
                        ## replacement from an urn which contains both black and white balls.
                        ## m: the number of white balls in the urn.
                        ## n: the number of black balls in the urn.
                        ## k: the number of balls drawn from the urn.
                        p <- phyper(q=intersect(SAfeats, arcFeats) %>% nrow,
                                    m=SAfeats %>% nrow,
                                    n=nrow(featsUniv) - nrow(SAfeats),
                                    k=arcFeats %>% nrow,
                                    lower.tail=F)
                        
                        foldEnrich <- nrow(intersect(SAfeats, arcFeats)) / expIntersect;
                        ## if ( p > .01/(nArchs * 5) ) { foldEnrich <- 0 }
                        ## return(foldEnrich)
                        return(p)
                    })
        ## Make sure the overlap with the closest super-archetype is significant
        cutOff <- .01 / (nArchs * 5)
        if ( min(arcScores) > cutOff ) {
            return(NA)
        } else {
            return(which.min(arcScores))
        }
    })

##################################################

## ## Match archetypes to super archetypes (minimal distance)

## if ( isMetabric ) {
##     SAMSigAll <-
##         read_csv("../TCGA/ALL_UCSC/MSigDBenrichment_continuous_All.csv")
## } else {
##     SAMSigAll <-
##         read_csv("../ALL_UCSC/MSigDBenrichment_continuous_All.csv")
## }

## MSigDBAll <- read_csv("MSigDBenrichment_continuous_All.csv")

## Amat <- inner_join(
##     tibble("Feature Name"=as.character(unlist(toShow))),
##     MSigDBAll) %>% select(`Feature Name`, `archetype #`, `Mean Difference`) %>%
##     spread(`archetype #`, `Mean Difference`)
## SAmat <- inner_join(
##     tibble("Feature Name"=as.character(unlist(toShow))),
##     SAMSigAll) %>% select(`Feature Name`, `archetype #`, `Mean Difference`) %>%
##     spread(`archetype #`, `Mean Difference`)

## ## make sure MSigDB are ordered the same way
## sum(Amat[,1] != SAmat[,1]) == 0 

## library(gtools)
## allPerms <- permutations(n=ncol(SAmat) - 1, r=ncol(Amat) - 1);
## allDists <-
##     apply(allPerms, 1, function(myPerm) {
##         sqrt(sum(( SAmat[,-1][,myPerm] - Amat[,-1] )^2))
##     })
## plot(sort(allDists))
## bestPerms <- order(allDists)[1:5]
## apply(allPerms[bestPerms,], 1, function(x) {
##     names(arcCols)[x]
## })

## ## SAmapping <- allPerms[which.min(allDists),]

##################################################

tCols <- arcCols[SAmapping]
if ( isMetabric ) {
    tCols[4] <- "grey";
    names(tCols)[4] <- "HER2";
}
names(tCols) %>% enframe %>% write_csv("archetypes_tasks.csv")

##################################################
## Represent the polyhedron in 3D

## 2 / diff(range(archProj)) is about 0.007
## 15 / diff(range(archProj)) is about 0.05 on BLCA

posInPCspaceV <- posInPCspace; ## for visualization
archProjV <- archProj;
if ( nArchs <= 3 ) {
    posInPCspaceV[,3] <- 0;
    archProjV[,3] <- 0;
}

open3d()
pp <- dget("characterizeArchetypes_3DplotSetup.R")
par3d(pp)
plot3d(posInPCspaceV[,1:3],
       ## alpha=.35, col="grey",
       col="black", type="s", radius=diff(range(archProjV)) * 0.007, shininess=100,
       ## alpha=.5,
       box=F, axes=F,
       ## aspect=F,
       ## xlab="PC1", ylab="PC2", zlab="PC3", axes=T, box=F,
       xlab="", ylab="", zlab="",
       xlim=range(c(archProj[,1:3], posInPCspace[,1:3])),
       ylim=range(c(archProj[,1:3], posInPCspace[,1:3])),
       zlim=range(c(archProj[,1:3], posInPCspace[,1:3])))
       ## xlim=range(c(archProjV[,1], posInPCspace[,1])),
       ## ylim=range(c(archProjV[,2], posInPCspace[,2])),
       ## zlim=range(c(archProjV[,3], posInPCspace[,3])))
## spheres3d(hPt, col="green", radius=5)
## spheres3d(archProjV, radius=15, col=tCols)


## text3d(archProjV[,1], archProjV[,2], archProjV[,3],
##        ## allPerms[[topPerm]],
##        1:nrow(archProjV),
##        adj=2.5)
sapply(1:(nrow(archProjV)-1), function(i) {
    sapply(seq(i+1, nrow(archProjV)), function(j) {
        segments3d(archProjV[c(i,j),1], archProjV[c(i,j),2],
                   archProjV[c(i,j),3], col="grey", lwd=10)
    })
})
## spheres3d(archProjV, radius=diff(range(archProjV)) * 0.05, col="red")
spheres3d(archProjV, radius=diff(range(archProjV)) * 0.05, col=tCols)
## ## Save the scene:
## pp <- par3d(no.readonly=T)
## pp$userMatrix[1:3,1:3] <-
##     diag(c(1,1,1))
##     matrix(c(0, 1, 0,
##              -1, 0, 0,
##              0, 0, 1), ncol=3, byrow=T)
##     matrix(c(0, -1, 0,
##              -1, 0, 0,
##              0, 0, -1), ncol=3, byrow=T)
## dput(pp, file="characterizeArchetypes_3DplotSetup.R")
pp <- dget("characterizeArchetypes_3DplotSetup.R")
par3d(pp)
rgl.snapshot("characterizeArchetypes_3Dplot.png", fmt="png")
rgl.close()

## Show randomized version of the same plot

## rndE <- E;
## for (i in 1:ncol(rndE)) {
##     rndE[,i] <- sample(rndE[,i])
## }
## ## for (i in 1:nrow(rndE)) {
## ##     rndE[i,] <- sample(rndE[i,])
## ## }
## ## dudi2rnd <- dudi.pca(t(rndE), scale=F, scannf=F, nf=nPCs)
## dudi2rnd <- dudi.pca(t(rndE[,-rev(order(apply(E, 2, sd)))[1:5]]),
##                      scale=F, scannf=F, nf=nPCs)

## rndPosInPCspace <- as.matrix(dudi2rnd$li);

rndPosInPCspace <- posInPCspaceV[,1:3]
for (i in 1:ncol(rndPosInPCspace)) {
    rndPosInPCspace[,i] <- sample(rndPosInPCspace[,i])
}

## plot(rndPosInPCspace[,1], rndPosInPCspace[,2])
## plot(rndPosInPCspace[,1], rndPosInPCspace[,3])
## which(abs(rndPosInPCspace[,1]) > 50)
## which(abs(rndPosInPCspace[,2]) > 50)
## which(abs(rndPosInPCspace[,3]) > 50)
## plot(apply(E, 2, sd))
## plot(apply(E, 1, var))
## rev(order(apply(E, 2, sd)))[1:10]
## hist(apply(E, 2, var), 20)
## hist(apply(E, 1, var), 20)

## plot(density(E[,1]), type="n", ylim=c(0,1))
## for ( i in sample(1:ncol(E), 5) ) {
##     lines(density(E[,i], bw=.2))
## }
## for ( i in rev(order(apply(E, 2, sd)))[1:5] ) {
##     lines(density(E[,i], bw=.2), col="red")
## }
## lines(density(E[,270]), col="red")
## lines(density(E[,5]), col="red")
## plot(apply(E, 2, mean), 
##      apply(E, 2, sd))
## which(apply(E, 2, sd) > 2)

## plot3d(rnorm(1e3), rnorm(1e3), rnorm(1e3)/10, aspect=T,
##        xlim=c(-4,4), ylim=c(-4, 4), zlim=c(-4,4))
## rgl.close()

## plot3d(rndPosInPCspace[,1:3],
##        ## alpha=.35, col="grey",
##        col="black", type="s", radius=diff(range(archProjV)) * 0.007, shininess=100,
##        ## alpha=.5,
##        ## aspect=F,
##        ## xlab="PC1", ylab="PC2", zlab="PC3", axes=T, box=F,
##        xlab="", ylab="", zlab="", box=F, axes=F,
##        xlim=range(c(archProjV[,1:3], posInPCspaceV[,1:3])),
##        ylim=range(c(archProjV[,1:3], posInPCspaceV[,1:3])),
##        zlim=range(c(archProjV[,1:3], posInPCspaceV[,1:3])))
## ## sapply(1:(nrow(archProjV)-1), function(i) {
## ##     sapply(seq(i+1, nrow(archProjV)), function(j) {
## ##         segments3d(archProjV[c(i,j),1], archProjV[c(i,j),2],
## ##                    archProjV[c(i,j),3], col="grey", lwd=10)
## ##     })
## ## })
## ## spheres3d(archProjV, radius=diff(range(archProjV)) * 0.05, col="red")
## ## pp <- dget("characterizeArchetypes_3DplotSetup.R")
## par3d(pp)
## ## rgl.snapshot("characterizeArchetypes_3Dplot_random.png", fmt="png")
## rgl.close()

if ( length(grep("LGG_UCSC", getwd())) == 1 ) {
    dType <- "original"
    for ( dType in c("original", "shuffled")) {
        if ( dType == "original" ) {
            myData <- posInPCspaceV[,1:2];
            myArch <- archProjV[,1:2];
        } else {
            myData <- rndPosInPCspace[,1:2];
            myArch <- archetypes::archetypes(myData, 3)$archetypes
        }

        pdf(sprintf("archFit_illust_%s.pdf", dType), height=5, width=5);
        plot(myData[,1], myData[,2],
             xlim=range(c(myData, myArch)),
             ylim=range(c(myData, myArch)),
             pch=20)
        ## lines(archProjV[c(1:nrow(archProjV), 1),1],
        ##       archProjV[c(1:nrow(archProjV), 1),2], lwd=2)
        lines(myArch[c(1:nrow(myArch), 1),1],
              myArch[c(1:nrow(myArch), 1),2], lwd=2)
        myChull <- grDevices::chull(myData)
        myChull <- c(myChull, myChull[1])
        lines(myData[myChull,], col="grey")
        dev.off()
    }
}


##################################################

contClin <- read.table("continuousClinicalData_reOrdered.tsv",
                       as.is=T, sep="\t", h=T)
if ( length(grep("ESTIMATE", colnames(contClin))) > 0 ) {
    ggplot(tibble(Axis1=posInPCspaceV[,1], Axis2=posInPCspaceV[,2],
                  Purity=contClin[isPT,"ESTIMATE"])) +
        geom_point(aes(x=Axis1, y=Axis2, color=Purity)) +
        geom_point(aes(x=x, y=y), color="red", size=5,
                   data=tibble(x=archProjV[,1], y=archProjV[,2]))
    ggsave("purity_PC1_PC2.pdf", height=4, width=5)
}

## ggplot(tibble(Axis1=posInPCspaceV[,1], Axis2=posInPCspaceV[,2],
##               nSNVs=contClin[isPT,"number_of_SNVs_frac"])) +
##     geom_point(aes(x=Axis1, y=Axis2, color=nSNVs)) +
##     geom_point(aes(x=x, y=y), color="red", size=5,
##                data=tibble(x=archProjV[,1], y=archProjV[,2]))
## ggsave("SNVsFrac_PC1_PC2.pdf", height=4, width=5)

## ggplot(tibble(Axis1=posInPCspaceV[,1], Axis2=posInPCspaceV[,2],
##               nCNAs=contClin[isPT,"number_of_CNAs_frac"])) +
##     geom_point(aes(x=Axis1, y=Axis2, color=nCNAs)) +
##     geom_point(aes(x=x, y=y), color="red", size=5,
##                data=tibble(x=archProjV[,1], y=archProjV[,2]))
## ggsave("CNAsFrac_PC1_PC2.pdf", height=4, width=5)

plot3d(posInPCspaceV[,1:3],
       col=heat.colors(100, alpha = 1)[round(contClin[isPT,"number_of_SNVs_frac"] * 100)],
       type="s", radius=diff(range(archProjV)) * 0.007, shininess=100,
       box=F, axes=F,
       xlab="", ylab="", zlab="",
       xlim=range(c(archProjV[,1:3], posInPCspaceV[,1:3])),
       ylim=range(c(archProjV[,1:3], posInPCspaceV[,1:3])),
       zlim=range(c(archProjV[,1:3], posInPCspaceV[,1:3])))
sapply(1:(nrow(archProjV)-1), function(i) {
    sapply(seq(i+1, nrow(archProjV)), function(j) {
        segments3d(archProjV[c(i,j),1], archProjV[c(i,j),2],
                   archProjV[c(i,j),3], col="grey", lwd=10)
    })
})
spheres3d(archProjV, radius=diff(range(archProjV)) * 0.05, col="red")
rgl.close()

plot3d(posInPCspaceV[,1:3],
       col=heat.colors(100, alpha = 1)[round(contClin[isPT,"number_of_CNAs_frac"] * 100)],
       type="s", radius=diff(range(archProjV)) * 0.007, shininess=100,
       box=F, axes=F,
       xlab="", ylab="", zlab="",
       xlim=range(c(archProjV[,1:3], posInPCspaceV[,1:3])),
       ylim=range(c(archProjV[,1:3], posInPCspaceV[,1:3])),
       zlim=range(c(archProjV[,1:3], posInPCspaceV[,1:3])))
sapply(1:(nrow(archProjV)-1), function(i) {
    sapply(seq(i+1, nrow(archProjV)), function(j) {
        segments3d(archProjV[c(i,j),1], archProjV[c(i,j),2],
                   archProjV[c(i,j),3], col="grey", lwd=10)
    })
})
spheres3d(archProjV, radius=diff(range(archProjV)) * 0.05, col="red")
rgl.close()

## ##################################################
## ## 2. Align by distance in 5-arch PC space
## ## center LGG archetypes
## arcsOrig0 <- apply(arcsOrig, 1, function(x) { x - geneExprAvg } )
## ## load 5-arch PC
## pc5archs <-
##     read_csv("../ALL_UCSC/projOrig_varsXdims.csv", col_names=F) %>% as.data.frame
## dim(pc5archs)
## ## pc5archs[1:5,1:5]
## rownames(pc5archs) <- colnames(superArcs)
## pc5archs <- pc5archs[commonGenes,1:4];
## SAPC5 <- as.matrix(superArcs[,commonGenes]) %*% as.matrix(pc5archs)
## arcsPC5 <- t(arcsOrig0[commonGenes,]) %*% as.matrix(pc5archs)
## plot3d(0,0,0)
## spheres3d(SAPC5[,1:3], radius=10, col="red")
## spheres3d(arcsPC5[,1:3], radius=10, col="blue")
## text3d(SAPC5[,1:3], texts=1:5, adj=2.5, col="red")
## text3d(arcsPC5[,1:3], texts=1:4, adj=2.5, col="blue")
## ## Archetypes are rows
## dim(SAPC5)
## dim(arcsPC5)
## sapply(1:nrow(arcsPC5), function(arcsIdx) {
##     sapply(1:nrow(SAPC5), function(SAidx) {
##         ## sqrt(sum((SAPC5[SAidx,] - arcsPC5[arcsIdx,])^2))
##         ## cor(SAPC5[SAidx,], arcsPC5[arcsIdx,])
##         sum(SAPC5[SAidx,] * arcsPC5[arcsIdx,]) /
##             sqrt(sum(SAPC5[SAidx,]^2) * sum(arcsPC5[arcsIdx,]^2))
##     }) ## %>% rank
## })

##################################################

## saMemb <- read_rds("../saMemb_GMhallmarks.rds")

## ## Super-tasks according to clustering
## arcSAclust <-
##     tibble(archetype=sprintf("%s %d", cancerType, 1:nArchs)) %>%
##     inner_join(list_to_tibble(saMemb, key="Task", value="archetype"))

## tCols <- arcCols[unlist(arcSAclust %>% select(Task))]

## ## Re-label invasion into invasion or immune
## idx <- which(names(tCols) == "Invasion & signaling")[1]
## for ( idx in which(names(tCols) == "Invasion & signaling") ) {
##     newSAidx <- which.min(matchSAtab[,idx]);
##     cat(sprintf("Renaming '%s' into '%s'\n", names(tCols[idx]),
##                 names(arcCols[newSAidx])))
##     tCols[idx] <- arcCols[newSAidx]
##     names(tCols[idx]) <- names(arcCols[newSAidx])
## }

##################################################
## Collect direction vectors to archetypes, edges, faces

vRefs <- list()

## directions to archetypes:
vRefs[["archs"]] <- 
    sapply(1:nrow(archProj), function(i) {
        ## vRef <- archProj[i,] - archProj[healthyArch,]
        vRef <- archProj[i,]
        ## normalize archetypes so distance doesn't matter, only
        ## direction
        vRef <- vRef / sqrt(sum(vRef^2))        
        return(vRef)
    })
colnames(vRefs[["archs"]]) <- sprintf("archetype: %s", names(tCols))
## sum(vRefs[,2]^2)

## directions to edges
vRefs[["edges"]] <- matrix(NA, ncol(archProj), 0);
if ( nArchs >= 3 ) {
    for (i in 1:(nrow(archProj)-1)) {
        for (j in (i+1):nrow(archProj)) {
            vRef <- archProj[i,] + archProj[j,]
            vRef <- vRef / sqrt(sum(vRef^2))
            vRefs[["edges"]] <- cbind(vRefs[["edges"]], vRef)
            colnames(vRefs[["edges"]])[ncol(vRefs[["edges"]])] <-
                sprintf("edge: %s - %s", names(tCols)[i], names(tCols)[j])
        }
    }
}

## directions to faces
vRefs[["faces"]] <- matrix(NA, ncol(archProj), 0);
if ( nArchs >= 4 ) {
    for ( i in 1:(nrow(archProj)-2) ) {
        for ( j in (i+1):(nrow(archProj)-1) ) {
            for ( k in (j+1):nrow(archProj) ) {
                vRef <- archProj[i,] + archProj[j,] + archProj[k,]
                vRef <- vRef / sqrt(sum(vRef^2))
                vRefs[["faces"]] <- cbind(vRefs[["faces"]], vRef)
                colnames(vRefs[["faces"]])[ncol(vRefs[["faces"]])] <-
                    sprintf("face: %s - %s - %s",
                            names(tCols)[i], names(tCols)[j], names(tCols)[k])
            }
        }
    }
}

vRefs <- cbind(vRefs$archs, vRefs$edges, vRefs$faces)

## mutOrientToArchs <-
##     Xproj[unlist(allMuts %>% 
##                  select(alteration)),] %*%
##     vRefs
## rownames(Xproj) <- paste(rownames(Xproj), "=1", sep="")

## mutation vectors for aligned drivers
mutOrientToArchs <-
    ## as_tibble(Xproj) %>% mutate(alteration=rownames(Xproj)) %>%
    ## filter(alteration %in%
    ##        (allMuts %>% filter(isAln & isDrivingAlteration) %>% 
    ##         select(alteration))
    ##        )
    as_tibble(Xproj) %>% mutate(alteration=rownames(Xproj)) %>%
    inner_join(
        allMuts %>% filter(isAln & isDrivingAlteration) %>% 
        select(alteration)
    )
tmp <- mutOrientToArchs %>% select(-alteration) %>% as.matrix
rownames(tmp) <-
    mutOrientToArchs %>% select(alteration) %>% unlist %>% as.character
mutOrientToArchs <- tmp
rm(tmp)
## ## Same, but if there is only one driver, the matrix collapses into
## ## a vector
## mutOrientToArchs <-
##     Xproj[allMuts %>% filter(isAln & isDrivingAlteration) %>% 
##           select(alteration) %>% unlist %>% as.character,]

## set irrelevant dimensions to 0
mutOrientToArchs[,setdiff(1:ncol(mutOrientToArchs), 1:(nArchs-1))] <- 0
## normalize directions
mutOrientToArchs <- 
    apply(mutOrientToArchs, 1, function(x) { x / sqrt(sum(x^2))})
## compute scalar products
180 * acos(sort(
          t(mutOrientToArchs) %*% vRefs
          )) / pi
mutOrientToArchs <- t(mutOrientToArchs) %*% vRefs
## rownames(mutOrientToArchs) <-
##     unlist(allMuts %>% filter(isAln & isDrivingAlteration) %>% 
##                  select(alteration))

mutOrientToArchs <-
    t(apply(mutOrientToArchs, 1, function(x) {
        x[x < max(x)] <- 0
        return(x)
    }))
tibble(
    cancer=cancerType,
    alteration=rownames(mutOrientToArchs),
    orientation=colnames(mutOrientToArchs)[
        apply(mutOrientToArchs, 1, which.max)]) %>%
    inner_join(mutsAlongFront$Front %>% select(alteration, frequency)) %>% 
    write_csv(sprintf("%s_aln_driver_orient.csv", altType))

sapply(1:nrow(archProj), function(i) {
    list(rev(sort(mutOrientToArchs[,i]))[1:5])
})

mutAngleToArchs <-
    ## ( t(apply(Xproj[unlist(allMuts %>% filter(isAln & isDrivingAlteration) %>% 
    ( t(apply(Xproj[unlist(allMuts %>% filter(isDrivingAlteration) %>% 
                           select(alteration)),], 1, function(x) { x / sqrt(sum(x^2)) })) %*%
      vRefs %>% acos() ) * 360 / (2 * pi)
## We compute the angle in the plane spanned by the mutation and
## archetype vector; a plane has three archetypes which puts these 120
## degrees apart; thus, a mutation can never be off an archetype by
## more than 60 degrees. 
ggplot(enframe(apply(mutAngleToArchs, 1, min))) +
    geom_histogram(aes(x=value), breaks=seq(0, 60, by=5)) +
    coord_cartesian(xlim=c(0, 60)) +
    labs(x="angle to archetype [degrees]", y="number of driver mutations") +
    geom_vline(aes(xintercept=60/2)) + theme_grey()
ggsave(sprintf("angleToArchetypes_%s.pdf", altType), height=4, width=4)
sum(apply(mutAngleToArchs, 1, min) < 30)
sum(apply(mutAngleToArchs, 1, min) > 30)
sum(apply(mutAngleToArchs, 1, min) < 30) / nrow(mutAngleToArchs)

m <- "GATA3=1"
m <- "BRAF=1"
m <- "TP53=1"
m <- "SMAD4=-1";
m <- "ERBB2=1";
m <- "BRCA1=1";
m <- "CCND3=2"
m <- "CCNE1=2"
m <- "CTNNB1=1"
m <- "CDH1=1"
m <- "MYC=1"
m <- "EGFR=1"
i <- 1
i <- 2
i <- 3
i <- 4

## save.image(sprintf("alt3Dplots_%s.rda", altType))
## load(sprintf("alt3Dplots_%s.rda", altType))

sapply(1:nrow(archProj), function(i) {
    m <- names(rev(sort(mutOrientToArchs[,i]))[1])
    if ( mutOrientToArchs[m,i] == 0 ) { return(); }
    ## x <- Xproj[m,1:3]
    ## v <- x - hPt;
    v <- Xproj[m,1:3]

    ptCol <- rep("grey", ncol(Mp))
    ## if ( isMetabric ) {
    ##     ptCol[Mp[gsub("=1", "", m),]==1] <- "red"
    ## } else {
    ##     ptCol[Mp[m,]==1] <- "red"
    ## }
    ptCol[Mp[m,]==1] <- "red"

    cat(sprintf("%s - arch %d\n", m, i))
    open3d()
    pp <- dget("mutVec_plotSetup.R")
    par3d(pp)

    greyAlpha <- ifelse(nArchs<=3, .5, .1);
    plot3d(posInPCspaceV[,1:3],
           col=ptCol, type="s", shininess=100,
           radius=diff(range(archProjV)) * 0.007,
           alpha=(ptCol == "red") * (1-greyAlpha) + greyAlpha,
           ## alpha=.5,
           box=F, axes=F, aspect=F, xlab="", ylab="", zlab="",
           ## box=T, axes=T, aspect=T, xlab="PC1", ylab="PC2", zlab="PC3",
           xlim=range(c(archProj[,1], posInPCspace[,1])),
           ylim=range(c(archProj[,2], posInPCspace[,2])),
           zlim=range(c(archProj[,3], posInPCspace[,3])))
    ## spheres3d(hPt, col="green", radius=5)
    spheres3d(archProjV, radius=diff(range(archProjV)) * 0.05, col=tCols)
    ## text3d(archProjV[,1], archProjV[,2], archProjV[,3],
    ##        allPerms[[topPerm]],
    ##        ## 1:nrow(archProjV),
    ##        adj=2.5)
    sapply(1:(nrow(archProjV)-1), function(i) {
        sapply(seq(i+1, nrow(archProjV)), function(j) {
            segments3d(archProjV[c(i,j),1], archProjV[c(i,j),2],
                       archProjV[c(i,j),3], col="grey", lwd=10)
        })
    })
    cmGrey <- apply(posInPCspaceV[ptCol=="grey",], 2, mean)
    if ( nArchs == 3 ) {
        v[3] <- 0;
        cmGrey[3] <- 0;
    }
    
    v <- (v / sqrt(sum(v^2))) * (diff(range(archProjV)) / 4)

    ## ## Save the scene:
    ## pp <- par3d(no.readonly=T)
    ## dput(pp, file="mutVec_plotSetup.R")
    pp <- dget("mutVec_plotSetup.R")
    par3d(pp)

    ## Compute the vector that points to the camera
    ## myAx <- solve(par3d(no.readonly=T)$userMatrix[1:3,1:3])[,1] * 100
    ## lines3d(c(0, myAx[1]), c(0, myAx[2]), c(0, myAx[3]), col="green")
    ## myAx <- solve(par3d(no.readonly=T)$userMatrix[1:3,1:3])[,2] * 100
    ## lines3d(c(0, myAx[1]), c(0, myAx[2]), c(0, myAx[3]), col="blue")
    myAx <- solve(par3d(no.readonly=T)$userMatrix[1:3,1:3])[,3] * 100
    ## lines3d(c(0, myAx[1]), c(0, myAx[2]), c(0, myAx[3]), col="red")
    
    arrow3d(cmGrey ## + diff(range(archProjV)) * myAx / 400
           ,
            cmGrey + v ## + diff(range(archProjV)) * myAx / 400
           ,
            width=diff(range(archProjV)) * 0.0007, shininess=100,
            color="black", type="rotation")
    ## arrow3d(c(0, 0, 0) + c(-50,100,-50),
    ##         v*s + c(-50,100,-50),
    ##         width=diff(range(archProjV)) * 0.0007)
    rgl.snapshot(sprintf("mutVec_%s.png", m), fmt="png")
    rgl.close()

    pdf(sprintf("mutVec_hist_%s.pdf", m), height=2, width=3, pointsize=16)
    par(mar=c(4, 4, 0, 0))
    h <- hist(90*(1-rratios), 20, freq=F, col="grey", main="", xlim=c(0,90),
              xlab="Angle to front\n(degrees)", ylab="", lab=c(3, 3, 3))
    mRatio <- mutsAlongFront[["Front"]] %>% filter(alteration == m) %>%
        select(fracOnFront) %>% as.numeric
    abline(v=90*(1-mRatio))
    dev.off();
    
    ## text3d(s * v[1], s * v[2], s * v[3],
    ##        sub("_CNA$", "", sub("=1$", "", m)),
    ##        adj=c(.5, 0))
})

if ( cancerType == "COAD" ) {
    ## Stain MSI-H samples
    discClin <- read.table("discreteClinicalData_reOrdered.tsv",
                           as.is=T, sep="\t", h=T)

    unique(discClin[,"CDE_ID_3226963"])
    ptCol <- (discClin[isPT,"CDE_ID_3226963"] == "MSI-H") + 1
    ptCol <- c("grey", "red")[ptCol]
    ## ptCol <- (discClin[,"CDE_ID_3226963"] == "MSS") + 1
    ## ptCol <- (discClin[,"CDE_ID_3226963"] == "MSI-L") + 1
    ## ptCol <- as.numeric(as.factor(discClin[,"CDE_ID_3226963"])) - 2

    open3d();
    plot3d(posInPCspaceV,
           col=ptCol, type="s", shininess=100,
           radius=diff(range(archProjV)) * 0.007,
           ## alpha=(ptCol == "red") * .9 + .1,
           ## alpha=1,
           box=F, axes=F, aspect=F, xlab="", ylab="", zlab="",
           ## box=T, axes=T, aspect=T, xlab="PC1", ylab="PC2", zlab="PC3",
           xlim=range(c(archProjV[,1], posInPCspaceV[,1])),
           ylim=range(c(archProjV[,2], posInPCspaceV[,2])),
           zlim=range(c(archProjV[,3], posInPCspaceV[,3])))
    spheres3d(archProjV, radius=diff(range(archProjV)) * 0.05, col=tCols)
    sapply(1:(nrow(archProjV)-1), function(i) {
        sapply(seq(i+1, nrow(archProjV)), function(j) {
            segments3d(archProjV[c(i,j),1], archProjV[c(i,j),2],
                       archProjV[c(i,j),3], col="grey", lwd=10)
        })
    })
    
    pp <- dget("mutVec_plotSetup.R")
    par3d(pp)
    rgl.snapshot("COAD_MSI-H.png", fmt="png")
    rgl.close()
}
