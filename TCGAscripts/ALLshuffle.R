rm(list=ls())

## What cancer types actually have arch. 1?

library(ade4)
library(rgl)
library(ggrepel)
library(stringr)
library(tidyverse)

## Import archetype color scheme
## source("../ParTI/TCGAscripts/hallmarkOmeter.inc.R")

samples <- read_csv("expMatrix.csv", col_names=F);
clinical <- read_delim("discreteClinicalData_reOrdered.tsv", delim="\t")
table(clinical$sample_type)
dim(samples)
samplesPT <- samples[clinical$sample_type == "Primary Tumor",]
samplesPT <- as.matrix(as.data.frame(samplesPT))

sum(is.na(samplesPT)) ## 7000 NAs
dim(samplesPT)
which(apply(is.na(samplesPT), 1, sum)>0) ## samples with NAs
## samplesPT <- samplesPT[,-which(apply(is.na(samplesPT), 2, sum)>0)]
which(apply(is.na(samplesPT), 2, sum)>0) ## genes with NAs
samplesPT <- samplesPT[,-which(apply(is.na(samplesPT), 2, sum)>0)]

arcs <- read_csv("arcsOrig_genes.csv", col_names=F);
dim(arcs)
nArchs <- nrow(arcs);
nPCs <- nArchs - 1;
if ( nPCs == 2 ) { nPCs == 3 }

save.image("ALLshuffle_prep.rda")

samplesPTr <- samplesPT;
for (i in 1:ncol(samplesPTr)) {
    samplesPTr[,i] <- sample(samplesPTr[,i])
}

samplesPT[1,1:5]
samplesPTr[1,1:5]
dim(samplesPT); dim(samplesPTr);
dudi1 <- dudi.pca(samplesPT, center=F, scale=F, scannf=F, nf=nPCs)
dudi1r <- dudi.pca(samplesPTr, center=F, scale=F, scannf=F, nf=nPCs)

## mutations <- read_csv("mutMatrix_reOrdered_booleanized_justData.csv",
##                       col_names=F, na = c("", "NA", "NaN"))
## mutNames <- read_csv("mutMatrix_reOrdered_booleanized_geneNames.list",
##                      col_names=F)
## colnames(mutations) <- sub("=1$", "", as.data.frame(mutNames)[,1])

## samples <- read_delim("projOrig_varsXdims.tsv", "   ", col_names=F)
## samples <- read_delim("tmp.tsv", "\t", col_names=F)

samples <- dudi1$li
samples <- dudi1r$li ## This randomization has no structure

## for (i in 1:ncol(samples)) {
##     samples[,i] <- sample(samples[,i])
## }

sRadius <- diff(range(samples)) * 2 / 325
## par3d(read_rds("characterizeArchetypes_3Dplot_setting.rds"))
open3d()
plot3d(samples[,1:3], ## alpha=1, 
       col="black", type="s", radius=sRadius,
       xlab="PC1", ylab="PC2", zlab="PC3", axes=T,
       ## xlab="", ylab="", zlab="", axes=F,
       box=F, shininess=100)
## xlim=range(arcs$X1), ylim=range(arcs$X2), zlim=range(arcs$X3))
## xlim=range(samples[,1:3]))
spheres3d(rep(min(samples[,1]), nrow(samples)),
          samples[,2],
          samples[,3],
          radius=sRadius, col="grey", shininess=100)
spheres3d(samples[,1],
          rep(max(samples[,2]), nrow(samples)),
          samples[,3],
          radius=sRadius, col="grey", shininess=100)
spheres3d(samples[,1],
          samples[,2],
          rep(min(samples[,3]), nrow(samples)),
          radius=sRadius, col="grey", shininess=100)

## spheres3d(arcs[,1:3], col=arcCols, radius=15)
## i <- 1;
## for (i in 1:(nrow(arcs)-1)) {    
##     j <- 2;
##     for (j in 2:nrow(arcs)) {
##         lines3d(unlist(arcs[c(i,j),1]),
##                 unlist(arcs[c(i,j),2]),
##                 unlist(arcs[c(i,j),3]), col="grey")
##     }
## }

## for(i in 1:nrow(arcs)) {
##     text3d(arcs[i,1], arcs[i,2], arcs[i,3], i, adj=3)
## }

## pp <- par3d(no.readonly=TRUE)
## write_rds(par3d()$userMatrix, "characterizeArchetypes_3Dplot_setting.rds")

## rgl.postscript("characterizeArchetypes_3Dplot.svg", fmt="svg");
rgl.viewpoint(userMatrix=read_rds("characterizeArchetypes_3Dplot_setting.rds"),
              fov=30, zoom=0.8)

x11(); plot(as.data.frame(samples[,1:3]))
