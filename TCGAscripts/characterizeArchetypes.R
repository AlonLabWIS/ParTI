rm(list=ls())

## What cancer types actually have arch. 1?

library(fdrtool)
library(ade4)
library(rgl)
library(ggrepel)
library(stringr)
library(tidyverse)

## Import archetype color scheme
source("../ParTI/TCGAscripts/hallmarkOmeter.inc.R")

arcs <- read_csv("arcs_dims.csv", col_names=F)
samples <- read_csv("pcsOrig_samplesXdims.csv", col_names=F)
clinical <- read_delim("discreteClinicalData_reOrdered.tsv", delim="\t")

mutations <- read_csv("mutMatrix_reOrdered_booleanized_justData.csv",
                      col_names=F, na = c("", "NA", "NaN"))
mutNames <- read_csv("mutMatrix_reOrdered_booleanized_geneNames.list",
                     col_names=F)
colnames(mutations) <- sub("=1$", "", as.data.frame(mutNames)[,1])

## samples <- read_delim("projOrig_varsXdims.tsv", "   ", col_names=F)
## samples <- read_delim("tmp.tsv", "\t", col_names=F)

nArchs <- nrow(arcs);
nPCs <- nArchs - 1;
if ( nPCs == 2 ) { nPCs == 3 }

samplesS <- samples;
for (i in 1:nPCs) {
    samplesS[,i] <- sample(samplesS[,i])
}
## par3d(read_rds("characterizeArchetypes_3Dplot_setting.rds"))
plot3d(samples[,1:3], ## alpha=1, 
       col="black",
       ## type="s", radius=2, shininess=100, ## sphere version
       xlab="PC1", ylab="PC2", zlab="PC3", axes=T,
       ## xlab="", ylab="", zlab="", axes=F,
       box=F, aspect=F,
       xlim=range(arcs$X1), ylim=range(arcs$X2), zlim=range(arcs$X3))
## xlim=range(samples[,1:3]))
points3d(rep(min(arcs[,1]), nrow(samples)),
         as.data.frame(samples)[,2],
         as.data.frame(samples)[,3],
         size=2, col="grey")
points3d(as.data.frame(samples)[,1],
         rep(max(arcs[,2]), nrow(samples)),
         as.data.frame(samples)[,3],
         size=2, col="grey")
points3d(as.data.frame(samples)[,1],
         as.data.frame(samples)[,2],
         rep(min(arcs[,3])-50, nrow(samples)),
         size=2, col="grey")
## spheres3d(rep(min(arcs[,1]), nrow(samples)),
##           as.data.frame(samples)[,2],
##           as.data.frame(samples)[,3],
##           radius=2, col="grey", shininess=100)
## spheres3d(as.data.frame(samples)[,1],
##           rep(max(arcs[,2]), nrow(samples)),
##           as.data.frame(samples)[,3],
##           radius=2, col="grey", shininess=100)
## spheres3d(as.data.frame(samples)[,1],
##           as.data.frame(samples)[,2],
##           rep(min(arcs[,3])-50, nrow(samples)),
##           radius=2, col="grey", shininess=100)

points3d(arcs[,1:3], col=arcCols, size=15)
## spheres3d(arcs[,1:3], col=arcCols, radius=15)
i <- 1;
for (i in 1:(nrow(arcs)-1)) {    
    j <- 2;
    for (j in 2:nrow(arcs)) {
        lines3d(unlist(arcs[c(i,j),1]),
                unlist(arcs[c(i,j),2]),
                unlist(arcs[c(i,j),3]), col="#808080")
    }
}
## for(i in 1:nrow(arcs)) {
##     text3d(arcs[i,1], arcs[i,2], arcs[i,3], i, adj=3)
## }


## pp <- par3d(no.readonly=TRUE)
## write_rds(par3d()$userMatrix, "characterizeArchetypes_3Dplot_setting.rds")

rgl.viewpoint(userMatrix=read_rds("characterizeArchetypes_3Dplot_setting.rds"),
              fov=30, zoom=0.8)
## rgl.postscript("characterizeArchetypes_3Dplot.svg", fmt="svg");
## rgl.snapshot("characterizeArchetypes_3Dplot.png", fmt="png")

noNormal <- (clinical$sample_type == "Primary Tumor")
clinical <- clinical[noNormal,]
mutations <- mutations[noNormal,]
clinical$cancer_type

cancerID <- clinical$cancer_type[1]
distToArchs <-
    sapply(1:nArchs, function(i) {
        sapply(1:nrow(samples), function(j) {
            sqrt(sum((samples[j,1:(nArchs-1)] - arcs[i,])^2))
        })
    })
binPct <- .05;
binPct <- .1;
archVsCancerTypes <- sapply(unique(clinical$cancer_type), function(cancerID) {
    sel <- which(clinical$cancer_type == cancerID)
    ## nTop <- round(nrow(samples)/(nArchs+1))
    nTop <- round(nrow(samples) * binPct)
    ## What fraction of tumors of a given cancer type are among the
    ## binPct% tumors closest to each archetype?
    sapply(1:nArchs, function(i) {
        sum(clinical$cancer_type[ order(distToArchs[,i])[1:nTop] ] ==
            cancerID) /
            sum(clinical$cancer_type == cancerID)
    })
    ## Break-down of cancer types among the binPct% tumors closest to each
    ## archetype?
    ## sapply(1:nArchs, function(i) {
    ##     sum(clinical$cancer_type[ order(distToArchs[,i])[1:nTop] ] ==
    ##         cancerID)
    ## })
    ## closestArch <- apply(distToArchs[sel,], 1, which.min)
    ## sapply(1:nArchs, function(i) {
    ##     sum(closestArch == i)
    ## }) / length(sel)
})

rownames(archVsCancerTypes) <- names(arcCols)
pdf("archVsCancerTypes.pdf", height=4, width=9);
par(mar=c(4, 5, 2, 1))
barplot(100*archVsCancerTypes, bes=T, legend=T, col=rep(arcCols, ncol(archVsCancerTypes)),
        ylab=sprintf(paste("fraction of tumors of a given type",
                           "among the %d%% tumors closest to archetype",
                           sep="\n"), binPct*100),
        ylim=c(0, 25))
abline(h=100*binPct, lty=2, lwd=2, col="grey")
dev.off();

ent <- apply(archVsCancerTypes, 2, function(v) {
    p <- v / sum(v)
    log2p <- log2(p)
    log2p[p==0] <- 0
    -sum(p * log2p)
})
pdf("cancerTypeVsArchetypeEntropies.pdf", height=4, width=9)
par(mar=c(4, 5, 2, 1))
bp <- barplot(ent, ylab="entropy (bits)", ylim=c(0,log2(nArchs)*1.2), col="lightblue")
abline(h=log2(nArchs), lty=2, lwd=2)
text(bp, 2.5, sprintf("%.1f", 2^ent))
dev.off();

## Number of effective archetypes
mean(ent) ## entropy of 2.1 bit on average
mean(2^ent) ## 4.3 archetypes per cancer type on average
round(range(2^ent), 1) ## between 3.4 and 5 archetypes per cancer type
mean(ent) / log2(nArchs) ## 90% of maximal entropy

## What archetypes are found in different cancer types?

cancerID <- "LGG"
pArchVsCancerTypes <- sapply(unique(clinical$cancer_type), function(cancerID) {
    ## sel <- which(clinical$cancer_type == cancerID)
    ## nTop <- round(nrow(samples)/(nArchs+1))
    nTop <- round(nrow(samples) * binPct)
    ## What fraction of tumors of a given cancer type are among the
    ## 10% tumors closest to each archetype?
    counts <-
        sapply(1:nArchs, function(i) {
            sum(clinical$cancer_type[order(distToArchs[,i])[1:nTop]] == cancerID)
        })
    pbinom(counts, size=sum(counts), prob=1/nArchs)
})
rownames(pArchVsCancerTypes) <- names(arcCols)
cancerTypesVsArchMat <- t(pArchVsCancerTypes > .01)
cancerTypesVsArchMatF <- cancerTypesVsArchMat
cancerTypesVsArchMatF[cancerTypesVsArchMat] <- "X"
cancerTypesVsArchMatF[!cancerTypesVsArchMat] <- ""
## colnames(cancerTypesVsArchMatF) <- names(arcCols)

write.table(t(cancerTypesVsArchMatF), file="cancerTypesVsArchMat.tsv", sep="\t");

apply(t(pArchVsCancerTypes > .01/prod(dim(pArchVsCancerTypes))), 1, sum)

##################################################
## Make matrix of tissue archetypes vs super-archetypes

SAMSig <- read_csv("MSigDBenrichment_continuous_significant.csv")
## arcsMSig <- arcsMSig %>% filter(`P value (Mann-Whitney)` < 1e-3)
SAMSig <- SAMSig %>% filter(`Mean Difference` > 0.1)

cancerIDs <-
    read_tsv("../TCGA_frac_nArchs.tab", col_names=F) %>% .[,1] %>%
    unlist %>% as.character %>%
    setdiff(c("HNSC", "LUAD", "BRCA"))
arcsMSigFiles <-
    c(paste("../", cancerIDs,
            "_UCSC/MSigDBenrichment_continuous_significant.csv",
            sep=""),
      "~/work/cancerTaskAtlas/brca_metabric/MSigDBenrichment_continuous_significant.csv")
cancerIDs <- c(cancerIDs, "BRCA")
cbind(cancerIDs, arcsMSigFiles)

cancerIdx <- 1
arcsMSigL <-
    map(1:length(arcsMSigFiles), function(cancerIdx) { # iterate over tissues
        arcsMSig <-
            read_csv(arcsMSigFiles[cancerIdx]) %>%
            filter(`Mean Difference` > 0.1)
        featsUniv <- union(SAMSig %>% select(`Feature Name`),
                           arcsMSig %>% select(`Feature Name`))

        arcIdx <- 1;
        ## Match each archetype to the super archetype with highest overlap in
        ## enriched MSigDBs
        SAmapping <- # iterate over tissue archetypes
            map(unlist(arcsMSig %>% select(`archetype #`) %>% unique), function(arcIdx) {
                SAidx <- 1;
                arcScores <- # iterate over super-archetypes
                    map(unlist(SAMSig %>% select(`archetype #`) %>% unique),
                        function(SAidx) {
                            SAfeats <- SAMSig %>%
                                filter(`archetype #` == SAidx) %>%
                                select(`Feature Name`)
                            arcFeats <- arcsMSig %>%
                                filter(`archetype #` == arcIdx) %>%
                                select(`Feature Name`)
                            ## nrow(SAfeats)
                            ## nrow(arcFeats)
                            
                            expIntersect <- nrow(SAfeats) * nrow(arcFeats) / nrow(featsUniv)
                            
                            ## phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
                            ## q,x: number of white balls drawn without
                            ## replacement from an urn which contains
                            ## both black and white balls.
                            ## m: the number of white balls in the urn.
                            ## n: the number of black balls in the urn.
                            ## k: the number of balls drawn from the urn.
                            p <- phyper(q=intersect(SAfeats, arcFeats) %>% nrow,
                                        m=SAfeats %>% nrow,
                                        n=nrow(featsUniv) - nrow(SAfeats),
                                        k=arcFeats %>% nrow,
                                        lower.tail=F)
                            
                            foldEnrich <- nrow(intersect(SAfeats, arcFeats)) / expIntersect;
                            ## return(foldEnrich)
                            ## return(p)

                            cutOff <- .01 /
                                ((arcsMSig %>% select(`archetype #`) %>%
                                  unique %>% nrow) *
                                 (SAMSig %>% select(`archetype #`) %>%
                                  unique %>% nrow))

                            return(c("obs"=nrow(intersect(SAfeats, arcFeats)),
                                     "exp"=expIntersect,
                                     "p"=p,
                                     "isSignif"=p<cutOff))
                        })
                arcScoresT <-
                    sapply(arcScores, function(x) { x }) %>% t %>% as.data.frame %>% 
                    rownames_to_column() %>%
                    mutate(tissueArch=arcIdx) %>%
                    rename("univArch"=rowname)
                return(arcScoresT)
            })
        SAmappingT <-
            SAmapping %>% bind_rows %>%
            mutate(cancer=cancerIDs[cancerIdx])
    })
arcsMSigT <- 
    arcsMSigL %>% bind_rows %>% 
    inner_join(
        tibble(univArch=paste("archetype #", 1:5, sep=""),
               cancerTask=names(arcCols))) %>%
    select(-univArch) %>%
    mutate(tissueArch=paste(cancer, tissueArch)) %>%
    select(-cancer) %>%
    mutate(enrichment=obs/exp, pp=-log10(p)) %>%
    as_tibble %>%
    mutate(isSignif=parse_logical(isSignif)) %>% 
    mutate(logEnrich=ifelse(enrichment<1 & !isSignif, 0, log2(enrichment))) %>%
    mutate(logEnrich=ifelse(logEnrich>3, 3, logEnrich)) %>% 
    mutate(pp=ifelse(isSignif, pp, 0)) %>%
    mutate(pp=ifelse(pp>80, 80, pp)) %>%
    arrange(desc(tissueArch)) %>%
    mutate(cancerTask=parse_factor(cancerTask,
                                   levels=rev(names(arcCols)[c(3,2,1,4,5)]),
                                   ordered=T))
    ## mutate(tissueArch=forcats::as_factor(tissueArch, levels=(sort(unique(tissueArch)))))

## Find best task for each archetype
bestTask <-
    sapply(arcsMSigT %>% select(tissueArch) %>% unique %>%
           unlist %>% as.character, function(arch) {
               signifSA <-
                   arcsMSigT %>% inner_join(tibble(tissueArch=!!arch)) %>%
                   filter(isSignif)
               if ( nrow(signifSA) == 0 ) { return(NA) }
               SAidx <- signifSA %>% select(p) %>% unlist %>% as.numeric %>% which.min
               signifSA[SAidx,"cancerTask"] %>% unlist %>% as.character
           }) %>% enframe %>%
    mutate(value=parse_factor(value,
                              levels=rev(names(arcCols)[c(3,2,1,4,5)]),
                              ordered=T))
       

ggplot(arcsMSigT) +
    geom_tile(aes(y=cancerTask, x=tissueArch, fill=pp)) +
    scale_fill_distiller(palette="GnBu") +
    geom_point(aes(x=name, y=value), shape=21, data=bestTask, size=3, fill="grey") +
    ## scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(fill="-log10 p", y="universal cancer task",
         x="tissue archetypes")
ggsave("matrixUnivsalVsTissueArchs.pdf", height=3, width=14)

#################################

## For each universal archetype, list tissue specific gene expression

SAidx = 1
tsMSigDBT = map(SAMSig %>% pull(`archetype #`) %>% unique, function(SAidx) { ## iterate over univ tasks
  defTasks = SAMSig #%>% filter(`archetype #` == !!SAidx) ## these are defining genes of this task
  tTask = arcCols[SAidx] %>% names

  tissueArchs = bestTask %>% filter(value == !!tTask) %>% pull(name)

  ta = tissueArchs[1]
  tsMSigDBsT = map(tissueArchs, function(ta) {
    pta = ta %>% str_split(" ") %>% .[[1]]
    tibble(`Tissue-specific MSigDB genes enriched at archetype`=
             arcsMSigFiles %>% .[which(cancerIDs == pta[1])] %>% 
             read_csv() %>% 
             filter(`Mean Difference` > 0.1) %>% 
             filter(`archetype #` == pta[2]) %>% 
             pull(`Feature Name`) %>% 
             setdiff(defTasks %>% pull(`Feature Name`))) %>% 
      mutate(archetype=!!tTask, tissue=!!pta[1]) %>% .[,c(2,3,1)]
  }) %>% bind_rows
}) %>% bind_rows

write_csv(tsMSigDBT, "univArchetypes_tissueSpecific_MSigDB.csv")

##################################################

## Can we find outliers by finding two regimes in sorted distances to
## center of gravity?

dists <- rev(sort(sqrt(apply(samples[,1:4]^2, 1, sum))));
plot(dists, type="l",
     ylab="distance to center of mass")
plot(diff(dists), type="l",
     ylab="distance to center of mass")

plot(diff(dists)[1:200], type="l",
     ylab="distance to center of mass")
abline(h=mean(diff(dists[100:200])))
abline(h=mean(diff(dists[100:200])) - 1.96*sd(diff(dists[100:200])), lty=2)
abline(v=50, lty=2)

toKeep <- rev(order(dists))[-(1:50)]
## samples <- samples[toKeep,]

##################################################

## plot3d(samples[,c(1, 2, 4)], alpha=.5, col="grey",
##        xlab="", ylab="", zlab="",
##        box=F, axes=F)
## spheres3d(arcs[,c(1, 2, 4)], col="blue", radius=10)
## i <- 1;
## for (i in 1:(nrow(arcs)-1)) {    
##     j <- 2;
##     for (j in 2:nrow(arcs)) {
##         lines3d(unlist(arcs[c(i,j),1]),
##                 unlist(arcs[c(i,j),2]),
##                 unlist(arcs[c(i,j),4]), col="blue")
##     }
## }
## for(i in 1:nrow(arcs)) {
##     text3d(arcs[i,1], arcs[i,2], arcs[i,4], i, adj=3)
## }

##################################################

library(gtools)

## Generate a perfect simplex
plot3d(rdirichlet(5e3, alpha=c(1,1,1,1)) +
       matrix(rnorm(5e3 * 4, 0, .2), ncol=4))
sPlots <-
    as_tibble(rdirichlet(5e3, alpha=c(1,1,1,1)) +
              matrix(rnorm(5e3 * 4, 0, .1), ncol=4)) %>%
    mutate(isArc=F, label=NA) %>% 
    bind_rows(tribble(~V1, ~V2, ~V3, ~V4, ~isArc, ~label,
                      0, 0, 0, 0, T, 1,
                      1, 0, 0, 0, T, 2,
                      0, 1, 0, 0, T, 3,
                      0, 0, 1, 0, T, 4,
                      0, 0, 0, 1, T, 5)) %>%
    rename(X1=V1, X2=V2, X3=V3, X4=V4) %>% 
    mutate(PC3=ntile(X3, n=3), PC4=ntile(X4, n=3))
arcs <- 
    (filter(sPlots, isArc) %>% select(X1, X2, X3, X4))

ggplot(sPlots) +
    geom_point(aes(x=X1, y=X2), shape=1, color="grey") +
    geom_density2d(aes(x=X1, y=X2)) +
    ## geom_point(aes(x=X1, y=X2), shape=1, color="blue") +
    geom_point(aes(x=X1, y=X2), data=filter(sPlots, isArc)) +
    geom_label_repel(aes(x=X1, y=X2, label=label),
                     data=filter(sPlots, isArc)) +
    labs(x="PC1", y="PC2") +
    facet_grid(PC4~PC3)
ggsave("samplesThirdiles_PC1_vs_PC2_synth.pdf", height=5, width=5)

## Let's represent the 4D data in thirdiles of PC3 and PC4
arcs <- read_csv("arcs_dims.csv", col_names=F)
sPlots <-
    samples %>% select(X1, X2, X3, X4) %>% mutate(isArc=F, label=NA) %>%
    bind_rows(mutate(arcs, isArc=T, label=1:5)) %>% 
    mutate(PC3=ntile(X3, n=3), PC4=ntile(X4, n=3))
arcsPlots <- sPlots %>% filter(isArc) %>% mutate(arcCol=arcCols[label]);

ggplot(sPlots) +
    geom_point(aes(x=X1, y=X2), shape=16, size=3, color="grey") +
    geom_density2d(aes(x=X1, y=X2), alpha=.4) +
    geom_point(aes(x=X1, y=X2, color=factor(label)), data=arcsPlots) +
    geom_label_repel(aes(x=X1, y=X2, label=label), data=arcsPlots) +
    labs(x="PC1", y="PC2") +
    scale_color_manual(values=arcCols) +
    facet_grid(PC4~PC3) +
    theme(legend.position="none")
ggsave("samplesThirdiles_PC1_vs_PC2.pdf", height=5, width=5)

## Let's project the 4D onto all 5 3D faces

## i <- 1;
## for (i in 1:nArchs) {
##     ## We remove archetype i and project on the space by the remaining archetypes
##     remainingArcs <- (1:5)[-i]
##     V <- matrix(NA, ncol(arcs), length(remainingArcs)-1)
##     j <- 2
##     for (j in 2:length(remainingArcs)) {
##         V[,j-1] <- 
##             as.numeric(arcs[remainingArcs[j],] - arcs[remainingArcs[1],])
##     }
##     V <- GramSchmidt(V[,1], V[,2], V[,3])

##     myProj <- apply(as.matrix(sPlots[,1:4]), 1, function(x) {
##         V %*% matrix(x - as.numeric(arcs[remainingArcs[1],]), ncol=1)
##     })
##     mySamples <- t(myProj)[!sPlots$isArc,];
##     myArcs <- t(myProj)[sPlots$isArc,];

##     open3d()
##     plot3d(mySamples, alpha=.5, col="grey",
##            xlab="", ylab="", zlab="",
##            box=F, axes=F)
##     spheres3d(myArcs, col="blue", radius=diff(range(myArcs))/20)
##     j <- 1;
##     for (j in 1:(nrow(myArcs)-1)) {    
##         k <- 2;
##         for (k in 2:nrow(myArcs)) {
##             lines3d(unlist(myArcs[c(j,k),1]),
##                     unlist(myArcs[c(j,k),2]),
##                     unlist(myArcs[c(j,k),3]), col="blue")
##         }
##     }
##     for(j in 1:nrow(arcs)) {
##         text3d(myArcs[j,1], myArcs[j,2], myArcs[j,3], j, adj=3)
##     }    
## }

##################################################
## Let's project the 4D onto all 10 2D faces

## Highlight BC and LGG
library(forcats)
ctHighlight <- clinical$cancer_type
## ctHighlight[ctHighlight != "BRCA" & ctHighlight != "LGG"] <- "other"
ctHighlight[ctHighlight != "BRCA" & ctHighlight != "LGG" & ctHighlight != "LIHC"] <- "other"
ctHighlight <- as_factor(ctHighlight)

## par(mfrow=c(2,5))
i <- 1;
## i <- 2;
for (i in 1:(nArchs-1)) {
    j <- 3;
    for (j in (i+1):nArchs) {
        ## We remove archetype i and project on the space by the remaining archetypes
        remainingArcs <- (1:5)[-c(i,j)]
        V <- matrix(NA, ncol(arcs), length(remainingArcs)-1)
        k <- 2
        for (k in 2:length(remainingArcs)) {
            V[,k-1] <- 
                as.numeric(arcs[remainingArcs[k],] - arcs[remainingArcs[1],])
        }
        V <- GramSchmidt(V[,1], V[,2], V[,2])[1:2,]
        
        ## V[,1] <- V[,1] / sqrt(sum(V[,1]^2))
        ## V[,2] <- V[,2] - sum(V[,1] * V[,2]) * V[,1]
        ## V[,2] <- V[,2] / sqrt(sum(V[,2]^2))       

        myProj <- apply(as.matrix(sPlots[,1:4]), 1, function(x) {
            V %*% matrix(x - as.numeric(arcs[remainingArcs[1],]), ncol=1)
        })
        mySamples <- t(myProj)[!sPlots$isArc,];
        myArcs <- t(myProj)[sPlots$isArc,];
        myArcs <- as_tibble(myArcs) %>% mutate(label=1:nrow(myArcs));

        ptSize <- rep(.5, length(ctHighlight))
        ## ptSize[ctHighlight == "other"] <- .4;

        ggplot(as_tibble(mySamples) %>%
               mutate(cType=ctHighlight, ptSize=ptSize) %>%
               .[sample(1:nrow(mySamples), size=nrow(mySamples)),]
               ) +
            geom_point(aes(x=V1, y=V2, color=cType, size=ptSize), alpha=1, shape=16) +
            geom_path(aes(x=V1, y=V2), color="grey",
                      data=myArcs[c(remainingArcs, remainingArcs[1]),]) +
            geom_point(aes(x=V1, y=V2, color=factor(label)),
                       data=myArcs[remainingArcs,], size=5) +
            scale_color_manual(values=c(as.character(arcCols[remainingArcs]),
                                        "#ffba00", "#ff8000",
                                        "#ff2d00", "grey")) +
            theme_void() +
            theme(legend.position="none")
        ggsave(sprintf("faceProjections_%d_%d.pdf", i, j), height=4, width=4)
        
        ggplot(as_tibble(mySamples) %>%
                 mutate(cType=ctHighlight, ptSize=ptSize) %>%
                 .[sample(1:nrow(mySamples), size=nrow(mySamples)),]
        ) +
          geom_point(aes(x=V1, y=V2, color="black", size=ptSize), alpha=1, shape=16) +
          geom_density2d(aes(x=V1, y=V2)) +
          geom_path(aes(x=V1, y=V2), color="grey",
                    data=myArcs[c(remainingArcs, remainingArcs[1]),]) +
          geom_point(aes(x=V1, y=V2, color=factor(label)),
                     data=myArcs[remainingArcs,], size=5) +
          scale_color_manual(values=c(as.character(arcCols[remainingArcs]),
                                      "#ffba00", "#ff8000",
                                      "#ff2d00", "grey")) +
          theme_void() +
          theme(legend.position="none")
        ggsave(sprintf("faceProjectionsDensity_%d_%d.pdf", i, j), height=4, width=4)

        ## ggplot(as_tibble(mySamples)) +
        ##     geom_point(aes(x=V1, y=V2), color="grey", size=3) +
        ##     ## geom_density2d(aes(x=X1, y=X2), color="blue", alpha=.4) +
        ##     ## geom_bin2d(aes(x=X1, y=X2)) +
        ##     geom_path(aes(x=V1, y=V2), color="grey",
        ##               data=myArcs[c(remainingArcs, remainingArcs[1]),]) +
        ##     geom_point(aes(x=V1, y=V2, color=factor(label)),
        ##                data=myArcs, size=5) +
        ##     ## geom_label(aes(x=X1, y=X2, label=label),
        ##     ##            nudge_x=30, data=myArcs) +
        ##     scale_color_manual(values=as.character(arcCols)) +
        ##     theme_void() +
        ##     theme(legend.position="none")
        
        library(forcats)
        if ( i == 1 && j == 3 ) {
            mycType <- "LIHC"
            sapply(unique(clinical$cancer_type), function(mycType) {
                dfFace <- data.frame(
                    mySamples,
                    cType=factor(clinical[,"cancer_type"] == mycType)) %>%
                    as_tibble %>%
                    mutate(cType=fct_recode(cType, Yes="TRUE", No="FALSE"))
                ## unlist(dfFace[,"cType"])
                ggplot(dfFace) +
                    geom_point(aes(x=X1, y=X2), color="grey", size=3,
                               data=filter(dfFace, cType == "No")) +
                    geom_point(aes(x=X1, y=X2), color="red", size=3, alpha=.3,
                               data=filter(dfFace, cType == "Yes")) +
                    geom_path(aes(x=V1, y=V2), color="grey",
                              data=myArcs[c(remainingArcs, remainingArcs[1]),]) +
                    geom_point(aes(x=V1, y=V2, color=factor(label)),
                               data=myArcs, size=5) +
                    scale_color_manual(values=c(as.character(arcCols), "grey", "#ff5555")) +
                theme_void() +
                    theme(legend.position="none")
                ggsave(sprintf("faceProjections_%d_%d_%s.pdf", i, j, mycType),
                       height=4, width=4)
            })
        }
        
        ## plot(mySamples, col="grey", xlab="", ylab="",
        ##      xlim=range(myArcs[,1]), ylim=range(myArcs[,2]))
        ## points(myArcs, pch=20, col="blue")
        ## k <- 1;
        ## for (k in 1:(nrow(myArcs)-1)) {    
        ##     l <- 2;
        ##     for (l in (k+1):nrow(myArcs)) {
        ##         lines(myArcs[c(k,l),1],
        ##               myArcs[c(k,l),2],
        ##               col="blue")
        ##     }
        ## }
        ## for(k in 1:nrow(arcs)) {
        ##     text(myArcs[k,1], myArcs[k,2], k, adj=c(1,1))
        ##     ## text(myArcs[k,1], myArcs[k,2], k, adj=c(1,-1/2))
        ## }
    }
}

## Let's project the 4D object onto all 4 3D PC projections

## i <- 1;
## for (i in 1:ncol(arcs)) {
##     remainingDims <- (1:ncol(arcs))[-i]

##     mySamples <- filter(sPlots, !isArc)[,remainingDims];
##     myArcs <- filter(sPlots, isArc)[,remainingDims];

##     open3d()
##     plot3d(mySamples, alpha=.5, col="grey",
##            xlab="", ylab="", zlab="",
##            box=F, axes=F)
##     spheres3d(myArcs, col="blue", radius=diff(range(myArcs))/50)
##     j <- 1;
##     for (j in 1:(nrow(myArcs)-1)) {    
##         k <- 2;
##         for (k in 2:nrow(myArcs)) {
##             lines3d(unlist(myArcs[c(j,k),1]),
##                     unlist(myArcs[c(j,k),2]),
##                     unlist(myArcs[c(j,k),3]), col="blue")
##         }
##     }
##     for(j in 1:nrow(arcs)) {
##         text3d(myArcs[j,1], myArcs[j,2], myArcs[j,3], j, adj=3)
##     }
## }

## Let's project the 4D data onto all 6 2-PC spaces

i <- 1;
for (i in 1:(ncol(arcs)-1)) {
    j <- 2;
    for (j in (i+1):ncol(arcs)) {
        ## We remove dimension i and j to project on the space by the
        ## remaining archetypes
        remainingArcs <- (1:ncol(arcs))[-c(i,j)]

        mySamples <- filter(sPlots, !isArc)[,remainingArcs]
        colnames(mySamples) <- c("V1", "V2")
        myArcs <- filter(sPlots, isArc)[,remainingArcs]
        colnames(myArcs) <- c("V1", "V2")
        
        ggplot(as_tibble(mySamples)) +
            geom_point(aes(x=V1, y=V2), color="grey", size=3) +
            geom_density2d(aes(x=V1, y=V2)) +
            geom_point(aes(x=V1, y=V2), data=as_tibble(myArcs)) +
            geom_label_repel(aes(x=V1, y=V2, label=arc),
                             data=as_tibble(myArcs) %>%
                                 mutate(arc=1:5)) +
            xlab(sprintf("PC %d", remainingArcs[1])) +
            ylab(sprintf("PC %d", remainingArcs[2]))
        ggsave(sprintf("PCprojections_%d_%d.pdf", i, j), height=5, width=5)
        
        ## plot(mySamples, col="grey", xlab="", ylab="",
        ##      xlim=range(myArcs[,1]), ylim=range(myArcs[,2]))
        ## points(myArcs, pch=20, col="blue")
        ## k <- 1;
        ## for (k in 1:(nrow(myArcs)-1)) {    
        ##     l <- 2;
        ##     for (l in (k+1):nrow(myArcs)) {
        ##         lines(myArcs[c(k,l),1],
        ##               myArcs[c(k,l),2],
        ##               col="blue")
        ##     }
        ## }
        ## for(k in 1:nrow(arcs)) {
        ##     text(myArcs[k,1], myArcs[k,2], k, adj=c(1,1))
        ##     ## text(myArcs[k,1], myArcs[k,2], k, adj=c(1,-1/2))
        ## }
    }
}

##################################################

## Show TP53 vector with the 5-simplex

dimSel <- c(1, 2, 4);
dimSel <- c(1, 2, 4);
dimSel <- c(2, 3, 4);
ptCol <- rep("grey", nrow(mutations))
ptCol[mutations[,"TP53"]==1] <- "red"
axLabs <- sprintf("PC%d", dimSel)
plot3d(samples[,dimSel], alpha=.95,
       col=ptCol, type="p", size=3,
       ## axes=T, xlab=axLabs[1], ylab=axLabs[2], zlab=axLabs[3],
       axes=F, xlab="", ylab="", zlab="",
       box=F, shininess=100,
       xlim=range(arcs[,dimSel[1]]),
       ylim=range(arcs[,dimSel[2]]),
       zlim=range(arcs[,dimSel[3]]))
spheres3d(arcs[,dimSel], col=arcCols, radius=10)
i <- 1;
for (i in 1:(nrow(arcs)-1)) {    
    j <- 2;
    for (j in 2:nrow(arcs)) {
        lines3d(unlist(arcs[c(i,j),dimSel[1]]),
                unlist(arcs[c(i,j),dimSel[2]]),
                unlist(arcs[c(i,j),dimSel[3]]), col="grey")
    }
}
for(i in 1:nrow(arcs)) {
    text3d(arcs[i,dimSel[1]], arcs[i,dimSel[2]], arcs[i,dimSel[3]], i, adj=5)
}

vec <-
    apply(samples[ptCol=="red",dimSel], 2, mean) -
    apply(samples[ptCol=="grey",dimSel], 2, mean)
offset <- c(-100,-100,100)
arrow3d(offset, offset+vec*3)

## par(mfrow=c(2,5))
i <- 2; j <- 5;
i <- 1; j <- 2;

## We remove archetype i and j and project on the space by the
## remaining archetypes
remainingArcs <- (1:5)[-c(i,j)]
V <- matrix(NA, ncol(arcs), length(remainingArcs)-1)
k <- 2
for (k in 2:length(remainingArcs)) {
    V[,k-1] <- 
        as.numeric(arcs[remainingArcs[k],] - arcs[remainingArcs[1],])
}
V <- GramSchmidt(V[,1], V[,2], V[,2])[1:2,]

 ## V[,1] <- V[,1] / sqrt(sum(V[,1]^2))
 ## V[,2] <- V[,2] - sum(V[,1] * V[,2]) * V[,1]
 ## V[,2] <- V[,2] / sqrt(sum(V[,2]^2))

myProj <- apply(as.matrix(sPlots[,1:4]), 1, function(x) {
    V %*% matrix(x - as.numeric(arcs[remainingArcs[1],]), ncol=1)
})
mySamples <- t(myProj)[!sPlots$isArc,];
myArcs <- t(myProj)[sPlots$isArc,];

mySamples <- bind_cols(as_tibble(mySamples), mutations[,"TP53"]) %>%
    mutate(TP53fac=factor(TP53))

vec <-
    unlist(10 *
           ( filter(mySamples, TP53 == 1) %>% summarize(x=mean(V1), y=mean(V2)) -
             filter(mySamples, TP53 == 0) %>% summarize(x=mean(V1), y=mean(V2)) ))
offset <- c(x=50, y=20)
mutVec <-
    tribble(~x, ~y, ~xend, ~yend,
            offset["x"], offset["y"], offset["x"] + vec["x"], offset["y"] + vec["y"])

ggplot(mySamples) +
    ## geom_point(aes(x=V1, y=V2), color="grey") +
    geom_density2d(aes(x=V1, y=V2, group=TP53fac, color=TP53fac),
                   data=drop_na(mySamples)) +
    ## geom_bin2d(aes(x=V1, y=V2)) +
    geom_path(aes(x=V1, y=V2), color="grey",
              data=as_tibble(myArcs)[c(remainingArcs, remainingArcs[1]),]) +
    geom_point(aes(x=V1, y=V2), data=as_tibble(myArcs),
               color="blue", size=5) +
    geom_label_repel(aes(x=V1, y=V2, label=arc),
                     data=as_tibble(myArcs) %>%
                         mutate(arc=1:5)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend), data=mutVec,
                 arrow=arrow(angle = 30, length = unit(0.25, "inches"))) +
    theme_void()
ggsave(sprintf("TP53_faceProjections_%d_%d.pdf", i, j), height=5, width=5)

## plot(mySamples, col="grey", xlab="", ylab="",
##      xlim=range(myArcs[,1]), ylim=range(myArcs[,2]))
## points(myArcs, pch=20, col="blue")
## k <- 1;
## for (k in 1:(nrow(myArcs)-1)) {    
##     l <- 2;
##     for (l in (k+1):nrow(myArcs)) {
##         lines(myArcs[c(k,l),1],
##               myArcs[c(k,l),2],
##               col="blue")
##     }
## }
## for(k in 1:nrow(arcs)) {
##     text(myArcs[k,1], myArcs[k,2], k, adj=c(1,1))
##     ## text(myArcs[k,1], myArcs[k,2], k, adj=c(1,-1/2))
## }



## Now examine gene expression, clinical features and mutations for
## those cancer types.

##################################################

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
              group=`Feature Name`), size=1.5)

source("../ParTI/TCGAscripts/hallmarkOmeter.inc.R")
toShow <-
    as_tibble(unlist(sapply(names(toShow), function(x) {
        sprintf("%s,%s", x, toShow[[x]])
    }))) %>% separate(value, into=c("Group", "Feature Name"), sep=",")

MSigDBord <- read_rds("../MSigDBord.rds");
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
    geom_hline(aes(yintercept=0), linetype=2)
ggsave("hallmarkOmeter.pdf", height=10, width=14);

## What's really different between archetypes 4 and 5?

MSigDBdiff <-
    full_join(
        MSigDB %>%
        filter(`archetype #` == 4 &
               `Significant after Benjamini-Hochberg correction?` == 1) %>%
        select(`Feature Name`, `Mean Difference`) %>%
        rename(mean4=`Mean Difference`),
        MSigDB %>%
        filter(`archetype #` == 5 &
               `Significant after Benjamini-Hochberg correction?` == 1) %>% 
        select(`Feature Name`, `Mean Difference`) %>%
        rename(mean5=`Mean Difference`)) %>%
    replace_na(list(mean4=0, mean5=0)) %>% 
    mutate(mean4=ifelse(mean4<0, 0, mean4)) %>%
    mutate(mean5=ifelse(mean5<0, 0, mean5)) %>%
    mutate(mean4vs5=mean4-mean5)
## MSigDBdiff %>% arrange(mean4vs5) %>% View()
MSigDBdiff %>% arrange(desc(mean4vs5)) %>% View()

##################################################

archProps[,"Zdiff"] <-
    archProps %>%
    select(-Group) %>% 
    dplyr::rename(arch=`archetype #`,
                  feat=`Feature Name`,
                  diff=`Mean Difference`) %>%
    pmap(function(arch, feat, diff) {
        ## browser()
        Zattr <-
            archProps %>% filter(`Feature Name` == feat) %>%
            summarize(sigma=sd(`Mean Difference`),
                      mu=mean(`Mean Difference`))
        (diff - Zattr[["mu"]])  / Zattr[["sigma"]]
    }) %>% unlist()

ggplot(archProps %>%
       mutate(`Feature Name`=factor(`Feature Name`, levels=MSigDBord))) +
    geom_col(aes(x=`Feature Name`, y=Zdiff)) +
    labs(y="log expression (Z-score)", x="MSigDB gene group")+ 
    coord_flip() +
    ## geom_errorbarh(aes(y=MSigDB,
    ##                    xmin=expression - SE,
    ##                    xmax=expression + SE)) +
    facet_grid(. ~ `archetype #`)
ggsave("hallmarkOmeter_byArch.pdf", height=6, width=14)

##################################################

## Show enrichment at different archetypes of MSigDBs representative
## of hallmarks
source("../ParTI/TCGAscripts/hallmarkOmeter.inc.R")

toShowFig <-
    as_tibble(unlist(sapply(names(toShowFig), function(x) {
        sprintf("%s,%s", x, toShowFig[[x]])
    }))) %>% separate(value, into=c("Group", "Feature Name"), sep=",")
featOrd <-
    unlist(toShowFig[,"Feature Name"]) %>%
    fmtMSigDBnames(maxWidth=40)

archProps <-
    inner_join(MSigDB, toShowFig) %>%
    ## mutate(`-log10 p (Mann-Whitney)`=
    ##            ifelse(`Significant after Benjamini-Hochberg correction?`,
    ##                   -log10(`P value (Mann-Whitney)`), 1)) %>% 
    mutate(`-log10 p (Mann-Whitney)`=
               ifelse(`Significant after Benjamini-Hochberg correction?` &
                      `Is first bin maximal?`,
                      -log10(`P value (Mann-Whitney)`), 1)) %>% 
    mutate(`Mean Difference`=
               ifelse(`Significant after Benjamini-Hochberg correction?` == 1 &
                     `Mean Difference` > 0,
                      `Mean Difference`, 0.01)) %>% 
    ## mutate(`-log10 p (Mann-Whitney)`=-log10(`P value (Mann-Whitney)`)) %>% 
    ## select(`archetype #`, `Feature Name`, `P value (Mann-Whitney)`, Group) %>%
    mutate(`Feature Name`=fmtMSigDBnames(`Feature Name`, maxWidth=40)) %>% 
    mutate(`Feature Name`=factor(`Feature Name`, levels=rev(featOrd))) %>% 
    inner_join(enframe(arcCols) %>%
               mutate(`archetype #`=1:length(arcCols)) %>%
               mutate(name=factor(name, levels=names(arcCols)))
               )

## mutate(`Feature Name`=factor(`Feature Name`, levels=MSigDBord))

ggplot(archProps) +
    geom_col(aes(x=`Feature Name`, y=`-log10 p (Mann-Whitney)`)) +
    ## geom_col(aes(x=`Feature Name`, y=`Mean Difference`)) +
    xlab("Functional gene group (MSigDB)") +
    ylab(expression(paste(-log[10], " enrichment p-value"))) +
    ## ylab("median fold-change") +
    coord_flip() +
    facet_grid(.~name)
ggsave("hallmarkOmeterSummarized.pdf", height=4, width=9);

##################################################

## Look at lipid metabolism genes

geneEnrich <- read_csv("geneEnrichment_continuous_All.csv")

source("~/work/jeanLib.R")
source("../ParTI/TCGAscripts/hallmarkOmeter.inc.R")
lipLevs <- as.character(unlist(lipidShow))
lipidShow <-
    list_to_tibble(lipidShow, key="Group", value="Feature Name")

archProps <-
    inner_join(geneEnrich, lipidShow) %>% 
    select(`archetype #`, `Feature Name`, `Mean Difference`, Group) %>%
    ## mutate(`Feature Name`=fmtMSigDBnames(`Feature Name`)) %>%
    mutate(`Feature Name`=factor(`Feature Name`, levels=rev(lipLevs))) %>% 
    mutate(`Group`=factor(Group))
##     ## complete(`archetype #`, `Feature Name`) %>%
## replace_na(list(`Mean Difference`=0))


ggplot(archProps) +
    geom_line(aes(x=`archetype #`, y=`Mean Difference`)) +
    ylab(expression(paste(log[2], " fold change"))) +
    facet_wrap(Group ~ `Feature Name` , scales="free_y") +
    geom_hline(aes(yintercept=0), linetype=2)

ggplot(archProps) +
    geom_col(aes(x=`Feature Name`, y=`Mean Difference`, fill=Group)) +
    labs(y="log expression", x="Gene")+ 
    coord_flip() +
    ## geom_errorbarh(aes(y=MSigDB,
    ##                    xmin=expression - SE,
    ##                    xmax=expression + SE)) +
    facet_grid(. ~ `archetype #`)

ggsave("lipids_byArch.pdf", height=3, width=15)
