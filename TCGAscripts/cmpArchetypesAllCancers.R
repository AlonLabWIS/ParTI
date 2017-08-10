rm(list=ls())

library(gplots)
library(stringr)
library(mclust)
library(tidyverse)
source("ParTI/TCGAscripts/hallmarkOmeter.inc.R")

zMat <- function(mat, maxZ=2.5) {
    ## z-scores a matrix by row, thresholding the max z-score
    ## this helps compare archetypes, giving the same weight to
    ## different MSigDB groups (other the immune domminates
    ## everything), and prevents outliers (like the bile genes
    ## dominating one liver archetype, presumably one with low grade)

    mat <-
        t(apply(mat, 1,
                function(x) { (x - mean(x)) / sd(x) }))
    mat[mat > maxZ] <- maxZ;
    mat[mat < (-maxZ)] <- (-maxZ);
    return(mat)
}

cancersTab <-
    read_tsv("TCGA_frac_nArchs.tab", col_names=F) %>%
    rename(cancerType=X1, frac=X2, nArch=X3)

cancerIDs <- as.character(unlist(select(cancersTab, cancerType)))

ct <- cancerIDs[1]
MSigDBsig <-
    lapply(cancerIDs, function(ct) {
        read_csv(sprintf("%s_UCSC/MSigDBenrichment_continuous_significant.csv", ct))
    })
names(MSigDBsig) <- cancerIDs;

multip <- sort(table(unlist(sapply(MSigDBsig, function(tab) {
    tab %>% select(`Feature Name`) %>% unique()
}))))
plot(multip)
table(multip)

## If we take MSigDBs that are enriched in at least one arch with
## p<pCut in all cancers, what should pCut be to have ~25 MSigDBs left
## in the end
log10pCuts <- seq(-15, -10, len=10);
plot(log10pCuts,
     sapply(log10pCuts,
            function(log10pCut) {
                ## cat(sprintf("log10pCut = %.1f", log10pCut))
                isSignifInNcancers <-
                    sapply(MSigDBsig, function(tab) {
                        tab %>% filter(`P value (Mann-Whitney)` < 10^log10pCut) %>%
                            select(`Feature Name`) %>% unique()
                    }) %>% unlist() %>% table()
                sum(isSignifInNcancers == length(MSigDBsig))
            }), xlab="log10p cut-off", ylab="# MSigDBs categories", type="l")

isSignifInNcancers <-
    sapply(MSigDBsig, function(tab) {
        tab %>% filter(`P value (Mann-Whitney)` < 10^(-14)) %>%
            select(`Feature Name`) %>% unique() 
    }) %>% unlist() %>% table() 
isSignifInNcancers[isSignifInNcancers == length(MSigDBsig)]

MSigDBall <-
    lapply(cancerIDs, function(ct) {
        read_csv(
            sprintf("%s_UCSC/MSigDBenrichment_continuous_All.csv", ct)) %>%
            select(`archetype #`, `Feature Name`, `Mean Difference`) %>%
            mutate(cancerType=ct)
    })
names(MSigDBall) <- cancerIDs;

multip <- sort(table(unlist(sapply(MSigDBall, function(tab) {
    select(tab, `Feature Name`) %>% unique()
}))))
plot(multip)
length(multip)

tMSigDBall <-
    bind_rows(MSigDBall) %>%
    mutate(archName=sprintf("%s %d", cancerType, `archetype #`))

ggplot(tMSigDBall) +
    geom_tile(aes(x=archName, y=str_trunc(`Feature Name`, width=40),
                  fill=`Mean Difference`));

## Prepare table for heatmap
GOvsArchs <-
    tMSigDBall %>%
    select(archName, `Feature Name`, `Mean Difference`) %>%
    spread(archName, `Mean Difference`) %>%
    mutate(featName=`Feature Name` %>% fmtMSigDBnames())
GOvsArchsMat <-
    as.matrix(select(GOvsArchs, -featName, -`Feature Name`))
rownames(GOvsArchsMat) <- as.data.frame(GOvsArchs[,"Feature Name"])[,1]
colnames(GOvsArchsMat) <- colnames(GOvsArchs)[c(-1,-ncol(GOvsArchs))]

GOvsArchsMat0 <- zMat(GOvsArchsMat)
range(GOvsArchsMat)
range(GOvsArchsMat0)

## library(ade4)
## library(rgl)
## dudi1 <- dudi.pca(t(GOvsArchsMat0), center=F, scale=F, nf=18, scan=F)
## cumsum(dudi1$eig) / sum(dudi1$eig)
## plot3d(dudi1$li)
## plot(hclust(dist(dudi1$li), method="ward.D2"))
## plot(Mclust(dudi1$li, model="VVI", G=1:6))

pdf("cmpArchetypesAllCancers_allMSigDB.pdf", height=6, width=8);
## h0 <- heatmap(GOvsArchsMat, scale="row", labRow=F,
##               hclustfun=function(d) { hclust(d, method="ward.D2") })
h0 <- heatmap(GOvsArchsMat0,
              scale="none", labRow=F,
              hclustfun=function(d) { hclust(d, method="ward.D2") })
dev.off();

## GOvsArchsSmaller <-
##     GOvsArchsMat[rev(order(apply(GOvsArchsMat, 1, sd)))[1:20],]
## heatmap(GOvsArchsSmaller, scale="none",
##         hclustfun=function(d) { hclust(d, method="ward.D") })

toShow <-
    as_tibble(unlist(sapply(names(toShow), function(x) {
        sprintf("%s,%s", x, toShow[[x]])
    }))) %>%
    separate(value, into=c("Group", "Feature Name"), sep=",")

GOvsArchsSmaller <- GOvsArchsMat[as.data.frame(toShow)[,"Feature Name"],]
rownames(GOvsArchsSmaller) <-
    toShow %>%
    mutate(featName=`Feature Name` %>% fmtMSigDBnames()

           ) %>%
    select(featName) %>% unlist()

## heatmap(GOvsArchsSmaller, scale="row",
##         hclustfun=function(d) { hclust(d, method="ward.D") })

rowColNum <- as.numeric(as.factor(
    as.data.frame(toShow)[,"Group"]));
rowColPal <- rainbow(length(unique(rowColNum)))

GOvsA0 <- zMat(GOvsArchsSmaller)
## ## Z-transform by cancer type
## ct <- "BLCA";
## GOvsA0 <- 
##     lapply(unique(gsub(" [1-5]$", "", colnames(GOvsArchsSmaller))), function(ct) {
##         ret <-
##             zMat(GOvsArchsSmaller[,grep(ct,
##                                         colnames(GOvsArchsSmaller))])
##         as_tibble(t(ret))
##         ## rownames(ret) <- rownames(GOvsArchsSmaller)
##     }) %>% bind_rows %>% as.data.frame
## rownames(GOvsA0) <- colnames(GOvsArchsSmaller)
## GOvsA0 <- t(GOvsA0);
## GOvsA0[1:5,1:5]

range(GOvsArchsSmaller)
range(GOvsA0)

## library(ade4)
## library(rgl)
## ## dudi1 <- dudi.pca(t(GOvsArchsMat0), center=F, scale=F, nf=4, scan=F)
## dudi1 <- dudi.pca(t(GOvsA0), center=F, scale=F, nf=4, scan=F)
## cumsum(dudi1$eig) / sum(dudi1$eig)
## plot3d(dudi1$li)
## h <- hclust(dist(dudi1$li), method="ward.D2")
## plot(Mclust(dudi1$li, model="EII", G=1:6)$BIC)

## ## mat4dist <- t(GOvsA0);
## plot(Mclust(GOvsA0, model="EII", G=1:6)$BIC)
## featClust <- Mclust(GOvsA0, model="EII", G=4)
## ## View(enframe(featClust$classification) %>% arrange(value))
## i <- 1
## mat4dist <- sapply(unique(featClust$classification), function(i) {
##     apply(GOvsA0[names(featClust$classification)[
##         featClust$classification == i],], 2, mean)
## })
## plot(Mclust(mat4dist, model="VII", G=1:6)$BIC)
## arcClust <- Mclust(mat4dist, model="VII", G=4)
## order(arcClust$classification)
## table(arcClust$classification)

## h1 <- heatmap.2(GOvsA0[,order(arcClust$classification)],
##                 scale="none",
##                 dendrogram="row",
##                 Colv=F,
##                 hclustfun=function(d) { hclust(d, method="ward.D2") },
##                 ## dendrogram="column", Rowv=NA, RowSideColors=rowColPal[rowColNum],
##                 trace="none",
##                 col="heat.colors",
##                 ## col="heat.colors",
##                 ## breaks=seq(-max(abs(GOvsArchsSmaller)),
##                 ##            max(abs(GOvsArchsSmaller)), len=256),
##                 margins=c(5, 15))

## h <- hclust(dist(mat4dist), method="ward.D2")

## ## mat4dist[mat4dist<=0] <- 0;
## ## dudi1 <- dudi.pca(mat4dist, center=F, scale=F, nf=3, scan=F)
## ## cumsum(dudi1$eig) / sum(dudi1$eig)
## ## plot3d(dudi1$li)
## ## plot(Mclust(dudi1$li, model="VVI", G=1:6)$BIC)
## ## h <- hclust(dist(dudi1$li), method="ward.D2")

## plot(h)

h <- hclust(dist(t(GOvsA0)), method="ward.D2")

pdf("cmpArchetypesAllCancers.pdf", height=6, width=8)
h1 <- heatmap.2(GOvsA0,
                scale="none",
                Colv=as.dendrogram(h),
                hclustfun=function(d) { hclust(d, method="ward.D2") },
                ## dendrogram="column", Rowv=NA, RowSideColors=rowColPal[rowColNum],
                trace="none",
                col="heat.colors",
                ## col="heat.colors",
                ## breaks=seq(-max(abs(GOvsArchsSmaller)),
                ##            max(abs(GOvsArchsSmaller)), len=256),
                margins=c(5, 15))
dev.off()

write_rds(rownames(GOvsA0)[h1$rowInd], "MSigDBord.rds")

pdf("cmpArchetypesAllCancers_byFunctions.pdf", height=6, width=8)
h2 <- heatmap.2(GOvsA0,
                scale="none",
                hclustfun=function(d) { hclust(d, method="ward.D2") },
                dendrogram="column", Rowv=NA, RowSideColors=rowColPal[rowColNum],
                trace="none",
                col="heat.colors",
                ## col="heat.colors",
                ## breaks=seq(-max(abs(GOvsArchsSmaller)),
                ##            max(abs(GOvsArchsSmaller)), len=256),
                margins=c(5, 15))
dev.off()

## Determine the number of clusters:

## see how many clusters it finds to be optimal, set it to search for
## at least 1 model and up 10.
d_clust <- Mclust(t(as.matrix(GOvsArchsMat0)),
                  model="EII", G=1:6)
plot(d_clust, what="BIC") # 3-4 clusters on all MSigDBs

d_clust <- Mclust(t(as.matrix(GOvsA0)), ## model="VVI",
                  model="EII", G=1:7)

pdf("superArchetypes_Mclust_BIC_hallmarks.pdf", height=4, width=4)
plot(d_clust, what="BIC") # 4 clusters
dev.off();

myG <- which.max(d_clust$BIC)
myG <- 4

saBy <- "GMall"
saBy <- "hclust"
saBy <- "GMhallmarks"
for ( saBy in c("GMall", "GMhallmarks", "hclust")) {
    if ( str_detect(saBy, "^GM") ) {
        if ( saBy == "GMall" ) {
            d_clust <- Mclust(t(as.matrix(GOvsArchsMat0)), model="VVI", G=myG)
        } else {
            d_clust <- Mclust(t(as.matrix(GOvsA0)), model="VVI", G=myG)
        }
        saMemb <-
            sapply(unique(d_clust$classification), function(x) {
                names(d_clust$classification)[x == d_clust$classification]
            })
        archNames <-
            c("Cell division", #
              "Biomass & energy", #
              "Invasion & signaling", #
              ## "Immune interaction", #
              "Peroxisome activity") #
        names(saMemb) <- archNames
        write_rds(saMemb, sprintf("saMemb_%s.rds", saBy))

        saProfMean <- 
            sapply(saMemb, function(idcs) { apply(GOvsA0[,idcs], 1, mean) }) %>%
            as_tibble() %>%
            mutate(MSigDB=rownames(GOvsA0)) %>%
            gather(-MSigDB, key="archetype", value="expression")
        saProfSE <-
            sapply(saMemb, function(idcs) { apply(GOvsA0[,idcs], 1, function(x) {
                sd(x) / sqrt(length(x)) }) }) %>%
            as_tibble() %>%
            mutate(MSigDB=rownames(GOvsA0)) %>%
            gather(-MSigDB, key="archetype", value="SE")
    } else {
        nArchPerClust <- c(4, 4 + 9, 5, 4+5)
        archNames <- c("Biomass & energy",
                       "Peroxisome activity",
                       ## "Unknown",
                       "Invasion & signaling",
                       "Cell division")
        sum(nArchPerClust) == length(h1$colInd)
        
        ## Make hall-mark like profiles for each super-archetype
        saMemb <-
            apply(tibble(begin=c(1, cumsum(nArchPerClust) + 1)[1:length(nArchPerClust)],
                         end=cumsum(nArchPerClust)), 1, function(x) { x["begin"]:x["end"] })
        saProfMean <- 
            saMemb %>% 
            purrr::map(function(idcs) { apply(GOvsA0[,h1$colInd[idcs]], 1, mean) }) %>%
            setNames(archNames) %>%
            as_tibble() %>%
            mutate(MSigDB=rownames(GOvsA0)) %>%
            gather(-MSigDB, key="archetype", value="expression")
        saProfSE <- 
            saMemb %>% 
            purrr::map(function(idcs) {
                apply(GOvsA0[,h1$colInd[idcs]], 1, function(x) {
                    sd(x) / sqrt(length(x)) }) }) %>%
            setNames(archNames) %>%
            as_tibble() %>%
            mutate(MSigDB=rownames(GOvsA0)) %>%
            gather(-MSigDB, key="archetype", value="SE")
    }
    saProf <- inner_join(saProfMean, saProfSE) %>%
        mutate(MSigDB=factor(MSigDB,
                             levels=rownames(GOvsA0)[h1$rowInd]))
    ggplot(saProf) +
        geom_col(aes(x=MSigDB, y=expression)) +
        geom_errorbar(aes(x=MSigDB,
                          ymin=expression - SE,
                          ymax=expression + SE)) +
        coord_flip() +
        ## geom_errorbarh(aes(y=MSigDB,
        ##                    xmin=expression - SE,
        ##                    xmax=expression + SE)) +
        facet_grid(. ~ archetype)
    ggsave(sprintf("superArchetypeProfiles_%s.pdf", saBy),
           height=6, width=14)

    ## Make a table of cancer type vs super-archetype
    if ( saBy == "hclust" ) {
        saArchToArch <-
            saMemb %>%
            purrr::map(function(idcs) {
                ## inClust <- 
                colnames(GOvsA0)[h1$colInd[idcs]]
            }) %>%
            setNames(archNames)
    } else {
        saArchToArch <- saMemb
    }
    names(saArchToArch) %>%
        purrr::map(function(x) {
            tibble(sa=rep(x, length(saArchToArch[[x]])),
                   ca=gsub(" [1-5]$", "", saArchToArch[[x]])) }) %>%
        bind_rows %>% unique %>% mutate(mb="X") %>%
        spread(ca, mb, fill="") %>%
        as.data.frame() %>% 
        write_tsv(sprintf("superArchetype_vs_cancerType_%s.tsv", saBy))
}

## ## Maybe the Unknown archetype is a superposition of DNA replicat ion
## ## & mitosis and peroxisome activity!?

## ## What MSigDBs are high in the unknown cluster compared to other
## ## archetypes?
## saArchToArch[["Unknown"]]
## setdiff(flatten(saArchToArch), saArchToArch[["Unknown"]])

## g <- rownames(GOvsArchsMat)[1]
## unkDeltaP <-
##     t(sapply(rownames(GOvsArchsMat),
##              function(g) {
##                  x <- GOvsArchsMat[g,saArchToArch[["Unknown"]]]
##                  y <- GOvsArchsMat[g,setdiff(unlist(saArchToArch),
##                                              saArchToArch[["Unknown"]])]
##                  c(delta=mean(x) - mean(y), p=wilcox.test(x, y)$p.value)
##              }))
## head(unkDeltaP)
## ggplot(as_tibble(unkDeltaP) %>% mutate(log10p=log10(p))) +
##     geom_point(aes(delta,log10p))

## View(as_tibble(unkDeltaP) %>%
##      mutate(MSigDB=str_trunc(rownames(unkDeltaP), width=40)) %>% arrange(desc(delta)))
## View(as_tibble(unkDeltaP) %>%
##      mutate(MSigDB=str_trunc(rownames(unkDeltaP), width=50)) %>%
##      filter(delta>0) %>% arrange(p))

## ## Can't really find a group of genes that are specifically high in
## ## the unknown archetype.

## ## The unknown archetype could be a DNA replication & mitosis
## ## specialist which avoids triggering immune response -> earlier
## ## cancers?
## ## Clinical properties: in PRAD 2, PRAD 3, LGG 3, we find more
## ## tumor cells (5-10%) in this archetype. That's about it.
## saArchToArch[["Unknown"]]


## show table is similar when grouping by all genes / MSigDBs
## DONE

##################################################

## Visualize archetypes from different cancer types in 3D PCA
MSigDBall <-
    lapply(c("ALL", cancerIDs), function(ct) {
        read_csv(
            sprintf("%s_UCSC/MSigDBenrichment_continuous_All.csv", ct)) %>%
            select(`archetype #`, `Feature Name`, `Mean Difference`) %>%
            mutate(cancerType=ct)
    })
names(MSigDBall) <- c("ALL", cancerIDs);

multip <- sort(table(unlist(sapply(MSigDBall, function(tab) {
    select(tab, `Feature Name`) %>% unique()
}))))
plot(multip)
length(multip)

tMSigDBall <-
    bind_rows(MSigDBall) %>%
    mutate(archName=sprintf("%s %d", cancerType, `archetype #`))

ggplot(tMSigDBall) +
    geom_tile(aes(x=archName, y=str_trunc(`Feature Name`, width=40),
                  fill=`Mean Difference`));

## Prepare table for heatmap
GOvsArchs <-
    tMSigDBall %>%
    select(archName, `Feature Name`, `Mean Difference`) %>%
    spread(archName, `Mean Difference`) %>%
    mutate(featName=`Feature Name` %>% fmtMSigDBnames())
GOvsArchsMat <-
    as.matrix(select(GOvsArchs, -featName, -`Feature Name`))
rownames(GOvsArchsMat) <- as.data.frame(GOvsArchs[,"Feature Name"])[,1]
colnames(GOvsArchsMat) <- colnames(GOvsArchs)[c(-1,-ncol(GOvsArchs))]

## PCA
library(ade4)
library(rgl)
GOvsArchsMat[1:5,1:5]
dudi1 <- dudi.pca(t(GOvsArchsMat), center=F, scale=F, scannf=F, nf=3)
levels(as.factor(gsub(" .", "", rownames(dudi1$li))))

plot3d(dudi1$li, col=as.numeric(as.factor(gsub(" .", "", rownames(dudi1$li)))))
spheres3d(dudi1$li, col=as.numeric(as.factor(gsub(" .", "", rownames(dudi1$li)))))
text3d(dudi1$li, texts=gsub("^[^ ]*", "", rownames(dudi1$li)), adj=2)

cType <- levels(as.factor(gsub(" .", "", rownames(dudi1$li))))[1]
for (cType in levels(as.factor(gsub(" .", "", rownames(dudi1$li))))) {
    idcs <- grep(paste("^", cType, sep=""), rownames(dudi1$li))
    for (i in 1:(length(idcs)-1)) {
        for (j in 2:length(idcs)) {
            lines3d(
                dudi1$li[idcs[c(i,j)], 1],
                dudi1$li[idcs[c(i,j)], 2],
                dudi1$li[idcs[c(i,j)], 3],
                col=which(cType == levels(as.factor(gsub(" .", "", rownames(dudi1$li)))))
            )
        }
    }
}


plot(rep(0, length(levels(as.factor(gsub(" .", "",
                                         rownames(dudi1$li)))))),
     1:length(levels(as.factor(gsub(" .", "",
                                    rownames(dudi1$li))))),
     pch=20, cex=2, col=1:length(levels(as.factor(gsub(" .", "",
                                    rownames(dudi1$li))))))

text(rep(0.3, length(levels(as.factor(gsub(" .", "", rownames(dudi1$li)))))),
     1:length(levels(as.factor(gsub(" .", "", rownames(dudi1$li))))),
     levels(as.factor(gsub(" .", "", rownames(dudi1$li)))),
     col=1:length(levels(as.factor(gsub(" .", "", rownames(dudi1$li))))))

pdf("cmpArchetypesAllCancers_PCA2D.pdf", height=7, width=7);
plot(dudi1$li[,1:2],
     col=as.numeric(as.factor(gsub(" .", "", rownames(dudi1$li)))),
     pch=20, cex=2)
for (cType in levels(as.factor(gsub(" .", "", rownames(dudi1$li))))) {
    idcs <- grep(paste("^", cType, sep=""), rownames(dudi1$li))
    for (i in 1:(length(idcs)-1)) {
        for (j in 2:length(idcs)) {
            lines(
                dudi1$li[idcs[c(i,j)], 1],
                dudi1$li[idcs[c(i,j)], 2],
                lty=2,
                col=which(cType == levels(as.factor(gsub(" .", "", rownames(dudi1$li)))))
            )
        }
    }
}

text(dudi1$li[,1:2], labels=gsub("^[^ ]*", "", rownames(dudi1$li)),
     adj=2)
legend("topleft",
       levels(as.factor(gsub(" .", "", rownames(dudi1$li)))),
       col=1:length(levels(as.factor(gsub(" .", "", rownames(dudi1$li))))),
       pch=20)
dev.off()
