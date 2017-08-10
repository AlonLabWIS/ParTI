rm(list=ls())

library(fdrtool)
library(gdata)
library(gplots)
library(ade4)
library(rgl)
library(tidyverse)

TCGAfracNarchs <- read.table("TCGA_frac_nArchs.tab", as.is=T, row.names=1)
colnames(TCGAfracNarchs) <- c("quantile", "nArchetypes");
TCGAfracNarchs["ALL",] <- c(0, 5)
cancerIDs <- rownames(TCGAfracNarchs);

cancerID <- cancerIDs[2]
cancerID <- "BRCA"
treatPerCancer <-
    sapply(cancerIDs, function(cancerID) {
        discrete <-
            read_tsv(
                sprintf("%s_UCSC/discreteClinicalData_reOrdered.tsv",
                        cancerID))
        treats <-
            discrete %>% filter(sample_type == "Primary Tumor") %>%
            .[,grep("^treat.", colnames(discrete))] %>%
            select(-treat.patientBestResp)
        apply(treats, 2, function(x) { mean(x, na.rm=T) })
         
        ## discrete <- 
        ##     read.table(
        ##         sprintf("%s_UCSC/treatTab.tsv",
        ##                 cancerID), as.is=T, h=T, sep="\t")[,c(-1, -2)]
        ## isPT <-
        ##     unlist(read_tsv(
        ##         sprintf("%s_UCSC/discreteClinicalData_reOrdered.tsv",
        ##                 cancerID))[,"sample_type"]) == "Primary Tumor"
        ## apply(discrete[isPT,], 2, function(x) { mean(x, na.rm=T) })
    })
rownames(treatPerCancer) <- gsub("^treat.", "", rownames(treatPerCancer))
##

## exprCancers <- 
##     lapply(cancerIDs, function(cancerID) {
##         cat(paste("Loading", cancerID, "expression\n"))
##         geneList <-
##             read_tsv(
##                 sprintf("%s_UCSC/geneListExp.list",
##                         cancerID), col_names=F)[,1]

##         expr <-
##             read_csv(sprintf("%s_UCSC/expMatrix.csv", cancerID),
##             col_names=F, col_types=paste(rep("d", nrow(geneList)), collapse="")) %>%
##             as.data.frame()
##         colnames(expr) <- as.character(unlist(geneList));
##         isPT <-
##             unlist(read_tsv(
##                 sprintf("%s_UCSC/discreteClinicalData_reOrdered.tsv",
##                         cancerID))[,"sample_type"]) == "Primary Tumor"
        
##         geneListAfterFiltering <-
##             read_tsv(
##                 sprintf("%s_UCSC/geneNamesAfterExprFiltering.list",
##                         cancerID), col_names=F)[,1]
        
##         return(expr[isPT,as.character(unlist(geneListAfterFiltering))])
##     })
## names(exprCancers) <- cancerIDs;
## save(exprCancers, file="exprCancers.rda")
load("exprCancers.rda")

##

exprPerCancer <-
    sapply(cancerIDs, function(cancerID) {
        return(
            apply(exprCancers[[cancerID]], 2, mean)
        )
    })

## archPerCancerList <-
##     sapply(cancerIDs, function(cancerID) {
##         cat(paste("Loading", cancerID, "archetypes\n"))
##         arcsOrig <-
##             t(as.data.frame(read_csv(
##                 sprintf("%s_UCSC/arcsOrig_genes.csv", cancerID), col_names=F)))
##         geneList <-
##             as.character(unlist(read_tsv(
##                 sprintf("%s_UCSC/geneNamesAfterExprFiltering.list",
##                         cancerID), col_names=F)))
##         rownames(arcsOrig) <- geneList;
##         return(arcsOrig)
##     })
## ## archPerCancerList <-
## ##     sapply(cancerIDs, function(cancerID) {
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
## save(archPerCancerList, file="archPerCancerList.rda")
load("archPerCancerList.rda")

## Keep only genes that are present in all the archetypes
geneTab <- table(unlist(sapply(archPerCancerList, function(x) { rownames(x) })))
## These genes are found in all archetypes of all cancer
ubiqGenes <- names(geneTab)[geneTab == length(cancerIDs)] 

alnArchs <-
    matrix(NA, length(ubiqGenes), 
           sum(sapply(archPerCancerList, function(a) { ncol(a) })));
rownames(alnArchs) <- ubiqGenes;
n <- 0;
for (i in 1:length(archPerCancerList)) {
    ## Normalize out abs gene expression inside each cancer
    ## archPerCancerList[[i]] <-
    ##     t(apply(archPerCancerList[[i]], 1, function(x) {x - mean(x)}))
    for (j in 1:ncol(archPerCancerList[[i]]) ) {
        n <- n + 1;
        alnArchs[,n] <- archPerCancerList[[i]][ubiqGenes,j]
    }
}

a <- names(archPerCancerList)[1]
colnames(alnArchs) <-
    unlist(sapply(names(archPerCancerList), function(a) {
        sprintf("%s.%d",
                rep(a, ncol(archPerCancerList[[a]])),
                1:ncol(archPerCancerList[[a]])
                )
    }))

## Compute optimal cut-off on patent response
respRanking <- c("Clinical Progressive Disease",
                 "Stable Disease", "Partial Response",
                 "Complete Response")

pdf("respCutOffs.pdf", height=3*4, width=4*4)
par(mfrow=c(3,4))
cancerID <- "BRCA"
respCutOffs <-
    sapply(cancerIDs, function(cancerID) {
        cat(paste(cancerID,"\n"))
        
        treatments <- 
            read.table(
                sprintf("%s_UCSC/discreteClinicalData_reOrdered.tsv",
                        cancerID), as.is=T, h=T, sep="\t")
        isPT <- treatments[,"sample_type"] == "Primary Tumor"
        treatments <- treatments[isPT,]
        
        treatments <- treatments[,c(which("treat.patientBestResp" == colnames(treatments)),
                                    grep("treat.target.", colnames(treatments)))]
        colnames(treatments) <- gsub("^treat.target.", "", colnames(treatments))
        colnames(treatments)[1] <- "best.resp";

        best.respVec <-
            sapply(treatments[,"best.resp"], function(x) {
                if ( is.na(x) ) { return(NA) }
                which(respRanking == x)
            })
        treatments[,"best.resp"] <- best.respVec;

        if ( ! any(!is.na(best.respVec)) ) { return(NA) }
        
        bp <- barplot(table(as.numeric(best.respVec)), main=cancerID)
        respCutOff <- median(best.respVec, na.rm=T);
        respCutOff <- respCutOff +
            which.min(abs(c(sum(best.respVec < respCutOff - .5, na.rm=T) / sum(!is.na(best.respVec)),
                            sum(best.respVec < respCutOff + .5, na.rm=T) /
                            sum(!is.na(best.respVec))) - .5)) - 1.5
        ## if ( respCutOff == 1 ) {
        ##     respCutOff <- 1.5
        ## } else if ( respCutOff == 4 ) {
        ##     respCutOff <- 3.5
        ## }
        abline(v=(respCutOff - 1)/(nrow(bp)-1) * diff(range(bp)) + min(bp), lty=2)
        return(respCutOff)
    })
dev.off();

save.image("drugVsCancers_init.rda")

## Done Init

load("drugVsCancers_init.rda")

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
## cancerID <- "SYNT";
EVsL <- list();
## load("allCancers_EVs.rda")

showEverything <- F;
responsesL <- list();

## treatPerCancerArch <- list();

cancerID <- "ALL";
cancerID <- "BLCA";

pdf("treatArchResponse.pdf", height=8, width=8);
for ( cancerID in setdiff(cancerIDs, "SYNT") ) {
    cat(paste(cancerID,"\n"))
    par(mfrow=c(2,2))

    expr <- exprCancers[[cancerID]];
    exprU <- t(expr[,rownames(archPerCancerList[[cancerID]])])
    
    treatments <- 
        read.table(
            sprintf("%s_UCSC/discreteClinicalData_reOrdered.tsv",
                    cancerID), as.is=T, h=T, sep="\t")
    isPT <- treatments[,"sample_type"] == "Primary Tumor"
    treatments <- treatments[isPT,]
    ## treatments <- treatments[,c(which("treat.patientBestResp" == colnames(treatments)),
    ##                             grep("treat.target.", colnames(treatments)))]
    ## colnames(treatments) <- gsub("^treat.target.", "", colnames(treatments))
    ## colnames(treatments)[1] <- "best.resp";
    if ( cancerID == "BRCA" ) {
        isERpos <- treatments[,"ER_Status_nature2012"] == "Positive"
    }
    treatments <- treatments[,grep("^treat", colnames(treatments))];
    colnames(treatments) <- gsub("^treat.", "", colnames(treatments))
    colnames(treatments)[grep("patientBestResp", colnames(treatments))] <- "best.resp";
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
        t(sapply(1:nrow(exprU), function(i) { #for all genes
            exprU[i,] - centerVec[i] })) ## set the mean of all genes to 0
    archs0 <- #same for archetypes
        t(sapply(1:nrow(exprU), function(i) {
            myArchs[i,] - centerVec[i]
        }))

    nPCs <- ncol(archs0) - 1;
    EVs <- svd(exprU0, nu=nPCs, nv=nPCs)$u[,1:nPCs];

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
    topFrac <- 1 / (ncol(archs0)+1)
    topPts <- round(nrow(exprProj) * topFrac)

    distsToArchs <-
        sapply(1:nrow(archProj), function(j) {
            sapply(1:nrow(exprProj), function(i) {
                sqrt(sum((exprProj[i,] - archProj[j,])^2))
            })
        })
    
    j <- 1;
    closestPtss <-
        sapply(1:nrow(archProj), function(j) {
            ## dists <- sapply(1:nrow(exprProj), function(i) {
            ##     sqrt(sum((exprProj[i,] - archProj[j,])^2))
            ## })
            closestPts <- order(distsToArchs[,j])[1:topPts];
            return(closestPts)
        })
    
    ## archFac <- rep(NA, nrow(exprProj))    
    ## for (i in 1:ncol(closestPtss)) {
    ##     ## FIXME the last archetypes snatched patients from other
    ##     ## archetypes
    ##     archFac[closestPtss[,i]] <- i;
    ## }
    ## archFac <- as.factor(archFac);
    ## summary(archFac)

    archMat <- matrix(F, nrow(exprProj), nrow(archProj))
    for (i in 1:ncol(closestPtss)) {
        archMat[closestPtss[,i],i] <- T
    }
    colnames(archMat) <- sprintf("closestTo%d", 1:ncol(archMat))

    minOcc <- (-1);
    treatOcc <- sort(apply(treatments[,-1], 2, function(x) { sum(x, na.rm=T) }))
    keptTreats <- names(which(treatOcc > minOcc))
    if ( length(keptTreats) == 0 ) { next; }
    treatFact <- as.data.frame(treatments[,keptTreats])
    ## for ( i in 1:ncol(treatFact) ) {
    ##     treatFact[,i] <- as.factor(treatFact[,i], levels=c(0,1))
    ## }
    ## table(treatments[,"best.resp"])

    respCutOff <- respCutOffs[cancerID]

    if ( cancerID == "BRCA" ) {
        ## Positive control for BRCA
        treated <- treatments[,"target.sex.hormone.signalling.inhibition"]
        treated[is.na(treated)] <- F
        plot(exprProj[,1], exprProj[,2],
             xlim=range(c(exprProj[,1], archProj[,1])),
             ylim=range(c(exprProj[,2], archProj[,2])),
             col=isERpos+2,
             main=cancerID)
        sel <- treatments[,"best.resp"] > respCutOff &
            treatments[,"target.sex.hormone.signalling.inhibition"]
        sel[is.na(sel)] <- F
        points(exprProj[sel,1], exprProj[sel,2],
               pch=20, col="red")
        points(archProj, pch=20, col="blue", cex=2)
        text(archProj,
             labels=c("DNA\nmitosis", "immune\ninvade", "peroxisome"),
             pos=1)

        ## People who are ER- have 45% chance to get sex hormone
        ## inhibitors vs 80% of ER+ people:
        table(treatments[,"target.sex.hormone.signalling.inhibition"], isERpos)
    }
    
    preGlmData <-
        data.frame(
            best.resp=as.numeric(treatments[,"best.resp"] > respCutOff),
            treatFact,
            ## arch=archFac,
            archMat)

    ## myTreat <- "tamoxifen"
    ## pValMat <- 
    ##     sapply(setdiff(colnames(treatFact), "best.resp"), function(myTreat) {
    ##         withTreat <- !is.na(preGlmData[,myTreat]) &
    ##             preGlmData[,myTreat] == 1
    ##         withBP <- !is.na(preGlmData[,"best.resp"]) &
    ##             preGlmData[,"best.resp"] == 1
    ##         notWithBP <- !is.na(preGlmData[,"best.resp"]) &
    ##             preGlmData[,"best.resp"] == 0
    ##         if ( sum(withTreat & withBP) < 5 ||
    ##              sum(withTreat & notWithBP) < 5 ) {
    ##             return(rep(NA, ncol(distsToArchs)))
    ##         }
    ##         sapply(1:ncol(distsToArchs), function(archIdx) {
    ##             ## wilcox.test(distsToArchs[withTreat & withBP, archIdx],
    ##             ##             distsToArchs[withTreat & !withBP, archIdx])
    ##             ## median(distsToArchs[withTreat & withBP, archIdx])
    ##             ## median(distsToArchs[withTreat & !withBP, archIdx])
    ##             wilcox.test(distsToArchs[withTreat & withBP, archIdx],
    ##                         distsToArchs[withTreat & notWithBP, archIdx],
    ##                         alt="less")$p.value
    ##         })
    ##     })
    ## pValMat[is.na(pValMat)] <- 1
    ## pValMat<.01
    ## sort(pValMat) < .1*(1:length(pValMat))/length(pValMat)
    ## ## Only one p-value left in the end. Probably not a good test to use.

    ## log OR test
    noNA <- !apply(is.na(preGlmData), 1, any)
    glmData <- preGlmData[noNA,]
    write_rds(glmData, sprintf("glmData_%s.rds", cancerID))

    summary(glmData)
    i <- "nucleotide.depletion";
    i <- "DNA.damaging";
    i <- "angiogenesis.signalling.inhibitor";
    i <- "microtubule";
    i <- "topoisomerase.inhibitor";
    
    responses <- lapply(colnames(treatFact), function(i) {
        cat(paste("Treatment", i, "\n"))
        ## tgd <- glmData[glmData[,i] == 1, c("best.resp","arch")]
        tgd <-
            glmData[glmData[,i] == 1,
                    c(grep("best.resp", colnames(glmData)),
                      grep("closestTo", colnames(glmData)))]

        logORSE <- matrix(NA, 3, nrow(archProj));
        rownames(logORSE) <- c("logOR", "SE", "p");
        ## colnames(logORSE) <- 1:ncol(logORSE)

        archIdx <- 1
        logORSE <-
            sapply(1:nrow(archProj), function(archIdx) {
                ## Focus on each archetype
                archLab <- sprintf("closestTo%d", archIdx);
                mytgd <- tgd[,c("best.resp", archLab)]
                if ( any(length(table(mytgd[,1])) < 2) |
                     any(length(table(mytgd[,2])) < 2) ) {
                    ## If patients who got treatment i are close to the
                    ## archetype, or all are far, we can't test if the
                    ## treatment has a different effect close to the
                    ## archetype. Same if all patients respond to
                    ## treatment.
                    return(c(logOR=NA, SE=NA, p=NA))
                }

                isClosest <- mytgd[,archLab]
                pc <- 1; #pseudo-count
                
                ORs <- c(
                  (sum(mytgd[isClosest,"best.resp"])+pc) /
                  (sum(!mytgd[isClosest,"best.resp"])+pc),
                  (sum(mytgd[!isClosest,"best.resp"])+pc) /
                  (sum(!mytgd[!isClosest,"best.resp"])+pc)
                )

                myct <- table(mytgd[,c("best.resp", archLab)]);
                
                return(
                    c(logOR=log(ORs[1] / ORs[2]),
                      SE=sqrt(sum(1/(myct+pc))),
                      p=fisher.test(myct)$p.value)
                )
            })

        ## If everything is NA, we're done.
        if ( all(is.na(logORSE["logOR",])) ) { return(logORSE) }
        
        if ( any(logORSE["p",] < .05, na.rm=T) || showEverything ) {
            bp <- barplot(logORSE["logOR",],
                          names=paste("A", 1:ncol(logORSE), sep=""),
                          main=sprintf("%s", i),
                          ylab=sprintf("log OR (n=%d)", nrow(tgd)))
            sapply(1:ncol(logORSE), function(j) {
                arrows(bp[j,1], logORSE["logOR",j] + 1.96 * logORSE["SE",j],
                       bp[j,1], logORSE["logOR",j] - 1.96 * logORSE["SE",j],
                       angle=90, code=3, length=0.05)
            })
        }
        
        return(logORSE);
    })
    names(responses) <- colnames(treatFact)
    responsesL[[cancerID]] <- responses;
}
dev.off();

## Print out results of enrichment analysis
sapply(names(responsesL), function(cancerID) {
    responsesL[[cancerID]]

    archIdx <- 1;
    sapply(1:ncol(responsesL[[cancerID]][[1]]), function(archIdx) {
        drug <- "carboplatin"
        sapply(names(responsesL[[cancerID]]), function(drug) {
            p <- responsesL[[cancerID]][[drug]]["p",archIdx];
            logOR <- responsesL[[cancerID]][[drug]]["logOR",archIdx];
            if ( !is.na(p) & p < .05) {
                neg <- "";
                if ( logOR < 0 ) { neg <- "in"; }                
                cat(sprintf("%s.%d is %s-sensitive to %s (logOR=%.1f, p=%.3f)\n", 
                            cancerID, archIdx, neg, drug, logOR, p))
            }
        })
    })
})

## Combine patients close to different archetypes in different cancer
## types
saMemb <- read_rds("saMemb_GMhallmarks.rds")
cancerID <- "BRCA"
cancerID <- "BLCA"
library(stringr)
allGlmData <-
    lapply(setdiff(cancerIDs, "ALL"), function(cancerID) {
        glmData <- as_tibble(read_rds(sprintf("glmData_%s.rds", cancerID)))
        sa <- "DNA replication & mitosis"
        for ( sa in names(saMemb) ) {
            locArcs <-
                saMemb[[sa]][str_detect(saMemb[[sa]], sprintf("^%s ", cancerID))]
            glmData[,sa] <- 
                apply(glmData[,gsub("^[A-Z]* ", "closestTo", locArcs)], 1, any)
        }
        glmData <- glmData[,-grep("^closestTo", colnames(glmData))]
        return(glmData)
    }) %>% bind_rows()

## Now analyse log odds ratios
allTreats <- setdiff(colnames(allGlmData), c("best.resp", names(saMemb)));
i <- "target.DNA.damaging"
pdf("treatArchResponse_polled.pdf", height=8, width=8);
par(mfrow=c(2,2))
responses <-
    lapply(allTreats,
           function(i) {
               cat(paste("Treatment", i, "\n"))
               tgd <-
                   allGlmData[allGlmData[,i] == 1,
                              c("best.resp", names(saMemb))]
               
               logORSE <- matrix(NA, 3, length(saMemb));
               rownames(logORSE) <- c("logOR", "SE", "p");
               colnames(logORSE) <- names(saMemb);
               
               archLab <- names(saMemb)[1]
               logORSE <-
                   sapply(names(saMemb), function(archLab) {
                       ## Focus on each archetype
                       mytgd <- tgd[,c("best.resp", archLab)] %>%
                           as.data.frame();
                       if ( any(length(table(mytgd[,1])) < 2) |
                            any(length(table(mytgd[,2])) < 2) ) {
                           ## If patients who got treatment i are close to the
                           ## archetype, or all are far, we can't test if the
                           ## treatment has a different effect close to the
                           ## archetype. Same if all patients respond to
                           ## treatment.
                           return(c(logOR=NA, SE=NA, p=NA))
                       }
                       
                       isClosest <- mytgd[,archLab]
                       pc <- 1; #pseudo-count
                       
                       ORs <- c(
                       (sum(mytgd[isClosest,"best.resp"])+pc) /
                       (sum(!mytgd[isClosest,"best.resp"])+pc),
                       (sum(mytgd[!isClosest,"best.resp"])+pc) /
                       (sum(!mytgd[!isClosest,"best.resp"])+pc)
                       )
                       
                       myct <- table(mytgd[,c("best.resp", archLab)]);
                       
                       return(
                           c(logOR=log(ORs[1] / ORs[2]),
                             SE=sqrt(sum(1/(myct+pc))),
                             p=fisher.test(myct)$p.value)
                       )
                   })

               ## If everything is NA, we're done.
               if ( all(is.na(logORSE["logOR",])) ) { return(logORSE) }
               
               if ( any(logORSE["p",] < .05, na.rm=T) || showEverything ) {
                   bp <- barplot(logORSE["logOR",],
                                 names=paste("A", 1:ncol(logORSE), sep=""),
                                 main=sprintf("%s", i),
                                 ylab=sprintf("log OR (n=%d)", nrow(tgd)))
                   sapply(1:ncol(logORSE), function(j) {
                       arrows(bp[j,1], logORSE["logOR",j] + 1.96 * logORSE["SE",j],
                              bp[j,1], logORSE["logOR",j] - 1.96 * logORSE["SE",j],
                              angle=90, code=3, length=0.05)
                   })
               }
               
               return(logORSE);
           })
names(responses) <- allTreats
dev.off();

pollSummary <-
    sapply(names(saMemb), function(sa) {
        sapply(names(responses), function(drug) {
            p <- responses[[drug]][["p", sa]]
            logOR <- responses[[drug]]["logOR",sa]
            if ( !is.na(p) & p < .05) {
                neg <- "";
                if ( logOR < 0 ) { neg <- "in"; }
                cat(sprintf("%s is %s-sensitive to %s (logOR=%.1f, p=%.3f)\n", 
                            sa, neg, drug, logOR, p))
            }
            (!is.na(p) & p < .05) * sign(logOR)
        })
    })


##################################################

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
dTreatN <- dTreat / median(dTreat, na.rm=T)
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
## load("allCancers_EVs.rda")
showEverything <- T;

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
    sapply(setdiff(names(respCutOffs)[!is.na(respCutOffs)], "SYNT"),
           function(cancerID) {
        cat(paste(cancerID,"\n"))
        par(mfrow=c(2,2))
        
        expr <- exprCancers[[cancerID]];
        myQuantile <- TCGAfracNarchs[cancerID,"quantile"];
        if ( is.na(myQuantile) || myQuantile == 0 ) {
            exprU <- t(expr);
        } else {
            ## meanGeneExpr <- apply(expr, 2, mean);
            ## minExpr <- quantile(meanGeneExpr, myQuantile);
            ## exprU <- t(expr[,meanGeneExpr >= minExpr])
            exprU <- t(expr[,rownames(archPerCancerList[[cancerID]])])
        }
        rm(expr);
        
        ## discreteClinicalData_reOrdered_withTreatment.tsv, select
        ## only treat.target
        treatmentsFull <- 
            read.table(
                sprintf("%s_UCSC/discreteClinicalData_reOrdered.tsv",
                        cancerID), as.is=T, h=T, sep="\t")
        isPT <- treatmentsFull[,"sample_type"] == "Primary Tumor"
        ## if ( any(colnames(treatmentsFull) == "pathologic_M") ) {
        ##     Mstatus <- treatmentsFull[,"pathologic_M"]
        ##     ## Nstatus <- treatments[,"pathologic_N"]
        ## } else {
        ##     cat(sprintf("No metastasis information for %s\n", cancerID))
        ##     Mstatus <- rep("M0", nrow(treatmentsFull))
        ## }
        treatments <-
            treatmentsFull[,c(which("treat.patientBestResp" == colnames(treatmentsFull)),
                              grep("treat.target.", colnames(treatmentsFull)))]
        colnames(treatments) <- gsub("^treat.target.", "", colnames(treatments))
        colnames(treatments)[1] <- "best.resp";

        best.respVec <-
            sapply(treatments[,"best.resp"], function(x) {
                if ( is.na(x) ) { return(NA) }
                which(respRanking == x)
            })
        treatments[,"best.resp"] <- best.respVec;
        treatments <- treatments[isPT,]
        
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
        nPCs <- ncol(archs0) - 1;
        ## if ( is.null(EVsL[[cancerID]]) ) {
            EVsL[[cancerID]] <- svd(exprU0, nu=nPCsInit, nv=nPCsInit)$u[,1:nPCsInit];
            save(EVsL, file="allCancers_EVs.rda")
        ## } else {
        ##     cat("Reusing previously computed eigenvectors.\n")
        ## }
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

        ## plot3d(exprProj[,1], exprProj[,2], exprProj[,3],
        ##        size=3, type="p", alpha=.5,
        ##        col="red", xlab="", ylab="", zlab="")
        ## points3d(archProj[,1], archProj[,2], archProj[,3],
        ##          col="blue", size=8)

        ## Compute distance from each tumor to each archetype
        distTtoA <-
            t(apply(exprProj, 1, function(tumor) {
                apply(archProj, 1, function(arch) {
                    sqrt(sum((tumor - arch)^2))
                })
            }))
        ## rownames(distTtoA) <-
        ##     read.table(sprintf("%s_UCSC/patientIDs.list", cancerID))[,1]
        write.csv(distTtoA,
                  file=sprintf("%s_dist_tumor_to_arch.csv", cancerID))
        
        ## isPRAD <- treatmentsFull[,"cancer_type"] == "BRCA"
        ## plot3d(exprProj[isPRAD,1], exprProj[isPRAD,2],
        ##        exprProj[isPRAD,3], size=.5, type="s",
        ##        col="red", xlab="", ylab="", zlab="")
        ## points3d(archProj[,1], archProj[,2], archProj[,3],
        ##          col="blue", size=8)
        ## rgl.points(exprProj[!isPRAD,1], exprProj[!isPRAD,2],
        ##            exprProj[!isPRAD,3], size=2, col="black", alpha=.25)
        
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

        ## if ( cancerID == "ALL" ) {
        ##     ## x <- unique(treatmentsFull[,"cancer_type"])[1]
        ##     ## isBestResp <- 
        ##     ##     unlist(
        ##     ##         sapply(unique(treatmentsFull[,"cancer_type"]), function(x) {
        ##     ##             ## sum(treatmentsFull[,"cancer_type"] == x)
        ##     ##             treatments[treatmentsFull[,"cancer_type"] == x,
        ##     ##                        "best.resp"] > respCutOffs[x]
        ##     ##         })
        ##     ##         )
        ##     isBestResp <- treatments[,"best.resp"] > median(respCutOffs, na.rm=T)
        ## } else {
            respCutOff <- respCutOffs[cancerID]
            isBestResp <- treatments[,"best.resp"] > respCutOff;
        ## }
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
## cancerID <- "GBM"
## sapply(colnames(pVals), function(cancerID) {
sapply(colnames(pVals), function(cancerID) {
    cat(paste(cancerID, "\n"))

    ## qMat <- pVals[["p", cancerID]]
    qMat <- pVals[["q",cancerID]]
    if ( is.null(qMat) ) { return() }
    
    i <- 1;
    sapply(1:nrow(qMat), function(i) {
        j <- 1;
        sapply(1:ncol(qMat), function(j) {
            if ( !is.na(qMat[i,j]) && qMat[i,j] < .1 ) {
                ## browser();
                ## cat(sprintf("%s: arch %d - %s (%.2f)\n",
                ##             cancerID, i, colnames(qMat)[j],
                ##             pVals[["q", cancerID]][i,j]))
                cat(sprintf("%s: arch %d - %s (dir=%.2f, p=%.3f)\n",
                            cancerID, i, colnames(qMat)[j],
                            pVals[["dir",cancerID]][i,j],
                            pVals[["p",cancerID]][i,j]))
            }

        })
        return()
    })
    return()
})
