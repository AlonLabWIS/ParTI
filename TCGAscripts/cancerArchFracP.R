rm(list=ls())

library(tidyverse)

tag <- "_20171113";
tag <- "_linSpace";
## tag <- "";

myTab <- read_tsv(sprintf("cancerArchFracP%s.tsv", tag), col_names=F) %>%
    rename(cancerType=X1, nArch=X2, quantile=X3, p=X4) %>%
    ## mutate(nArch=as.factor(nArch)) %>%
    ## filter(cancerType != "ALL_UCSC") %>% 
    filter(nArch != 2) %>%
    filter(nArch != 6)
    
myTab <-
    inner_join(myTab, 
               myTab %>% group_by(cancerType) %>%
               summarize(avgLogP=mean(log(p + 1e-3)))) %>%
    arrange(avgLogP, nArch, quantile)
myTab <-
    myTab %>%
    mutate(cancerType=toupper(gsub(" UCSC", "", gsub("_", " ", cancerType))))

ggplot(myTab %>% mutate(p=p+1e-3) %>% mutate(pp=-log10(p)) %>%
       ## filter(cancerType=="ALL" | nArch != 5) %>% 
       rename(`# archetypes`=nArch)) +
    geom_line(aes(x=quantile, y=p, col=`# archetypes`, group=`# archetypes`)) +
    scale_y_log10() +
    labs(x="fraction of low-expressed gene excluded",
         y="-log10(p)") +
        ## geom_abline(aes(intercept=-log10(.01), slope=0)) +
        ## geom_abline(aes(intercept=-log10(.05), slope=0), linetype="dashed") +
    facet_wrap(~ cancerType)
## ggsave(sprintf("cancerArchFracP%s.pdf", tag), height=6, width=10)

myTab %>%
    filter(p<.01/(10*3)) %>% group_by(cancerType, nArch) %>%
    summarize(n=n(),
              range=sprintf("[%.1f - %.1f]", min(quantile), max(quantile)),
              p=exp(mean(log(p+1e-3)))) %>% filter(n>=3) %>%
    arrange(cancerType, desc(n))

myTab %>% filter(p<.01/(10*3)) %>% group_by(cancerType, nArch) %>%
    summarize(n=n()) %>% filter(n>=3) %>% select(cancerType) %>% unique()

## At least once p<.01
myTab %>% filter(p<=.01) %>% filter(quantile<=.5&nArch!=5) %>%
    group_by(cancerType, nArch) %>%
    summarize(n=n()) %>% unique() %>% 
    filter(n>=1) %>%
    select(cancerType) %>% unique()

## At least once p<.05
myTab %>% filter(p<=.05) %>% filter(quantile<=.5&nArch!=5) %>%
    group_by(cancerType, nArch) %>%
    summarize(n=n()) %>% unique() %>% 
    filter(n>=1) %>% select(cancerType) %>% unique()

## p<.01 with all genes
myTab %>% filter(p<=.01) %>% filter(quantile==0) %>%
    group_by(cancerType, nArch) %>%
    summarize(n=n()) %>% unique() %>% 
    filter(n>=1) %>% select(cancerType) %>% unique()

## p<.05 with all genes
myTab %>% filter(p<=.05) %>% filter(quantile==0) %>%
    group_by(cancerType, nArch) %>%
    summarize(n=n()) %>% unique() %>% 
    filter(n>=1) %>% select(cancerType) %>% unique()

## Benjamini-Hochberg correction (fdr=10%)
myTabBH <-
    myTab %>% filter(quantile==0) %>%
    filter(cancerType != "ALL") %>% arrange(p)
library(forcats)
myTabBH <-
    myTabBH %>% mutate(pBHcutoff=.1*(1:nrow(myTabBH))/nrow(myTabBH)) %>%
    mutate(isSignifBH=factor(p<=pBHcutoff)) %>% 
    select(cancerType, nArch, p, isSignifBH)
myTabBH <-
    myTabBH %>%
    mutate(isSignifBH=fct_recode(isSignifBH,
                                 "Yes"="TRUE", "No"="FALSE")) %>% 
    arrange(cancerType, nArch)##  %>% 
    ## arrange(cancerType, p, desc(nArch))##  %>% 

## Final cancer types x nArch x p-value selection
myTabBH %>% filter(isSignifBH=="Yes")
configTab <-
    myTabBH %>% filter(isSignifBH=="Yes") %>%
    group_by(cancerType) %>%    
    summarize(nArch=min(nArch)) %>%
    inner_join(myTabBH)
configTab
## Write it out:
configTab %>% select(cancerType, nArch) %>%
    mutate(fracGenes=0) %>%
    filter(cancerType != "BRCA METABRIC") %>% 
    .[,c("cancerType", "fracGenes", "nArch")] %>%
    write_tsv(path="TCGA_frac_nArchs.tab", col_names=F)

## More stringent filter for cancers which we'll analyze together:
myTabBH %>% filter(isSignifBH=="Yes") %>%
    filter(p<.01) %>%
    select(cancerType) %>% unique()

fdrCutoff <- myTabBH %>% filter(isSignifBH=="Yes") %>% summarize(max(p))
ggplot(myTabBH %>% mutate(pp=-log10(p+1e-3)) %>%
       ## filter(cancerType=="ALL" | nArch != 5) %>% 
       rename(`# archetypes`=nArch)) +
    geom_point(aes(x=`# archetypes`, y=pp, col=isSignifBH)) +
    labs(x="Number of archetypes",
         y="-log10(p)",
         ## y="p-value (t-ratio test)",
         col="Significant at FDR<10%") +
    ## geom_abline(aes(intercept=fdrCutoff, slope=0), linetype="dashed") +
    ## geom_abline(aes(intercept=-log10(fdrCutoff), slope=0), linetype="dashed") +
    coord_cartesian(ylim=c(0,3)) + 
    ## scale_y_log10() +
    facet_wrap(~ cancerType)
ggsave(sprintf("cancerArchFracP%s_q0.pdf", tag), height=5, width=7)

##################################################
## Previous version of same script

## myTab <- read.table("cancerArchFracP.tsv", h=F, sep="\t", as.is=T)

## cancers <- unique(myTab[,1])
## nArchs <- unique(myTab[,2])

## sortedCancers <-
##     sort(sapply(cancers, function(x) {
##         ## mean(log(myTab[myTab[,3] == .5 & myTab[,1] == x, 4] + .001))
##         mean(log(myTab[myTab[,1] == x, 4] + .001))
##     }))

## cancers <- names(sortedCancers)
## cancer <- cancers[1];
## nArch <- nArchs[1];

## pdf("cancerArchFracP_noHealthy.pdf", height=3*3, width=3*6);
## par(mfrow=c(3, 6))
## sapply(cancers, function(cancer) {
##     plot(.1, .1, xlim=c(0, 100), ylim=c(1e-3, 1), log="y", type="n",
##          xlab="cut-off (%)", ylab="p-value", main=cancer)
##     sapply(nArchs, function(nArch) {
##         thisTab <- myTab[myTab[,1] == cancer & myTab[,2] == nArch, ];
##         lines(100*thisTab[,3], thisTab[,4] + 1e-3,
##               col=which(nArch == nArchs), lwd=2)
##     })
##     abline(h=.05, lty=2, col="grey")
## })
## legend("bottomright", sprintf("%d arch.", nArchs),
##        lty=1, col=1:length(nArchs), lwd=2);
## dev.off();

## ## Cancers for which we can find a simplex
## sCancers <- names(sortedCancers[1:10]);
## sMyTab <- merge(data.frame(sCancers), myTab, by=1)

## p <- .3
## nArch <- 3;
## avgPvals <- sapply(nArchs, function(nArch) {
##     sapply(sort(unique(sMyTab[,3])), function(p) {
##         exp(mean(log(sMyTab[sMyTab[,2] == nArch & sMyTab[,3] == p, 4] + 1e-3)))
##     })
## })
## colnames(avgPvals) <- nArchs
## barplot(t(avgPvals), bes=T, legend=T)

## pdf("cutOffVsavgPval.pdf", height=5, width=5);
## plot(.1, .1, xlim=c(0, 100), ylim=c(1e-3, 1), log="y", type="n",
##      xlab="cut-off (%)", ylab="p-value (geom average)")
## sapply(1:ncol(avgPvals), function(i) {
##     lines(100*seq(0, .9, .1), avgPvals[,i], 
##           col=i, lwd=2)
## })
## abline(h=.05, lty=2, col="grey")
## legend("bottomright", sprintf("%d arch.", nArchs),
##        lty=1, col=1:length(nArchs), lwd=2);
## dev.off();
