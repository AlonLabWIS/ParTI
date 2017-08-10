rm(list=ls())

library(tidyverse)

tag <- "_withNormal";
tag <- "_withoutNormal";
tag <- "_withoutNormal6";
tag <- "";

myTab <- read_tsv(sprintf("cancerArchFracP%s.tsv", tag), col_names=F) %>%
    rename(cancerType=X1, nArch=X2, quantile=X3, p=X4) %>%
    mutate(nArch=as.factor(nArch)) %>%
    filter(nArch != 2) %>%
    filter(nArch != 6) %>% filter(cancerType != "ALL_UCSC")
myTab <-
    inner_join(myTab, 
               myTab %>% group_by(cancerType) %>%
               summarize(avgLogP=mean(log(p + 1e-3)))) %>%
    arrange(avgLogP, nArch, quantile)
myTab <-
    myTab %>%
    mutate(cancerType=toupper(gsub(" UCSC", "", gsub("_", " ", cancerType))))

ggplot(myTab %>% mutate(p=p+1e-3) %>% rename(`# archetypes`=nArch)) +
    geom_line(aes(x=quantile, y=p, col=`# archetypes`, group=`# archetypes`)) +
    scale_y_log10() +
    labs(x="fraction of low-expressed gene excluded",
         y="p-value") +
    facet_wrap(~ cancerType)
ggsave(sprintf("cancerArchFracP%s.pdf", tag), height=6, width=10)

myTab %>%
    filter(p<.01/(10*3)) %>% group_by(cancerType, nArch) %>%
    summarize(n=n(),
              range=sprintf("[%.1f - %.1f]", min(quantile), max(quantile)),
              p=exp(mean(log(p+1e-3)))) %>% filter(n>=3) %>%
    arrange(cancerType, desc(n))

myTab %>% filter(p<.01/(10*3)) %>% group_by(cancerType, nArch) %>%
    summarize(n=n()) %>% filter(n>=3) %>% select(cancerType) %>% unique()

##################################################
## Previous version of same script

myTab <- read.table("cancerArchFracP.tsv", h=F, sep="\t", as.is=T)

cancers <- unique(myTab[,1])
nArchs <- unique(myTab[,2])

sortedCancers <-
    sort(sapply(cancers, function(x) {
        ## mean(log(myTab[myTab[,3] == .5 & myTab[,1] == x, 4] + .001))
        mean(log(myTab[myTab[,1] == x, 4] + .001))
    }))

cancers <- names(sortedCancers)
cancer <- cancers[1];
nArch <- nArchs[1];

pdf("cancerArchFracP_noHealthy.pdf", height=3*3, width=3*6);
par(mfrow=c(3, 6))
sapply(cancers, function(cancer) {
    plot(.1, .1, xlim=c(0, 100), ylim=c(1e-3, 1), log="y", type="n",
         xlab="cut-off (%)", ylab="p-value", main=cancer)
    sapply(nArchs, function(nArch) {
        thisTab <- myTab[myTab[,1] == cancer & myTab[,2] == nArch, ];
        lines(100*thisTab[,3], thisTab[,4] + 1e-3,
              col=which(nArch == nArchs), lwd=2)
    })
    abline(h=.05, lty=2, col="grey")
})
legend("bottomright", sprintf("%d arch.", nArchs),
       lty=1, col=1:length(nArchs), lwd=2);
dev.off();

## Cancers for which we can find a simplex
sCancers <- names(sortedCancers[1:10]);
sMyTab <- merge(data.frame(sCancers), myTab, by=1)

p <- .3
nArch <- 3;
avgPvals <- sapply(nArchs, function(nArch) {
    sapply(sort(unique(sMyTab[,3])), function(p) {
        exp(mean(log(sMyTab[sMyTab[,2] == nArch & sMyTab[,3] == p, 4] + 1e-3)))
    })
})
colnames(avgPvals) <- nArchs
barplot(t(avgPvals), bes=T, legend=T)

pdf("cutOffVsavgPval.pdf", height=5, width=5);
plot(.1, .1, xlim=c(0, 100), ylim=c(1e-3, 1), log="y", type="n",
     xlab="cut-off (%)", ylab="p-value (geom average)")
sapply(1:ncol(avgPvals), function(i) {
    lines(100*seq(0, .9, .1), avgPvals[,i], 
          col=i, lwd=2)
})
abline(h=.05, lty=2, col="grey")
legend("bottomright", sprintf("%d arch.", nArchs),
       lty=1, col=1:length(nArchs), lwd=2);
dev.off();
