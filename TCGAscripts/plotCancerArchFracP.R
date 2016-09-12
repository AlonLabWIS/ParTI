rm(list=ls())

myTab <- read.table("cancerArchFracP.tsv", h=F, sep="\t", as.is=T)
cancers <- unique(myTab[,1])
nArchs <- unique(myTab[,2])

sortedCancers <-
    sort(sapply(cancers, function(x) {
        mean(log(myTab[myTab[,3] == .5 & myTab[,1] == x, 4] + .001))
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
