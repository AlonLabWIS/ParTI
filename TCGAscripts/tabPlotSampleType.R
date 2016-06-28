rm(list=ls())
library(LSD)

keepGenes <- .80; #Keep 33% top expressed genes

geneExpression <- read.csv("expMatrix.csv", as.is=T, h=F)
discrAttr <- read.table("discreteClinicalData_reOrdered.tsv", as.is=T,
                        sep="\t", h=T)
## attrList <- c("sample_type", "location_in_lung_parenchyma");
attrList <- c("tumor_tissue_site","sample_type");

## pdf("tabPlotSampleType.pdf", height=4, width=4*(length(attrList)+2));
## par(mfrow=c(1, length(attrList)+2));
pdf("tabPlotSampleType.pdf", height=3*2, width=3*3);
par(mfrow=c(2, 3))

par(mar=c(10, 4, 1, 1))
barplot(table(discrAttr[,attrList]), las=3, bes=F, leg=T)
avgExpr <- apply(geneExpression, 2, mean)
write.table(avgExpr, file="avgExpr.tsv", quote=F, row.names=F, col.names=F)

## plot(avgExpr, apply(geneExpression, 2, sd))

heatscatter(avgExpr, apply(geneExpression, 2, sd),
            xlab="log expression", ylab="inter-patients sd", main="");
abline(v=quantile(avgExpr, 1-keepGenes), lty=2, col="grey")
heatscatter(1:length(avgExpr),
            apply(geneExpression, 2, sd)[order(avgExpr)],
            xlab="expression (rank)", ylab="inter-patients sd", main="");
abline(v=(1-keepGenes)*length(avgExpr), lty=2, col="grey")
## plot(apply(geneExpression, 2, sd)[order(avgExpr)])

keptGenes <- rev(order(avgExpr))[1:round(length(avgExpr)*keepGenes)]
## summary(avgExpr); summary(avgExpr[keptGenes]);
keptGeneExpression <- geneExpression[,keptGenes];

library(ade4)
dudi1 <- dudi.pca(keptGeneExpression, nf=3, scannf=F, scale=F, center=T)
plot(100*(dudi1$eig[1:10] / sum(dudi1$eig)),
     ylab="% variance", xlab="# PC",
     ylim=c(0, max(100*(dudi1$eig[1:10] / sum(dudi1$eig)))));
rndESV <- sapply(1:11, function(r) {
    rndExpr <-
        apply(keptGeneExpression, 2, function(x) {
            sample(x, nrow(keptGeneExpression), replace=F)
        })
    dudir <- dudi.pca(rndExpr, nf=3, scannf=F, scale=F, center=T)
    100*(dudir$eig[1:10] / sum(dudir$eig));
})
lines(apply(rndESV, 1, mean), lty=2, col="grey")
## apply(rndESV, 1, sd) ## / sqrt(ncol(rndESV))
legend("topright", "shuffling", lty=2, col="grey")
## abline(h=100 / min(dim(keptGeneExpression)), lty=2)
attr <- attrList[2]
for ( attr in attrList ) {
    plot(dudi1$li[,1], dudi1$li[,2], xlab="PC1", ylab="PC2",
         xlim=range(dudi1$li[,1:2]), ylim=range(dudi1$li[,1:2]));
    vals <- unique(discrAttr[,attr]);
    sapply(1:length(vals), function(i) {
        mySel <- discrAttr[,attr] == vals[i]
        points(dudi1$li[mySel,1], dudi1$li[mySel,2],
               pch=20, col=i+1);
    })
    legend("topleft", vals, col=1+(1:length(vals)), pch=20)
}
dev.off();

##################################################
## Experiment with PCA

## geneExpression <- cbind(runif(1000) - .5, rnorm(1000, 0, .1/2),
##                         rnorm(1000, 0, .1/2), rnorm(1000, 0, .1/2))
## plot(geneExpression, xlim=range(geneExpression), ylim=range(geneExpression))

## dudi1 <- dudi.pca(geneExpression, nf=3, scannf=F, scale=F, center=T)
## plot(100*(dudi1$eig[1:10] / sum(dudi1$eig)),
##      ylab="% variance", xlab="# PC")
## plot(dudi1$li[,1], dudi1$li[,2], xlab="PC1", ylab="PC2",
##      xlim=range(dudi1$li[,1], dudi1$li[,2]),
##      ylim=range(dudi1$li[,1], dudi1$li[,2]))
