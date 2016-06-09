rm(list=ls())

cancerIDs <- read.table("TCGA_frac_nArchs.tab", as.is=T)[,1];

cancerID <- cancerIDs[1]
treatPerCancer <-
    sapply(cancerIDs, function(cancerID) {
        discrete <-
            read.table(
                sprintf("%s_UCSC/discreteClinicalData_reOrdered_withTreatment.tsv",
                        cancerID), as.is=T, h=T, sep="\t")
        justTreat <-
            discrete[,grep("^treat.", colnames(discrete))]
        apply(justTreat[,-1], 2, function(x) {
            mean(x, na.rm=T) })
    })

rownames(treatPerCancer) <- gsub("treat.", "", rownames(treatPerCancer))

dim(treatPerCancer)
which(apply(treatPerCancer, 1, sum) == 0)
treatPerCancer <- treatPerCancer[apply(treatPerCancer, 1, sum) > 0,
                                 apply(treatPerCancer, 2, sum) > 0]

heatmap(log10(treatPerCancer), scale="none",
        hclustfun=function(x) { hclust(x, method="ward.D") })

hist(log10(apply(treatPerCancer, 1, mean)),
     xlab="log10 average usage frequency")
isFrequent <- apply(treatPerCancer, 1, mean) > .01
heatmap(treatPerCancer[isFrequent,], scale="none",
        hclustfun=function(x) { hclust(x, method="ward.D") })
