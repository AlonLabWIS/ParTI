rm(list=ls())

drugDic <- read.csv("drugDictionnary.csv", as.is=T)

drugTargets <- read.csv("drugTargets3.csv", as.is=T)

drugHits <-
    sapply(unique(drugDic[,2]), function(x) {
        sum(x == drugTargets[,1])
    })

rev(sort(drugHits))

table(drugTargets[,3])

unique(unlist(sapply(c("angiogenesis", "Aurora Kinase", "Kinase inhibitor", "RAS",
         "Signaling"), function(x) {
             drugTargets[x == drugTargets[,3], "Main.target"]
         })))

setdiff(drugTargets[,3], c("angiogenesis", "Aurora Kinase", "Kinase inhibitor", "RAS",
         "Signaling", "transcription", "transcription, translation",
                           "Radiation", "cell division",
                           "immunosuprresant", "enhances Immune", "Support "))

unique(unlist(sapply(c("transcription", "transcription, translation", "Radiation"), function(x) {
             drugTargets[x == drugTargets[,3], "Main.target"]
         })))

unique(unlist(sapply(c("cell division"), function(x) {
             drugTargets[x == drugTargets[,3], "Main.target"]
         })))

unique(unlist(sapply(c("Support "), function(x) {
             drugTargets[x == drugTargets[,3], "Main.target"]
         })))

unique(unlist(sapply(c("immunosuprresant", "enhances Immune"), function(x) {
    drugTargets[x == drugTargets[,3], "Main.target"]
})))

idcs <-
    unique(unlist(sapply(c("", "Appoptosis", "Alternative", "omnipotency",
                    "proteasome"), function(x) {
                        which(x == drugTargets[,3])
                    })))

drugTargets[idcs,1:3]

## DEal with trailing whitspaces in 1st and 2nd column
length(gsub("[ ]*$", "", drugDic[,1]))
sort(table(gsub("[ ]*$", "", drugDic[,2])))



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
treatPerCancer <- treatPerCancer[apply(treatPerCancer, 1, sum) > 0,
                                 apply(treatPerCancer, 2, sum) > 0]

heatmap(log10(treatPerCancer), scale="none",
        hclustfun=function(x) { hclust(x, method="ward.D") })

hist(log10(apply(treatPerCancer, 1, mean)),
     xlab="log10 average usage frequency")
isFrequent <- apply(treatPerCancer, 1, mean) > .01
heatmap(treatPerCancer[isFrequent,], scale="none",
        hclustfun=function(x) { hclust(x, method="ward.D") })
