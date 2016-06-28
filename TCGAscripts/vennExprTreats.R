rm(list=ls())

cancerIDs <- read.table("TCGA_frac_nArchs.tab", as.is=T)[,1]
cancerIDs <- cancerIDs[cancerIDs != "SYNT"]

cancerID <- cancerIDs[1]
venns <- 
    sapply(cancerIDs, function(cancerID) {
        withExpr <-
            read.table(sprintf("%s_UCSC/patientIDs.list", cancerID), as.is=T)[,1]
        withExpr <- gsub("-[0-9]*$", "", withExpr)

        treatTab <-
            read.table(sprintf("%s_UCSC/treatTab.tsv", cancerID), as.is=T,
                       sep="\t", h=T)
        withTreat <- treatTab[,1]
        
        c("both"=length(intersect(withExpr, withTreat)),
          "expr. only"=length(setdiff(withExpr,withTreat)),
          "treat. only"=length(setdiff(withTreat,withExpr)))
    })

pdf("vennExprTreats.pdf", height=5, width=7)
barplot(venns, leg=T, las=3, ylab="# patients",
        main="Treatment information is commonly missing from TCGA records")
dev.off();

cancerID <- cancerIDs[1]
venns <- 
    sapply(cancerIDs, function(cancerID) {
        withExpr <-
            read.table(sprintf("%s_UCSC/patientIDs.list", cancerID), as.is=T)[,1]
        withExpr <- gsub("-[0-9]*$", "", withExpr)

        treatTab <-
            read.table(sprintf("%s_UCSC/treatTab.tsv", cancerID), as.is=T,
                       sep="\t", h=T)
        withTreat <- treatTab[,1]
        withDocTreat <-
            treatTab[
                     apply(treatTab[,-1], 1, function(x) {
                         sum(is.na(x)) == 0 })
                     ,1]

        c("both"=length(intersect(withExpr, withDocTreat)),
          "expr. only"=length(setdiff(withExpr,withDocTreat)),
          "doc. treat. only"=length(setdiff(withDocTreat,withExpr)))
    })

pdf("vennExprDocTreats.pdf", height=5, width=7)
barplot(venns, leg=T, las=3, ylab="# patients",
        main="Treatment documentation is uncommon among TCGA records")
dev.off();

##################################################

## Typical treatment durations
pdf("typicalTreamentDurations.pdf", height=3*3, width=3*4)
par(mfrow=c(3,4))

cancerID <- cancerIDs[1]
treatDurations <-
    sapply(cancerIDs, function(cancerID) {
        tabLabels <-
            read.table(sprintf("%s_UCSC/drugFile.txt", cancerID),
                       h=F, sep="\t", as.is=T)[1:3,]
        drugs <- read.table(sprintf("%s_UCSC/drugFile.txt", cancerID),
                            h=F, sep="\t", skip=3, as.is=T)
        colnames(drugs) <- tabLabels[1,]
        ## is.factor()
        
        head(drugs)
        drugs[,"pharmaceutical_tx_started_days_to"] <-
            as.numeric(drugs[,"pharmaceutical_tx_started_days_to"])
        drugs[,"pharmaceutical_tx_ended_days_to"] <-
            as.numeric(drugs[,"pharmaceutical_tx_ended_days_to"])
        
        h <- hist(drugs[,"pharmaceutical_tx_ended_days_to"] -
                  drugs[,"pharmaceutical_tx_started_days_to"], 20,
                  xlab="treatment duration [d]", main=cancerID)
        ## Mode treatment duration
        typicalTreatDuration <- h$mids[which.max(h$counts)]
        abline(v=typicalTreatDuration, lty=2, lwd=2, col="red")
        abline(v=50, lty=2, lwd=2, col="blue")
        return(typicalTreatDuration)
    })

barplot(treatDurations, las=3, ylab="typical duration [d]", col="red")
abline(h=50, lwd=2, lty=2, col="blue")
dev.off()

## Typical treatment durations as a function of start time
pdf("treatStartDurations.pdf", height=3*3, width=3*4)
par(mfrow=c(3,4))

cancerID <- cancerIDs[1]
treatStartDuration <-
    sapply(cancerIDs, function(cancerID) {
        tabLabels <-
            read.table(sprintf("%s_UCSC/drugFile.txt", cancerID),
                       h=F, sep="\t", as.is=T)[1:3,]
        drugs <- read.table(sprintf("%s_UCSC/drugFile.txt", cancerID),
                            h=F, sep="\t", skip=3, as.is=T)
        colnames(drugs) <- tabLabels[1,]
        ## is.factor()
        
        head(drugs)
        drugs[,"pharmaceutical_tx_started_days_to"] <-
            as.numeric(drugs[,"pharmaceutical_tx_started_days_to"])
        drugs[,"pharmaceutical_tx_ended_days_to"] <-
            as.numeric(drugs[,"pharmaceutical_tx_ended_days_to"])

        durations <-
            drugs[,"pharmaceutical_tx_ended_days_to"] -
                drugs[,"pharmaceutical_tx_started_days_to"];
        plot(drugs[,"pharmaceutical_tx_started_days_to"],
             durations,
             xlim=c(0,500), ylim=c(0,500),
             ## xlim=range(c(0,durations), na.rm=T),
             xlab="start time [d]", 
             ylab="duration [d]", main=cancerID)

        abline(v=treatDurations[cancerID], lty=2, lwd=2, col="red")
        abline(h=treatDurations[cancerID], lty=2, lwd=2, col="red")
        ## abline(v=50, lty=2, lwd=2, col="blue")
        abline(0, 1, col="grey")
    })
dev.off()
