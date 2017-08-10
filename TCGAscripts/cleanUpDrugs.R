rm(list=ls())

doFocusOnInitialTreatment <- F;

tabLabels <-
    read.table("drugFile.txt", h=F, sep="\t", as.is=T)[1:3,]
drugs <- read.table("drugFile.txt", h=F, sep="\t", skip=3, as.is=T)
colnames(drugs) <- tabLabels[1,]
## is.factor()

head(drugs)
drugs[,"pharmaceutical_tx_started_days_to"] <-
    as.numeric(drugs[,"pharmaceutical_tx_started_days_to"])
drugs[,"pharmaceutical_tx_ended_days_to"] <-
    as.numeric(drugs[,"pharmaceutical_tx_ended_days_to"])

h <- hist(drugs[,"pharmaceutical_tx_ended_days_to"] -
          drugs[,"pharmaceutical_tx_started_days_to"], 20,
          xlab="treatment duration [d]", main="")
## Mode treatment duration
typicalTreatDuration <- h$mids[which.max(h$counts)]
abline(v=typicalTreatDuration, lty=2)

## Let's read in the days to collection information, so we can look
## up the most relevant treatment outcome and drugs?
clinicalData <- read.table("clinicalData_reOrdered.tsv", sep="\t",
                           h=T, as.is=T)
clinicalData[,c("X_PATIENT", "days_to_collection")]

## Let's read in the days to collection information, so we can look
## the most relevant treatment outcome and drugs?
radLabels <-
    read.table("radiationFile.txt", h=F, sep="\t", as.is=T)[1:3,]
radiations <- read.table("radiationFile.txt", h=F, sep="\t", skip=3, as.is=T)
colnames(radiations) <- radLabels[1,]
radiations[,"radiation_therapy_started_days_to"] <-
    as.numeric(radiations[,"radiation_therapy_started_days_to"])
radiations[,"radiation_therapy_ended_days_to"] <-
    as.numeric(radiations[,"radiation_therapy_ended_days_to"])

table(radiations[,"treatment_best_response"])
table(drugs[,"treatment_best_response"])

## Let's collapse [Not Applicable], [Not Available], [Unknown] into NA
for (i in c("[Not Applicable]", "[Not Available]", "[Unknown]", "[Discrepancy]")) {
    drugs[!is.na(drugs[,"treatment_best_response"]) &
          drugs[,"treatment_best_response"] == i,
          "treatment_best_response"] <- NA
    radiations[!is.na(radiations[,"treatment_best_response"]) &
          radiations[,"treatment_best_response"] == i,
          "treatment_best_response"] <- NA
}
radiations[!is.na(radiations[,"treatment_best_response"]) &
          radiations[,"treatment_best_response"] == "Radiographic Progressive Disease",
          "treatment_best_response"] <- "Clinical Progressive Disease";
table(drugs[,"treatment_best_response"])
table(radiations[,"treatment_best_response"])

mergeDrugRadiation <- function(thisRadiation, thisDrug) {
    thisTreat <-
        thisDrug[,c("pharmaceutical_therapy_drug_name",
                    "pharmaceutical_tx_started_days_to",
                    "pharmaceutical_tx_ended_days_to",
                    "treatment_best_response")]
    colnames(thisTreat) <-
        c("therapy_name", "started_days_to",
          "ended_days_to", "treatment_best_response");
    tmp <-
        data.frame(rep("radiotherapy", nrow(thisRadiation)),
                   thisRadiation[,c("radiation_therapy_started_days_to",
                                    "radiation_therapy_ended_days_to",
                                    "treatment_best_response")])
    colnames(tmp) <- colnames(thisTreat)
    thisTreat <- rbind(thisTreat, tmp)
    return(thisTreat)
}

## Are samples typically taken before or after treatment? After!
patient <- unique(clinicalData[,"X_PATIENT"])[1]
## patient <- "TCGA-E2-A155";
## patient <- "TCGA-B6-A0IN";
patient <- "TCGA-E9-A1RB";
patient <- "TCGA-XX-A899";
pdf("drugsAndCollTimelines.pdf", height=4, width=4);
daysTo <-
    sapply(unique(clinicalData[,"X_PATIENT"]), function(patient) {
        cat(sprintf("%s\n", patient))
        daysToColl <-
            clinicalData[clinicalData[,"X_PATIENT"] == patient,
                         "days_to_collection"]
        daysToDeath <- #you die only once...
            unique(clinicalData[clinicalData[,"X_PATIENT"] == patient,
                                "days_to_death"])

        thisRadiation <-
            radiations[radiations[,"bcr_patient_barcode"] == patient,]
        
        thisDrug <-
            drugs[drugs[,"bcr_patient_barcode"] == patient,]

        thisTreat <-
            mergeDrugRadiation(thisRadiation, thisDrug);

        ## thisDrug[,"pharmaceutical_tx_started_days_to"]
        ## thisDrug[,"pharmaceutical_tx_ended_days_to"]
        timeRange <-
            range(c(thisTreat[,"started_days_to"],
                    thisTreat[,"ended_days_to"],
                    0, daysToColl, daysToDeath), na.rm=T)
        
        plot(0,0, type="n", ylim=c(0,nrow(thisTreat)), xlim=timeRange,
             main=patient, xlab="time since diagnostic [d]", ylab="")
        sapply(daysToColl, function(x) {
            abline(v=x, lty=2, col="grey") })
        abline(v=daysToDeath, col="black")
        if ( nrow(thisTreat) == 0 ) { return() }
        i <- 1
        sapply(1:nrow(thisTreat), function(i) {
            lty <- 1;
            
            if ( is.na(thisTreat[i,"started_days_to"]) ) {
                thisTreat[i,"started_days_to"] <-
                    min(timeRange);
                lty <- 2;
            }
            
            if ( is.na(thisTreat[i,"ended_days_to"]) ) {
                if ( is.na(daysToDeath) ) {
                    thisTreat[i,"ended_days_to"] <-
                        max(timeRange);
                } else {
                    thisTreat[i,"ended_days_to"] <-
                        daysToDeath;
                }
                lty <- 2;
            }
            
            col <- "grey";
            if ( is.na(thisTreat[i,"treatment_best_response"]) ) {
                col <- "grey";
            } else if ( thisTreat[i,"treatment_best_response"] ==
                       "Complete Response" ) {
                col <- "#138737";
            } else if ( thisTreat[i,"treatment_best_response"] ==
                       "Partial Response" ||
                       thisTreat[i,"treatment_best_response"] ==
                       "Stable Disease" ) {
                col <- "yellow";
            } else if ( thisTreat[i,"treatment_best_response"] ==
                       "Clinical Progressive Disease" ) {
                col <- "red";
            }

            lines(c(thisTreat[i,"started_days_to"],
                    thisTreat[i,"ended_days_to"]),
                  c(i, i), col=col, lwd=2, lty=lty)
            
            text(mean(c(thisTreat[i,"started_days_to"],
                        thisTreat[i,"ended_days_to"])),
                 c(i,i), thisTreat[i,"therapy_name"],
                 pos=1, cex=.8)
        })
        ## legend("bottomright", c("complete response", "stable", "progression"),
        ##        col=c("#138737", "yellow", "red"), lty=1, lwd=2)
        return(
            c(treat_started_days_to=min(thisTreat[i,"started_days_to"]),
              death_days_to=daysToDeath,
              collection_days_to=min(daysToColl)));
    })
dev.off();

## plot(
##     unlist(sapply(daysTo,
##                   function(x) {
##                       x["treat_started_days_to"] })),
##     unlist(sapply(daysTo,
##                   function(x) { x["collection_days_to"] })),
##     xlab="days to treatment start", ylab="days to collection")
## abline(0,1)

drugDic <- read.csv("../drugDictionnary.csv", as.is=T);
head(drugDic)
sel <- (table(drugDic[,1]) > 1)
## Here are drug cocktails, with the number of drugs in them:
table(drugDic[,1])[sel]
                  
drugLists <-
    strsplit(
        tolower(
            gsub("/", "+",
                 gsub(" and ", "+",
                      gsub(", ", "+",
                           gsub(", and ", ", ",
                                gsub(" \\(.*\\)$", "",
                                     drugs[,"pharmaceutical_therapy_drug_name"])
                                )
                           )
                      ),
                 fixed=T)
            ),
        split="[ ]*\\+[ ]*", fix=F)

## Let's extract a list of all drugs given:
splitDrugs <- unique( unlist( drugLists ) )

stdDrugs <-
    sapply(splitDrugs, function(x) {
            sel <- drugDic[,1] == x
            if ( sum(sel) == 0 ) {
                cat(sprintf("Could not map '%s'\n", x))
                return()
            } else {
                return(drugDic[sel, 2])
            }
        })

## Which drugs couldn't be mapped?
sort(table(splitDrugs[sapply(stdDrugs, length) == 0]))

## We now build a table, where each row is a patient, and columns
## register the earliest / latest date of the treatment and a boolean
## vector of all drugs given

## if ( doFocusOnInitialTreatment ) {
##     isInitTreat <- !is.na(drugs[,"pharmaceutical_tx_started_days_to"]) &
##         drugs[,"pharmaceutical_tx_started_days_to"] < typicalTreatDuration;
##     drugsInit <- drugs[isInitTreat,]
## } else {
##     drugsInit <- drugs;
## }

## We iterate through all patients
drugLabels <- unique(drugDic[,2]);
drugLabels <- c(drugLabels, "radiotherapy");

patient <- unique(drugs[,"bcr_patient_barcode"])[2]

patientTreatments <- 
    sapply(unique(c(drugs[,"bcr_patient_barcode"],
                    radiations[,"bcr_patient_barcode"])),
           function(patient) {
               ## Fetch drugs and best response
               thisPatient <- drugs[,"bcr_patient_barcode"] == patient
               if ( doFocusOnInitialTreatment ) {
                   thisPatient <- thisPatient &
                       !is.na(drugs[,"pharmaceutical_tx_started_days_to"]) & 
                           drugs[,"pharmaceutical_tx_started_days_to"] <
                               typicalTreatDuration
               }
               thisPatient <- which(thisPatient);
               idx <- thisPatient[1]
               ## We iterate over all treatments given to this patient
               givenDrugs <-
                   unlist(
                       sapply(thisPatient, function(idx) {
                           ## A treatment may include several drugs, over which
                           ## we iterate
                           unlist(
                               sapply(drugLists[[idx]], function(x) {
                                   sel <- drugDic[,1] == x
                                   if ( sum(sel) == 0 ) {
                                       cat(sprintf("Could not map '%s'\n", x))
                                       return()
                                   } else {
                                       return(drugDic[sel, 2])
                                   }
                               })
                               )
                       })
                       )
               tVec <- rep(F, length(drugLabels))
               names(tVec) <- drugLabels;
               for ( x in givenDrugs ) {
                   tVec[x] <- T;
               }

               radioSel <- radiations[,"bcr_patient_barcode"] == patient;
               if ( doFocusOnInitialTreatment ) {
                   radioSel <- radioSel &
                       !is.na(radiations[,"radiation_therapy_started_days_to"]) &
                           radiations[,"radiation_therapy_started_days_to"] <
                               typicalTreatDuration
               }
               tVec["radiotherapy"] <- any(radioSel)
               radioResponse <- radiations[radioSel,"treatment_best_response"]
               treatResp <-
                   unique(c(drugs[thisPatient,"treatment_best_response"],
                            radiations[radioSel,"treatment_best_response"]))
               if ( any(is.na(treatResp)) ) {
                   treatResp <- treatResp[!is.na(treatResp)]
               }
               return(list(treatment_best_response=treatResp,
                           treatment_vector=tVec))
           })

patient <- colnames(patientTreatments)[2];
## patient <- "TCGA-A1-A0SG"
resMat <-
    sapply(colnames(patientTreatments), function(patient) {
        as.numeric(patientTreatments[["treatment_vector",patient]])
        ## cat(sprintf("%s\n", patient))
    })
rownames(resMat) <- drugLabels; 

head(resMat)
sort(apply(resMat, 1, sum)) #most common therapies
## as.logical(resMat)

##################################################
## Summarize by drug target

drugTargets <- read.csv("../drugTargets.csv", as.is=T)
drugTargets[is.na(drugTargets)] <- 0

d <- "radiotherapy"
alnDrugTargets <-
    sapply(rownames(resMat), function(d) {
        sel <- drugTargets[,1] == d;
        if ( sum(sel) != 1 ) { cat(sprintf("Trouble with %s\n", d)) }
        return(as.numeric(drugTargets[sel,c(-1, -2)]))
    })

rownames(alnDrugTargets) <-
    gsub("^", "target.", colnames(drugTargets)[c(-1,-2)])
resMatS <- alnDrugTargets %*% resMat
resMatS[resMatS > 1] <- 1

resMatF <- rbind(resMat, resMatS)

## Is looks as if the best therapy outcome is reported for each
## individual treatment but is actually the same within each
## patient. Let's confirm by comparing the radiotherapy best treatment
## outcomes to the drug best treatment outcomes when patients
## underwent both.

## patientId <- intersect(colnames(patientDrugs),
##                        radiations[,"bcr_patient_barcode"])[1]
## sapply(intersect(colnames(patientDrugs),
##                  radiations[,"bcr_patient_barcode"]),
patientId <- colnames(patientTreatments)[2]
patientBestResp <-
    sapply(colnames(patientTreatments),
           function(patientId) {
               allResps <- unlist(patientTreatments["treatment_best_response",patientId])
               allResps <- unique(allResps[!is.na(allResps)])
               if ( length(allResps) == 0 ) { allResps <- NA; }
               if ( length(allResps) > 1 ) {
                   bestOrdering <-
                       c("Complete Response",
                         "Partial Response",
                         "Stable Disease",
                         "Clinical Progressive Disease");
                   allResps <-
                       names(
                           sort(
                               sapply(allResps,
                                      function(x) {
                                          which(bestOrdering == x) })
                               )
                           )[1]
               }
               return(allResps)
           })
    
drugTab <-
    data.frame(colnames(resMatF),
               t(resMatF));
write.csv(drugTab, file="drugTab.csv", row.names=F)

treatTab <-
    merge(
        data.frame(names(patientBestResp),
                   patientBestResp),
        data.frame(colnames(resMatF),
                   t(resMatF)),
        by=1, all.x=T)
treatTab[,2] <- as.character(treatTab[,2])
## write.csv(treatTab, file="treatTab.csv", row.names=F)
write.table(treatTab, file="treatTab.tsv", row.names=F, sep="\t", quote=F)

samplesOrdered <-
    read.table("patientIDs.list", as.is=T, h=F)[,1]
sId <- samplesOrdered[2]
sId <- "TCGA-EL-A3ZS-11"
treatReordered <-
    t(sapply(samplesOrdered,
             function(sId) {
                 patId <- gsub("-[0-9]*$", "", sId);
                 sel <- treatTab[,1] == patId;
                 if ( sum(sel) == 1 ) {
                     return(unlist(treatTab[sel,-1]))
                 } else if ( sum(sel) == 0 ) {
                     return(rep(NA, ncol(treatTab) - 1))
                 } else {
                     stop("This shouldn't happen\n")
                 }
             }))

## sum(colnames(treatReordered) != colnames(treatTab)[-1])
## colnames(treatReordered) <- colnames(treatTab)[-1]
head(treatReordered)

## treatReordered <- as.data.frame(treatReordered)
## head(treatReordered)
## treatReordered[1:80,1]

## should be clinical progression in BRCA
## patientBestResp["TCGA-A8-A08O"]

write.table(treatReordered, file="treatments_reOrdered.tab",
            sep="\t", quote=F)

discreteFeats <-
    read.table("discreteClinicalData_reOrdered_withoutTreatment.tsv",
               h=T, as.is=T, sep="\t")

colnames(treatReordered) <- sprintf("treat.%s", colnames(treatReordered))
discTreat <- cbind(discreteFeats, treatReordered)
isNormal <- discTreat[,"sample_type"] == "Solid Tissue Normal";

## ncol(discTreat)
## ncol(discreteFeats)+seq(1,ncol(treatReordered))
## colnames(discTreat)[ncol(discreteFeats)+seq(1,ncol(treatReordered))]

discTreat[isNormal,ncol(discreteFeats)+seq(1,ncol(treatReordered))] <- NA;

write.table(discTreat,
            file="discreteClinicalData_reOrdered_withTreatment.tsv",
            sep="\t", quote=F, row.names=F)
