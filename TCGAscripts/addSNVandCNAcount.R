rm(list=ls())

library(tidyverse)
contDat <- read_tsv("continuousClinicalData_reOrdered.tsv")
isMetabric <- length(grep("metabric", getwd())) == 1

##################################################
## SNVs and CNAs

if ( isMetabric ) {
    snvDat <- read_tsv("mutMatrix_reOrdered.tsv", col_names=F) %>% select(-X1)
} else {
    snvDat <- read_tsv("mutMatrix_reOrdered_booleanized.tsv") %>% select(-sampleID)
    snvDat <- snvDat[,setdiff(1:ncol(snvDat), grep("=NaN", colnames(snvDat)))];
    ## colnames(snvDat)
}

snvCount <- apply(snvDat, 1, sum)

cnaDat <- read_tsv("copMatrix_reOrdered_booleanized.tsv") %>% select(-sampleID)
cnaDat <- cnaDat[,setdiff(1:ncol(cnaDat), grep("=NaN$", colnames(cnaDat)))];
cnaDat <- cnaDat[,setdiff(1:ncol(cnaDat), grep("=0$", colnames(cnaDat)))];
## colnames(cnaDat)
cnaCount <- apply(cnaDat, 1, sum)

contDat <- 
    contDat %>% mutate(number_of_SNVs=snvCount,
                       number_of_SNVs_frac=percent_rank(snvCount),
                       number_of_CNAs=cnaCount,
                       number_of_CNAs_frac=percent_rank(cnaCount))

##################################################
## Import estimate scores

library(stringr)
cType <- str_replace(getwd(), "_UCSC", "") %>% str_replace("^.*\\/", "")
estDat <-
    read_csv("../AranNatComm2015.csv") %>%
    select(`Sample ID`, `Cancer type`, ESTIMATE) %>%
    mutate(`Sample ID`=str_replace(`Sample ID`, "[A-Z]$", "")) %>%
    filter(`Cancer type`==cType) %>% 
    group_by(`Sample ID`) %>% summarize(ESTIMATE=mean(ESTIMATE, na.rm=T))

## Put ESTIMATE scores in the right order
patOrder <- read_tsv("patientIDs.list", col_names=c("Sample ID"))
estDat <- left_join(patOrder, estDat) %>%
    mutate(ESTIMATE=ifelse(is.na(ESTIMATE), NaN, ESTIMATE))

contDat <- 
    contDat %>% bind_cols( select(estDat, ESTIMATE) )
View(contDat)

##################################################

write_tsv(contDat,
          path="continuousClinicalData_reOrdered.tsv");
