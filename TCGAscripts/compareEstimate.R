rm(list=ls())
library(tidyverse)

cancers <-
    read_tsv("TCGA_frac_nArchs.tab",
             col_names=c("cancer", "fracGenes", "nArch"))

cType <- "BLCA"
estDat <-
    map(cancers %>% select(cancer) %>% unlist(), function(cType) {
        read_csv(sprintf("%s_UCSC/clinicalEnrichment_continuous_All.csv", cType)) %>%
            filter(`Feature Name` == "ESTIMATE") %>%
            mutate(cancer=cType)
    }) %>% bind_rows() %>%
    select(-`Feature Name`, -`Is first bin maximal?`) %>%
    mutate(arch=sprintf("%s %d", cancer, `archetype #`)) %>%
    mutate(`Mean Difference`=100 * `Mean Difference`) %>%
    mutate(signif=`Significant after Benjamini-Hochberg correction?` == 1)

ggplot(estDat) +
    geom_tile(aes(x=cancer, y=`archetype #`,
                  fill=`Mean Difference`)) +
    scale_fill_gradientn(colors=terrain.colors(256))

ggplot(estDat) +
    geom_col(aes(x=`archetype #`, y=`Mean Difference`, fill=signif)) +
    facet_wrap(~cancer, scales="free_x") +
    ylab("Mean difference in purity (%)") +
    theme(legend.position="none")
ggsave("compareEstimate.pdf", height=5, width=5)
