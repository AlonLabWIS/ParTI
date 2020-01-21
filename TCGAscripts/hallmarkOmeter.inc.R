## arcCols <- c("#4078fa", "#ffff00", "#e379e3", "#ca5e5e", "#28c928")
arcCols <- c("lipogenesis"="#ffff00",
	     "biomass & energy"="#019e59", 
	     "cell division"="#007eb1", 
	     "immune interaction"="#d02690", 
	     "invasion & signaling"="#111111")

## toShow <-
##     list(
##         ## "Replicative immortality"=c(
##         ##     ## "REACTOME PACKAGING OF TELOMERE ENDS",
##         ##     "REACTOME EXTENSION OF TELOMERES"),
##         "Resisting cell death"=c(## "REACTOME REGULATION OF APOPTOSIS",
##             "ANTI APOPTOSIS",
##             "REACTOME EXTENSION OF TELOMERES"
##             ## "REACTOME APOPTOTIC EXECUTION PHASE"
##         ),
##         "Evading growth suppressors"=c(
##             ## "KEGG DNA REPLICATION",
##             "REACTOME DNA REPLICATION",
##             "MITOSIS"## ,
##             ## "MITOTIC CELL CYCLE"
##         ),
##         "Interaction with immune system"=c(
##             "INFLAMMATORY RESPONSE",
##             ## "BIOCARTA NKT PATHWAY",
##             ## "BIOCARTA INFLAM PATHWAY",
##             ## "HUMORAL IMMUNE RESPONSE",
##             "KEGG ALLOGRAFT REJECTION"),
##         "Energetics"=c(
##             "REACTOME GLYCOLYSIS",
##             "REACTOME RESPIRATORY ELECTRON TRANSPORT",
##             "KEGG RIBOSOME",
##             "KEGG PROTEASOME",
##             "MITOCHONDRIAL RIBOSOME",
## 	    "KEGG PENTOSE PHOSPHATE PATHWAY"
## 	),
##         "DNA repair"=c(
##             "REACTOME DOUBLE STRAND BREAK REPAIR",
##             "KEGG NUCLEOTIDE EXCISION REPAIR"
##             ## "REACTOME DNA REPAIR",
##             ## "PID FANCONI PATHWAY",
##             ## "KEGG HOMOLOGOUS RECOMBINATION"
##         ),
##         "Angiogenesis"=c(## "REGULATION OF ANGIOGENESIS",
##             ## "TRANSFORMING GROWTH FACTOR BETA RECEPTOR SIGNALING PATHWAY",
##             ## "PID VEGFR1 PATHWAY",
##             ## "BIOCARTA VEGF PATHWAY",
##             ## "KEGG VEGF SIGNALING PATHWAY",
##             ## "BLOOD COAGULATION",
##             "ANGIOGENESIS"),
##         "Signaling"=c(
##             "PID WNT SIGNALING PATHWAY",
##             "REACTOME PI3K CASCADE",
##             "PID RAS PATHWAY",
##             ## "RAS GTPASE BINDING",            
##             ## "RAS GTPASE ACTIVATOR ACTIVITY",
##             ## "RAS PROTEIN SIGNAL TRANSDUCTION",
##             "BIOCARTA MTOR PATHWAY",
##             "TRANSFORMING GROWTH FACTOR BETA RECEPTOR SIGNALING PATHWAY"
##         ),
##         "Peroxisome"=c(
##             "REACTOME PEROXISOMAL LIPID METABOLISM",
##             "KEGG GLYCOSYLPHOSPHATIDYLINOSITOL GPI ANCHOR BIOSYNTHESIS",
##             "KEGG ASCORBATE AND ALDARATE METABOLISM",
##             "KEGG PENTOSE AND GLUCURONATE INTERCONVERSIONS"
##             ## "KEGG NITROGEN METABOLISM",
##             ## "GLUTAMATE RECEPTOR ACTIVITY",
##         ),
##         "Invasion & metastasis"=c(
##             "REACTOME DEGRADATION OF THE EXTRACELLULAR MATRIX",
##             "CELL MIGRATION"
##         ),
##         "Micro-environment"=c(
##             ## "AXON GUIDANCE",
##             "TRANSMISSION OF NERVE IMPULSE",
##             ## "REACTOME BILE ACID AND BILE SALT METABOLISM",
##             "COLLAGEN"
## 	    ## "RESPONSE TO HYPOXIA",
## 	    ## "PID HIF1 TFPATHWAY" #commented Sept 14 to get 4 super-arch
## 	    ## "REACTOME REGULATION OF HYPOXIA INDUCIBLE FACTOR HIF BY OXYGEN"
##         )
##     )

toShow <-
    list(
        ## "Replicative immortality"=c(
        ##     ## "REACTOME PACKAGING OF TELOMERE ENDS",
        ##     "REACTOME EXTENSION OF TELOMERES"),
        "Resisting cell death"=c(## "REACTOME REGULATION OF APOPTOSIS",
            "ANTI APOPTOSIS",
            "REACTOME EXTENSION OF TELOMERES"
            ## "REACTOME APOPTOTIC EXECUTION PHASE"
        ),
        "Evading growth suppressors"=c(
            ## "KEGG DNA REPLICATION",
            "REACTOME DNA REPLICATION",
            "M PHASE OF MITOTIC CELL CYCLE",
            "KEGG CELL CYCLE"
            ## "MITOSIS"## ,
            ## "MITOTIC CELL CYCLE"
        ),
        "Interaction with immune system"=c(
            ## "INFLAMMATORY RESPONSE",
            ## "BIOCARTA NKT PATHWAY",
            "BIOCARTA INFLAM PATHWAY",
            ## "HUMORAL IMMUNE RESPONSE",
            "BIOCARTA TCYTOTOXIC PATHWAY",
            "REACTOME INTERFERON GAMMA SIGNALING",
            "REACTOME PD1 SIGNALING",
            "BIOCARTA CTLA4 PATHWAY",
            "KEGG ALLOGRAFT REJECTION",
            "BIOCARTA LECTIN PATHWAY"),
        "Energetics"=c(
            "REACTOME GLYCOLYSIS",
            "REACTOME RESPIRATORY ELECTRON TRANSPORT",
            "KEGG RIBOSOME",
            "MITOCHONDRIAL RIBOSOME",
            "KEGG PROTEASOME"
	    ## "KEGG PENTOSE PHOSPHATE PATHWAY"
	),
        "DNA repair"=c(
            "REACTOME DOUBLE STRAND BREAK REPAIR",
            "KEGG NUCLEOTIDE EXCISION REPAIR"
            ## "REACTOME DNA REPAIR",
            ## "PID FANCONI PATHWAY",
            ## "KEGG HOMOLOGOUS RECOMBINATION"
        ),
        "Angiogenesis"=c(## "REGULATION OF ANGIOGENESIS",
            ## "TRANSFORMING GROWTH FACTOR BETA RECEPTOR SIGNALING PATHWAY",
            ## "PID VEGFR1 PATHWAY",
            ## "BIOCARTA VEGF PATHWAY",
            ## "KEGG VEGF SIGNALING PATHWAY",
            ## "BLOOD COAGULATION",
            "PID VEGF VEGFR PATHWAY",
            "POSITIVE REGULATION OF ANGIOGENESIS"),
        "Signaling"=c(
            ## "PID WNT SIGNALING PATHWAY",
            ## "REACTOME PI3K CASCADE",
            "PID HEDGEHOG 2PATHWAY",
            "INSULIN LIKE GROWTH FACTOR RECEPTOR BINDING",
            "REACTOME FGFR LIGAND BINDING AND ACTIVATION",
            "PID WNT SIGNALING PATHWAY",
            "KEGG TGF BETA SIGNALING PATHWAY",
            "PID RAS PATHWAY",
            ## "RAS GTPASE BINDING",            
            ## "RAS GTPASE ACTIVATOR ACTIVITY",
            ## "RAS PROTEIN SIGNAL TRANSDUCTION",
            "BIOCARTA MTOR PATHWAY"
            ## "TRANSFORMING GROWTH FACTOR BETA RECEPTOR SIGNALING PATHWAY"
        ),
        "Peroxisome"=c(
            "REACTOME PEROXISOMAL LIPID METABOLISM",
            "KEGG GLYCOSYLPHOSPHATIDYLINOSITOL GPI ANCHOR BIOSYNTHESIS",
            "KEGG ASCORBATE AND ALDARATE METABOLISM"
            ## "KEGG PENTOSE AND GLUCURONATE INTERCONVERSIONS"
            ## "KEGG NITROGEN METABOLISM",
            ## "GLUTAMATE RECEPTOR ACTIVITY",
        ),
        "Invasion & metastasis"=c(
            "REACTOME DEGRADATION OF THE EXTRACELLULAR MATRIX",
            "CELL MIGRATION"
        ),
        "Micro-environment"=c(
            ## "AXON GUIDANCE",
            "TRANSMISSION OF NERVE IMPULSE",
            "REGULATION OF NEUROGENESIS",
            ## "REACTOME BILE ACID AND BILE SALT METABOLISM",
            "COLLAGEN",
            "REACTOME SMOOTH MUSCLE CONTRACTION",
	    "RESPONSE TO HYPOXIA"
	    ## "PID HIF1 TFPATHWAY" #commented Sept 14 to get 4 super-arch
	    ## "REACTOME REGULATION OF HYPOXIA INDUCIBLE FACTOR HIF BY OXYGEN"
        )
    )

fmtMSigDBnames <- function(MSigDBnames, maxWidth=35) {
    MSigDBnames %>% 
        str_replace("^REACTOME ", "") %>%
        str_replace("^KEGG ", "") %>%
        str_replace("^PID ", "") %>%
        str_replace("^BIOCARTA ", "") %>%
        str_replace("THE ", " ") %>%
        str_to_title() %>%
        str_trunc(width=maxWidth)
}

lipidShow <-
    list(
        "Fatty acid synthesis"=c(
	    "SREBF1", #TF controlling FASN expression
            "ACLY", #cytosolic ATP cytrate lyase
            "ACACA", #Acetyl-CoA carboxylase
            "FASN", #fatty acid synthase
            "SCD" #stearoyl-CoA desaturase
        ),
	"anaerobic glycosis"=c(
	    "LDHA", #lactate dehydrogenase
	    "LDHB" #lactate dehydrogenase
	),
	"NADPH generation"=c(
	    "ME1", #malic enzyme
	    "ME2" #malic enzyme
	),
        "beta-oxidation"=c(
            "ACOX1",
            ## "ACOX3",
            "CPT1B"
        ),
        "key peroxisome proteins"=c(
            "CAT",
            "PEX19"
        ),
        "other"=c(
            "GPX4"
        )
    )

##################################################

toShowFig <-
    list(
        "Resisting cell death"=c(
            "ANTI APOPTOSIS"
        ),
        "Sustaining proliferative signaling"=c(
            "REACTOME PI3K CASCADE",
            "PID RAS PATHWAY"
        ),
        "Evading growth suppressors"=c(
            "REACTOME DNA REPLICATION",
            "MITOSIS"
        ),
        "Activating invasion & metastasis"=c(
            "REACTOME DEGRADATION OF THE EXTRACELLULAR MATRIX",
            "CELL MIGRATION"
            ## "TRANSFORMING GROWTH FACTOR BETA RECEPTOR SIGNALING PATHWAY"
        ),
        "Enabling replicative immortability"=c(
            "REACTOME EXTENSION OF TELOMERES"
        ),
        "Angiogenesis"=c(
            ## "ANGIOGENESIS"
            "PID VEGF VEGFR PATHWAY"
        ),
        "Deregulating cellular energetics"=c(
            "REACTOME GLYCOLYSIS",
            "REACTOME RESPIRATORY ELECTRON TRANSPORT",
            ## "KEGG PENTOSE PHOSPHATE PATHWAY",
            "REACTOME PEROXISOMAL LIPID METABOLISM"
        ),
        "Avoiding immune destruction"=c(
            ## "KEGG ALLOGRAFT REJECTION"
            "BIOCARTA CTLA4 PATHWAY"
        ),
        "Tumor-promoting inflammation"=c(
            ## "INFLAMMATORY RESPONSE"
            "BIOCARTA INFLAM PATHWAY"
        )
    )
