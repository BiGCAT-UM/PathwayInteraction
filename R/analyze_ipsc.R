# analyze_ipsc.R
# Analyze iPSC stem cell database
# and crete gene-level result (equivalent to limma)
# Used by neor_pw_interaction.R

library(dplyr)
library(DESeq2)

BASE_DIR <- "~/Desktop/BiGCAT/project_22q11ds/microarray_analysis/pathway_interaction"
DATA_DIR <- file.path(BASE_DIR, "data")
setwd("~/Desktop/BiGCAT/project_22q11ds/microarray_analysis/pathway_interaction")
# load expression data
stem_df <- read.csv(file=file.path(DATA_DIR, "GSE106589_geneCounts.csv"), 
                    row.names = 1,
                    stringsAsFactors=F)
dim(stem_df)
strsplit(colnames(stem_df), ".", fixed = TRUE)[[1]]
s <- sapply(colnames(stem_df), function(x) unlist(strsplit(x, ".", fixed=TRUE))[1])
