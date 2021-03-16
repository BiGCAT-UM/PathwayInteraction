# get_neuro_genes.R
# load and get neuro genes

library(data.table)
disgenet_df <- as.data.frame(read.table(file=file.path(DATA_DIR, "curated_gene_disease_associations.tsv"),
                                        sep="\t", header=TRUE,na.strings = "", fill=FALSE, quote=""))
schizo_df <- subset(disgenet_df, tolower(diseaseName)=="schizophrenia")
autism_df <- subset(disgenet_df, tolower(diseaseName)=="autistic disorder")
psycho_df <- subset(disgenet_df, tolower(diseaseName)=="psychotic disorders")

total_neuro_genes <- union(union(schizo_df$geneId,autism_df$geneId), psycho_df$geneId)
