# pw_ora.R
# Pathway over-representation analysis

# sub_limma <- subset(all_limma, weight <= 0.9 & weight >= 0.899)
# # hist(abs(sub_limma$logFC))
# # hist(sub_limma$P.Value)
# summary(abs(sub_limma$logFC))
# # [1] 0.1311532
# max(sub_limma$P.Value)
# # [1] 0.05379499

# neuro_limma 
# sub_neuro_limma <- subset(neuro_limma, weight <= 0.90 & weight >= 0.899)
# summary(abs(sub_neuro_limma$logFC))
# max(sub_neuro_limma$P.Value)
# # 
# 
# # nonneuro_limma
# sub_nonneuro_limma <- subset(nonneuro_limma, weight <= 0.9 & weight >= 0.8999)
# summary(abs(sub_nonneuro_limma$logFC))
# max(sub_nonneuro_limma$P.Value)

nrow(subset(neuro_limma, P.Value<0.05 & abs(logFC)>=0.15))
nrow(subset(neuro_limma, P.Value<0.05 & logFC>=0.15))
nrow(subset(nonneuro_limma, P.Value<0.05 & abs(logFC)>=0.15))
nrow(subset(nonneuro_limma, P.Value<0.05 & logFC>=0.15))

# Direct gene analysis
library(data.table)
disgenet_df <- as.data.frame(read.table(file=file.path(DATA_DIR, "curated_gene_disease_associations.tsv"),
                          sep="\t", header=TRUE,na.strings = "", fill=FALSE, quote=""))
schizo_df <- subset(disgenet_df, tolower(diseaseName)=="schizophrenia")
autism_df <- subset(disgenet_df, tolower(diseaseName)=="autistic disorder")
psycho_df <- subset(disgenet_df, tolower(diseaseName)=="psychotic disorders")

total_neuro_genes <- union(union(schizo_df$geneId,autism_df$geneId), psycho_df$geneId)
# all_neuro_genes <- union(union(subset(schizo_df, score>=0.4)$geneId,
#                                subset(autism_df, score>=0.4)$geneId),
#                          subset(psycho_df, score>=0.4)$geneId)

# Genes exclusively differentially expressed among neuro-patients with 22q11DS
neuro_gene_df <- subset(neuro_limma, P.Value<0.05 & abs(logFC)>=0.15 & entrez_id %in% total_neuro_genes)
nonneuro_gene_df <- subset(nonneuro_limma, P.Value<0.05 & abs(logFC)>=0.15 & entrez_id %in% total_neuro_genes)

excl_genes <- setdiff(neuro_gene_df$hgnc_symbol, nonneuro_gene_df$hgnc_symbol)
common_genes <- intersect(neuro_gene_df$hgnc_symbol, nonneuro_gene_df$hgnc_symbol)
out_genes <- setdiff(nonneuro_gene_df$hgnc_symbol, neuro_gene_df$hgnc_symbol)

excl_entrez <- unique(subset(disgenet_df, geneSymbol %in% excl_genes)$geneId)
out_entrez <- unique(subset(disgenet_df, geneSymbol %in% out_genes)$geneId)


excl_wpids <- character()
excl_descs <- character()
excl_num_genes <- numeric()
excl_norm_vals <- numeric()
excl_symbols <- character()
for (one_wpid in wpids){
  one_gene_set <- subset(wpid2gene, wpid==one_wpid)$gene
  one_wpid_desc <- subset(wpid2name, wpid==one_wpid)$name
  gene_intersect <- intersect(one_gene_set, excl_entrez)
  length_intersect <- length(gene_intersect)
  if (length_intersect>=1) {
    print(paste("Pathway", one_wpid, one_wpid_desc, "has intersect:", length_intersect))
    print(paste("Normalized value", length_intersect/length(one_gene_set)))
    excl_wpids <- c(excl_wpids, one_wpid)
    excl_descs <- c(excl_descs, one_wpid_desc)
    excl_num_genes <- c(excl_num_genes, length_intersect)
    excl_norm_vals <- c(excl_norm_vals, length(intersect)/length(one_gene_set))
    excl_symbol <- paste(unique(subset(disgenet_df, geneId %in% gene_intersect)$geneSymbol), collapse=";")
    excl_symbols <- c(excl_symbols, excl_symbol)
  }
}
excl_df <- data.frame(wpid=excl_wpids,
                      descs=excl_descs,
                      num_gene=excl_num_genes,
                      norm_val=excl_norm_vals,
                      gene=excl_symbols,
                      stringsAsFactors = F)
excl_df <- excl_df[order(excl_df$num_gene, excl_df$norm_val, decreasing=T),]



excl_genes

######################################################
# over-representation analysis

# Based on below, enrichment analysis criterion..
# TODO: define limma_result
limma_result <- nonneuro_limma
s <- subset(limma_result, P.Value<0.05 & abs(logFC)>=0.15)

print(paste("Number of total DEGs:", nrow(s)))
print(paste("Number of Down-regulated genes:", nrow(subset(s, logFC<0))))
print(paste("Number of Up-regulated genes:", nrow(subset(s, logFC>=0))))

# below will be the list of all genes

# all_N and all_R are common parameters
performORA <- function(limma_df, min_positive=1, sig_genes=NULL){
  if (length(sig_genes)>0){
    sig_limma <- subset(limma_df, P.Value<0.05 & abs(logFC)>=0.15 & entrez_id %in% sig_genes)
  } else {
    sig_limma <- subset(limma_df, P.Value<0.05 & abs(logFC)>=0.15)
  }
  
  all_N <- length(unique(wpid2gene$gene))
  all_R <- length(subset(sig_limma,
                         sig_limma$entrez_id %in% wpid2gene$gene)$entrez_id)
  wpids <- unique(wpid2gene$wpid)
  count <- 0
  z_scores <- numeric()
  wp_ids <- character()
  wp_descs <- character()
  n_measures <- numeric()
  r_positives <- numeric()
  p_values <- character()
  positive_genes <- character()
  for (one_wpid in wpids){
    one_wpid_desc <- subset(wpid2name, wpid==one_wpid)$name
    one_gene_set <- subset(wpid2gene, wpid==one_wpid)$gene
    positive_r_limma <- subset(sig_limma, entrez_id %in% one_gene_set)
    positive_r <- length(unique(positive_r_limma$entrez_id))
    measured_n <- nrow(subset(limma_df, entrez_id %in% one_gene_set))
    if (positive_r>=min_positive & measured_n>0){
      print(paste(one_wpid, "is measurable"))
      count <- count + 1
      one_z_score <- (positive_r - measured_n*all_R/all_N) / 
        sqrt( measured_n*(all_R/all_N)*(1-all_R/all_N)*(1-(measured_n-1)/(all_N-1)) )
      z_scores <- c(z_scores, one_z_score)
      
      one_pval <- pnorm(-one_z_score)
      one_pval_str <- ifelse(one_pval<0.001, "<0.001", sprintf("%.3f", one_pval))
      p_values <- c(p_values, one_pval_str)
      
      wp_ids <- c(wp_ids, one_wpid)
      # Make other column data
      r_positives <- c(r_positives, positive_r)
      n_measures <- c(n_measures, measured_n)
      wp_descs <- c(wp_descs, one_wpid_desc)
      p_genes <- positive_r_limma$hgnc_symbol
      positive_genes <- c(positive_genes, 
                          paste(p_genes[1:ifelse(length(p_genes)<=10,length(p_genes),10)], collapse=";")
      )
    }
  }
  res_ora <- data.frame(wpid=wp_ids,
                        description=wp_descs,
                        positive=r_positives,
                        measured=n_measures,
                        z_score=z_scores,
                        p_values=p_values,
                        positive_genes=positive_genes,
                        stringsAsFactors = F)
  res_ora <- res_ora[order(res_ora$z_score, decreasing = TRUE),]
  res_ora$z_score <- sprintf("%.3f", res_ora$z_score)
  rownames(res_ora) <- NULL
  return(res_ora)
}

excl_neuro_ora <- performORA(limma_df=neuro_limma, min_positive=1, sig_genes=excl_entrez)
out_nonneuro_ora <- performORA(limma_df=nonneuro_limma, min_positive=1, sig_genes=out_entrez)

excl_neuro_limma <- subset(neuro_limma, entrez_id %in% excl_entrez)



## TODO: Create a network of pathways that have 2 or more genes.
# We'll need pathway description as label in the resulting Cytoscape visualization
# genes within the same pathways are connected. 
# PPI not 'yet' added. 
top_ten_wp <- excl_neuro_ora[1:10,]$wpid
big_neuro_ora <- subset(excl_neuro_ora, positive>=2 | wpid %in% top_ten_wp)
big_neuro_edges <- subset(merged_edge, type=="gene_pathway" & 
                            target %in% sapply(excl_entrez, toString) & 
                            source %in% big_neuro_ora$wpid)
big_neuro_nodes <- subset(merged_node, name %in% unique(unlist(big_neuro_edges[,c("source", "target")])))
merged_neuro_nodes <- merge(big_neuro_nodes, wpid2name, by.x="name", by.y="wpid", all.x=TRUE)
merged_neuro_nodes$label <- ifelse(merged_neuro_nodes$type=="gene_product",
                                   merged_neuro_nodes$hgnc_symbol,
                                   merged_neuro_nodes$name.y)
ora_g <- graph_from_data_frame(big_neuro_edges,
                               directed=TRUE,
                               vertices=merged_neuro_nodes)
RCy3::createNetworkFromIgraph(ora_g,
                              title = "EXCL ORA Network",
                              collection = "Filtered Network with threshold 1.4 (for Paper)")
RCy3::loadTableData(neuro_limma, data.key.column = "ensembl_id", table.key.column = "ensembl_id")

