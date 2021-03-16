# get_ora_network.R
# Calculte neuropsychiatric subnetwork based on ORA

# Get the list of such genes (DEG & )
neuro_gene_df <- subset(neuro_limma, P.Value<0.05 & abs(logFC)>=0.15 & entrez_id %in% total_neuro_genes)
nonneuro_gene_df <- subset(nonneuro_limma, P.Value<0.05 & abs(logFC)>=0.15 & entrez_id %in% total_neuro_genes)

# hgnc IDs
excl_ora_genes <- setdiff(neuro_gene_df$hgnc_symbol, nonneuro_gene_df$hgnc_symbol)
common_ora_genes <- intersect(neuro_gene_df$hgnc_symbol, nonneuro_gene_df$hgnc_symbol)
out_ora_genes <- setdiff(nonneuro_gene_df$hgnc_symbol, neuro_gene_df$hgnc_symbol)

excl_ora_entrez <- unique(subset(disgenet_df, geneSymbol %in% excl_genes)$geneId)
out_ora_entrez <- unique(subset(disgenet_df, geneSymbol %in% out_genes)$geneId)

# Now, using the list of excl_ora_entrez, find pathways that include such genes
excl_ora_entrez

wpid_with_excl_entrez_list <- list()
for (one_wpid in unique(wpid2gene$wpid)){
  sub_wpid2gene <- subset(wpid2gene, wpid==one_wpid)
  wpid_with_excl_entrez <- subset(sub_wpid2gene, gene %in% excl_ora_entrez)
  if (nrow(wpid_with_excl_entrez)>=2){
    wpid_with_excl_entrez_list[[(length(wpid_with_excl_entrez_list)+1)]] <- wpid_with_excl_entrez
    #print(paste(one_wpid, "has two or more genes~"))
  }
}

excl_neuro_entrez_df <- bind_rows(wpid_with_excl_entrez_list)[,c("gene", "wpid")]

excl_neuro_nodes <- unique(unlist(excl_neuro_entrez_df))
excl_neuro_sub_df <- subset(neuro_merged_node, name %in% excl_neuro_nodes)
excl_neuro_node_df_merged <- merge(excl_neuro_sub_df, wpid2name, by.x="name", by.y="wpid", all.x=TRUE)
excl_neuro_node_df_merged$label <- c(excl_neuro_node_df_merged$hgnc_symbol[!is.na(excl_neuro_node_df_merged$hgnc_symbol)],
                                     excl_neuro_node_df_merged$name.y[!is.na(excl_neuro_node_df_merged$name.y)])
excl_neuro_node_df_full <- merge(excl_neuro_node_df_merged, neuro_limma[,c("ensembl_id", "logFC")], by="ensembl_id", all.x=TRUE)

excl_neuro_node_df <- excl_neuro_node_df_full[,c("label", "ensembl_id", "logFC")]
colnames(excl_neuro_node_df) <- c("name", "ensembl_id", "logFC")


excl_neuro_wpid_name <- merge(excl_neuro_entrez_df, wpid2name, by="wpid")
excl_neuro_edge_df <- merge(excl_neuro_wpid_name, neuro_merged_node, by.x="gene", by.y="name")[,c("hgnc_symbol", "name")]
rownames(excl_neuro_edge_df) <- NULL
colnames(excl_neuro_edge_df) <- c("source", "target")


excl_ora_g <- graph_from_data_frame(excl_neuro_edge_df,
                                    directed=FALSE,
                                    vertices=excl_neuro_node_df)

RCy3::cytoscapePing()
RCy3::createNetworkFromIgraph(excl_ora_g,
                              title = "Neuro ORA Network",
                              collection = "Filtered Network with threshold 1.4 (for Paper), 14 Feb 2021")
toggleGraphicsDetails()



## ----------------------------------

# Summarize pathways with gene-inclusion list (for supplementary materials)

wpids <- character()
pws <- character()
num_genes <- numeric()
name_genes <- character()
for (one_wpid in unique(wpid2gene$wpid)){
  sub_wpid2gene <- subset(wpid2gene, wpid==one_wpid)
  wpid_with_excl_entrez <- subset(sub_wpid2gene, gene %in% excl_ora_entrez)
  if (nrow(wpid_with_excl_entrez)>=1){
    wpids <- c(wpids, one_wpid)
    pws <- c(pws, subset(wpid2name, wpid==one_wpid)$name)
    num_genes <- c(num_genes, nrow(wpid_with_excl_entrez))
    hgnc_genes <- subset(neuro_merged_node, name %in% wpid_with_excl_entrez$gene)$hgnc_symbol
    name_genes <- c(name_genes, paste(hgnc_genes, collapse=";"))
    #print(paste(one_wpid, "has two or more genes~"))
  }
}
pw_containing_excl_genes <- data.frame(id=wpids,
                                       pathway=pws,
                                       number_of_genes=num_genes,
                                       name_of_genes=name_genes,
                                       stringsAsFactors = F)
res_excl_df <- pw_containing_excl_genes[order(pw_containing_excl_genes$number_of_genes, decreasing = T),]

library(xlsx)
write.xlsx(res_excl_df,
          file.path(SUPPL_DIR, "pathways_containing_exclusive_neuropsychiatric_genes.xlsx"),
          row.names = F)
