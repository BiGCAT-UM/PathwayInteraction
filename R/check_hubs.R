# check_hubs.R
# Check hub genes

genes_for_hub <- setdiff(V(one_subig)$name, "WP4657")
num_neighbors <- numeric()
for (one_gene in genes_for_hub){
  one_neigh <- unique(neighbors(one_subig, one_gene, "all"))$name
  len_neigh <- length(one_neigh)
  num_neighbors <- c(num_neighbors, len_neigh)
  if (length(one_neigh)>=10){
    print(paste(one_gene, "has", len_neigh, "neighbors!"))
    #print(one_neigh)
  }
}
gene_rank_df <- data.frame(gene=genes_for_hub,
                           neighbors=num_neighbors,
                           stringsAsFactors = F)
gene_rank_df <- gene_rank_df[order(gene_rank_df$neighbors, decreasing = T),]
rownames(gene_rank_df) <- NULL
