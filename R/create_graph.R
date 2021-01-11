# create_graph.R
# create graph, based on different limma_result


# Gene-pathway edge
sub_wpid <- subset(wpid2gene, gene %in% limma_result$entrez_id)
gp_edge <- data.frame(source=c(sub_wpid$wpid, sub_wpid$gene),
                      target=c(sub_wpid$gene, sub_wpid$wpid),
                      type="gene_pathway",
                      stringsAsFactors = F)
# PPI
#######################
# PPI interaction (900)
#######################
ppi <- read.table(
  file=file.path(DATA_DIR, "ppi_homo_sapiens_filtered900.txt"),
  sep=" ", header=TRUE, stringsAsFactors = F)
ppi$protein_a <- sapply(ppi$protein1, function(x) unlist(strsplit(x, "[.]"))[[2]])
ppi$protein_b <- sapply(ppi$protein2, function(x) unlist(strsplit(x, "[.]"))[[2]])
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

protein_a_ref <- getBM(attributes = c("ensembl_peptide_id", "entrezgene_id", "hgnc_symbol"),
                       filters = "ensembl_peptide_id",
                       values = unique(ppi$protein_a),
                       bmHeader = T,
                       mart = mart)
colnames(protein_a_ref) <- c("protein_a", "entrez_id_a", "hgnc_a")
protein_b_ref <- getBM(attributes = c("ensembl_peptide_id", "entrezgene_id", "hgnc_symbol"),
                       filters = "ensembl_peptide_id",
                       values = unique(ppi$protein_b),
                       bmHeader = T,
                       mart = mart)
colnames(protein_b_ref) <- c("protein_b", "entrez_id_b", "hgnc_b")

ppi_a <- merge(ppi[,c("protein_a", "protein_b")],
               protein_a_ref,
               by="protein_a")
ppi_ab <- merge(ppi_a,
                protein_b_ref,
                by="protein_b")
ppi_comp <- ppi_ab[complete.cases(ppi_ab),]
ppi_with_inv <- data.frame(source=c(as.character(ppi_comp$entrez_id_a),
                                    as.character(ppi_comp$entrez_id_b)),
                           target=c(as.character(ppi_comp$entrez_id_b),
                                    as.character(ppi_comp$entrez_id_a)),
                           type="protein_protein_interaction",
                           stringsAsFactors = F)

ppi_edge <- subset(ppi_with_inv, source %in% limma_result$entrez_id & 
                     target %in% limma_result$entrez_id)
ppi_edge <- ppi_edge[!duplicated(ppi_edge),]
rownames(ppi_edge) <- NULL
###################################################
# Combine gp and ppi networks
full_edge <- rbind(gp_edge, ppi_edge)
merged_edge <- merge(full_edge,
                     limma_result[,c("entrez_id", "weight")],
                     by.x = "target",
                     by.y = "entrez_id",
                     all.x = TRUE)
merged_edge <- merged_edge[!duplicated(merged_edge),]
merged_edge[is.na(merged_edge)] <- 0
merged_edge <- merged_edge[,c("source", "target", "type", "weight")]
full_node <- data.frame(name=unique(merged_edge$source),
                        stringsAsFactors = F)
merged_node <- merge(full_node,
                     limma_result[,c("entrez_id", "ensembl_id", "hgnc_symbol", "weight")],
                     by.x = "name",
                     by.y = "entrez_id",
                     all.x = TRUE)
merged_node <- merged_node[!duplicated(merged_node$name),]
merged_node$type <- ifelse(is.na(merged_node$hgnc_symbol),
                           "pathway",
                           "gene_product")
# merged_edge and merged_node will be used to create a graph
g <- graph_from_data_frame(merged_edge, directed=TRUE, vertices=merged_node)