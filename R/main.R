# main.R
# TODO 1: Run limma for autisom+schizophrenia, and find interaction network, and find drug targets network

#clean workspace
rm(list=ls())

library(biomaRt)
library(clusterProfiler)
# library(ComplexHeatmap)
library(DOSE)
library(dplyr)
library(httr)
library(igraph)
library(illuminaHumanv4.db)
library(limma)
library(PKNCA)
# library(propagate)
library(rstudioapi)
library(rWikiPathways)
library(RCy3)
library(SPARQL)
library(stringr)
library(tidyr)
BASE_DIR <- "~/Desktop/BiGCAT/project_22q11ds/microarray_analysis/pathway_interaction"
DATA_DIR <- file.path(BASE_DIR, "data")
RESULT_DIR <- file.path(BASE_DIR, "result")
setwd("~/Desktop/BiGCAT/project_22q11ds/microarray_analysis/pathway_interaction")
source("functions.R")

#######################
# expression data load (geometric mean)
#######################
probe_df <- read.csv(file=file.path(DATA_DIR, "GSE59216_normalized_ageOut.csv"), 
                     row.names = 1, stringsAsFactors=F)

# probe to Ensembl
probe_to_ensembl <- data.frame(Gene=unlist(mget(x = rownames(probe_df), envir = illuminaHumanv4ENSEMBL)))
probe_to_entrez <- data.frame(Gene=unlist(mget(x = rownames(probe_df), envir = illuminaHumanv4ENTREZID)))
# gene expression data, calculated by geometric mean
# load(file=file.path(DATA_DIR, "geometric_mean_probes.rda"))
# observation info, needed for limma (patients vs. controls)
obs_info <- read.csv(file=file.path(DATA_DIR, "obs_info_processed.csv"), 
                     stringsAsFactors=F)

# Choose healthy control + autism/schizophrenia patients
neuro_obs_info <- subset(obs_info, (is_patient=="Yes" & (psychotic_disorder=="Yes" | autism_spectrum_disorder=="Yes")) | 
                           is_patient=="No")
nonneuro_obs_info <- subset(obs_info, (is_patient=="Yes" & (psychotic_disorder=="No" & autism_spectrum_disorder=="No")) |
                              is_patient=="No")
rownames(neuro_obs_info) <- NULL
rownames(nonneuro_obs_info) <- NULL
neuro_exp_df <- probe_df[,neuro_obs_info$patient_id]
nonneuro_exp_df <- probe_df[,nonneuro_obs_info$patient_id]
all_exp_df <- probe_df[,obs_info$patient_id]
# Limma result using limma package. Transformed T-statistic (weight) is added. 
neuro_limma <- runLimma(df=neuro_exp_df, df_info=neuro_obs_info, dep_var="status")
nonneuro_limma <- runLimma(df=nonneuro_exp_df, df_info=nonneuro_obs_info, dep_var="status")
all_limma <- runLimma(df=all_exp_df, df_info=obs_info, dep_var="status")
#  Below 'limma_result' may be replaced with other limma results, 
# as long as data structure is same as neuro_limma
limma_result <- nonneuro_limma
limma_result <- neuro_limma
limma_result <- all_limma
write.csv(neuro_limma,
          file.path(RESULT_DIR, "neuro_limma.csv"),
          row.names = F,
          quote = F)
write.csv(nonneuro_limma,
          file.path(RESULT_DIR, "nonneuro_limma.csv"),
          row.names = F.
          quote = F)




# GP network
#######################
# Gene-Pathway netowork information
#######################
wp.hs.gmt <- rWikiPathways::downloadPathwayArchive(organism="Homo sapiens", format = "gmt")
wp2gene <- clusterProfiler::read.gmt(wp.hs.gmt)
#wp2gene <- clusterProfiler::read.gmt("wikipathways-20201110-gmt-Homo_sapiens.gmt")
wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
wpid2name <- unique(wp2gene %>% dplyr::select(wpid,name)) #TERM2NAME

source("getReactome.R")
wpid2gene <- rbind(wpid2gene, reactome2gene)
wpid2name <- rbind(wpid2name, reactome2name)


#######################
# Pathway Over-Representation Analysis & Gene Set Enrichment Analysis (Begin)
#######################
# For comparison, we want range of logFC and p value to match weight
# sort genes based on average path weight sum
# In order to run the script, need combined wpid2gene & wpid2name, and limma_result

## IMPORTANT: DEPENDING ON DIFFERENT LIMMA_RESULT
source("pw_ora.R")
print(head(enrich_wp[,c(1,2,5)]))
#######################
# Pathway Over-Representation Analysis & Gene Set Enrichment Analysis (End)
#######################


#######################
# create Graph (Begin)
#######################
# For comparison, we want range of logFC and p value to match weight
# sort genes based on average path weight sum
# In order to run the script, need combined wpid2gene & wpid2name, and limma_result
## IMPORTANT: DEPENDING ON DIFFERENT LIMMA_RESULT

limma_result <- nonneuro_limma
source("create_graph.R")
nonneuro_g <- g
nonneuro_merged_edge <- merged_edge
nonneuro_merged_node <- merged_node
#
limma_result <- neuro_limma
source("create_graph.R")
neuro_g <- g
neuro_merged_edge <- merged_edge
neuro_merged_node <- merged_node
#
limma_result <- all_limma
source("create_graph.R")
all_g <- g
all_merged_edge <- merged_edge
all_merged_node <- merged_node


#######################
# create Graph (End)
#######################


# We want to select probably about 100-150 pathways
full_paths_list <- list()
full_res_score <- list()
candidate_wp <- wpid2name[wpid2name$wpid != "WP4657",]$wpid
count <- 0
for (one_name in candidate_wp){
  if (count %% 100 == 0) print(paste("We are at:", count))
  count <- count + 1
  one_res <- calculatePathwayScore(input_g=nonneuro_g,
                                   another_pathway_node=one_name,
                                   WEIGHT_THRESHOLD=0.9,
                                   print_path = FALSE)
  # full_res_score is used to calculate overall distribution of path scores
  print(paste("pathway score of", one_name))
  print(one_res[[2]])
  if (length(one_res[[1]])>0){
    full_res_score[[(length(full_res_score)+1)]] <- one_res[[1]]
    one_path_df <- filterPath(one_res)
    full_paths_list[[(length(full_paths_list)+1)]] <- one_path_df
  }
}

nonneuro_full_paths_list <- full_paths_list
##########neuro_full_paths_list <- full_paths_list


res_rows <- sapply(full_paths_list, function(x) return(nrow(x)))
#full_paths_df <- bind_rows(full_paths_list)
# Calculate the total number of genes per pathway

#######################
# Calculate Pathway P-value (Begin)
#######################

# It is better to directly run the script, for neuro and non-neuro, respectively

# g <- neuro_g
# merged_edge <- neuro_merged_edge
# merged_node <- neuro_merged_node
# source("calculate_pval.R")
# g <- nonneuro_g
# merged_edge <- nonneuro_merged_edge
# merged_node <- nonneuro_merged_node
# source("calculate_pval.R")

# result
head(neuro_pval_df)
head(nonneuro_pval_df)

#######################
# Calculate Pathway P-value (End)
#######################


# TODO: Select pathways based on p-values
# Some kind of sensitivity analysis
neuro_pval_sensitivity <- numeric()
nonneuro_pval_sensitivity <- numeric()
both_pval_sensitivity <- numeric()
for (i in seq(0.001, 0.1, by=0.001)){
  neuro_pval_sensitivity <- c(neuro_pval_sensitivity, nrow(subset(neuro_pval_df, pvalue<i)))
  nonneuro_pval_sensitivity <- c(nonneuro_pval_sensitivity, nrow(subset(nonneuro_pval_df, pvalue<i)))
  both_pval_sensitivity <- c(both_pval_sensitivity,
                             nrow(subset(neuro_pval_df, pvalue<i))+nrow(subset(nonneuro_pval_df, pvalue<i)))
}
plot(seq(0.001, 0.1, by=0.001), both_pval_sensitivity)
#plot(seq(0.001, 0.1, by=0.001), neuro_pval_sensitivity)
#plot(seq(0.001, 0.1, by=0.001), nonneuro_pval_sensitivity)



# Values are selected based on p-value
neuro_sel_val <- rownames(subset(neuro_pval_df, pvalue<0.05))
nonneuro_sel_val <- rownames(subset(nonneuro_pval_df, pvalue<0.05))



# wp_len <- sapply(unique(unlist(sapply(full_res_score, names))), function(x) return(nrow(subset(wpid2gene, wpid==x))))
# norm_val <- res_rows / wp_len
# # 0.01 is a bit arbitrary
# sel_val <- norm_val[norm_val>0.01]



#######################
# Pathway Interaction Threshold Sensitivity Analysis (Begin)
#######################
# Increase path significance threshold 
# and see how they are collected. (use Imax=0.9)
#
# source("sensitivity_analysis.R")
# 
#######################
# Pathway Interaction Threshold Sensitivity Analysis (End)
#######################


# summarize result, extend with drug pathway
sel_val <- neuro_sel_val
g <- neuro_g
source("collect_final_result.R")
neuro_sel_path_df <- sel_path_df

sel_val <- nonneuro_sel_val
g <- nonneuro_g
source("collect_final_result.R")
nonneuro_sel_path_df <- sel_path_df

sel_path_dfs <- list(neuro_sel_path_df, nonneuro_sel_path_df)
#save(sel_path_dfs, file=file.path(DATA_DIR, "sel_path_dfs_005.rda"))

getNumGenes <- function(input_df, wp_tobe_rm=wpid2name$wpid){
  # Calculate number of genes from
  # sel_path_df data frame
  input_nodes <- unique(c(input_df$source, input_df$target))
  filtered_genes <- intersect(input_nodes, wp_tobe_rm)
  return(length(filtered_genes))
}
getNumGenes(neuro_sel_path_df)
getNumGenes(nonneuro_sel_path_df)

#######################
# Visualize Initial pathway result (Begin)
#######################
# Using filtered path, create a graph for visualization

sel_path_df <- neuro_sel_path_df
merged_node <- neuro_merged_node
limma_result <- neuro_limma
one_prefix <- "Neuro"
source("visualize_path.R")
neuro_interaction_g <- interaction_g
neuro_interaction_subg <- interaction_subg

sel_path_df <- nonneuro_sel_path_df
merged_node <- nonneuro_merged_node
limma_result <- nonneuro_limma
one_prefix <- "Non-neuro"
source("visualize_path.R")
nonneuro_interaction_g <- interaction_g
nonneuro_interaction_subg <- interaction_subg

#######################
# Visualize Initial pathway result (End)
#######################



#######################
# Check Hubs (Begin)
#######################

one_ig <- neuro_interaction_g
one_subig <- neuro_interaction_subg
source("check_hubs.R")
neuro_gene_rank_df <- gene_rank_df

one_ig <- nonneuro_interaction_g
one_subig <- nonneuro_interaction_subg
source("check_hubs.R")
nonneuro_gene_rank_df <- gene_rank_df

#######################
# Check Hubs (End)
#######################



##########################################
#######
# Test Some of characteristics (begin)
#######
neuro_path_nodes <- unique(c(neuro_sel_path_df$source, neuro_sel_path_df$target))
nonneuro_path_nodes <- unique(c(nonneuro_sel_path_df$source, nonneuro_sel_path_df$target))

# total_neuro_genes comes from pw_ora.R
neuro_path_genes <- intersect(total_neuro_genes, neuro_path_nodes)
nonneuro_path_genes <- intersect(total_neuro_genes, nonneuro_path_nodes)

neuro_path_entrez <- subset(neuro_limma, entrez_id %in% neuro_path_genes)$entrez_id
nonneuro_path_entrez <- subset(nonneuro_limma, entrez_id %in% nonneuro_path_genes)$entrez_id
neuro_path_hgnc <- subset(neuro_limma, entrez_id %in% neuro_path_genes)$hgnc_symbol
nonneuro_path_hgnc <- subset(nonneuro_limma, entrez_id %in% nonneuro_path_genes)$hgnc_symbol
excl_path_hgnc <- setdiff(neuro_path_hgnc, nonneuro_path_hgnc)
common_path_hgnc <- intersect(neuro_path_hgnc, nonneuro_path_hgnc)
out_path_hgnc <- setdiff(nonneuro_path_hgnc, neuro_path_hgnc)

#length(unique(neighbors(neuro_interaction_subg, excl_path_hgnc[[2]], "all")))

# excl_path_hgnc contains interactomes detected by neuro-limma and neuro-g
collected_nodes <- excl_path_hgnc
for (excl_g in excl_path_hgnc){
  print(paste("Our EXCL gene is", excl_g))
  selected_neighbors <- neighbors(neuro_interaction_g, excl_g, "all")$name
  print(selected_neighbors)
  collected_nodes <- c(collected_nodes, setdiff(selected_neighbors, wpid2name$wpid))
}
collected_nodes <- c(collected_nodes, "WP4657")
neuro_subgraph <- induced_subgraph(neuro_interaction_g, collected_nodes)
RCy3::createNetworkFromIgraph(neuro_subgraph,
                             title = "Neuro Excl. Subgraph (with first neighbors)",
                             collection = "Filtered Network with threshold 1.4 (for Paper)")
toggleGraphicsDetails()
RCy3::loadTableData(neuro_limma, data.key.column = "ensembl_id", table.key.column = "ensembl")



pith_subgraph <- induced_subgraph(neuro_interaction_g, c(excl_path_hgnc, 'CRKL', 'GRAP2', 'WP4657'))
undirect_pith_subgraph <- as.undirected(pith_subgraph, mode="collapse")
RCy3::createNetworkFromIgraph(undirect_pith_subgraph,
                              title = "Neuro Excl. Subgraph (with CRKL/GRAP2/WP4657)",
                              collection = "Filtered Network with threshold 1.4 (for Paper)")
toggleGraphicsDetails()
RCy3::loadTableData(neuro_limma, data.key.column = "ensembl_id", table.key.column = "ensembl")
#######
# Test Some of characteristics (End)
#######
##########################################






#####################
## Ranking Pathways Based on Hub-ness
one_g <- nonneuro_interaction_g
one_g_pathways <- V(one_g)$name[V(one_g)$type=="pathway"]
num_neighbors <- numeric()
one_wps <- character()
one_wp_descs <- character()
for (one_wp in one_g_pathways){
  num_neighbors_one_wp <- length(unique(neighbors(one_g, one_wp, "all")$name))
  one_wps <- c(one_wps, one_wp)
  num_neighbors <- c(num_neighbors, num_neighbors_one_wp)
  one_wp_descs <- c(one_wp_descs, wpid2name[wpid2name$wpid==one_wp,]$name)
}
wp_rank_df <- data.frame(wpid=one_wps,
                         description=one_wp_descs,
                         num_neighbors=num_neighbors,
                         stringsAsFactors = F)
wp_rank_df <- wp_rank_df[order(wp_rank_df$num_neighbors, decreasing = T),]
# write.csv(wp_rank_df,
#           file.path(RESULT_DIR, "neuro_pathway_verticies_by_neighbors.csv"),
#           row.names = F)
# write.csv(wp_rank_df,
#           file.path(RESULT_DIR, "nonneuro_pathway_verticies_by_neighbors.csv"),
#           row.names = F)


