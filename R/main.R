# main.R

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
WP_DIR <- "~/Desktop/BiGCAT/project_22q11ds/microarray_analysis/pathway_interaction/data/pathway_Feb2021/"
SUPPL_DIR <- "~/Desktop/BiGCAT/project_22q11ds/microarray_analysis/pathway_interaction/paper_supplemental"
setwd("~/Desktop/BiGCAT/project_22q11ds/microarray_analysis/pathway_interaction")
source("functions.R")

#######################
# expression data & Limma (Med Prob)
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
# Below 'limma_result' may be replaced with other limma results, 
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
          row.names = F,
          quote = F)
write.csv(all_limma,
          file.path(RESULT_DIR, "all_limma.csv"),
          row.names = F,
          quote = F)
#######################
# Limma data load (so above not necessary)
#######################
neuro_limma <- read.csv(file.path(RESULT_DIR, "neuro_limma.csv"))
nonneuro_limma <- read.csv(file.path(RESULT_DIR, "nonneuro_limma.csv"))
all_limma <- read.csv(file.path(RESULT_DIR, "all_limma.csv"))


# GP network
#######################
# Gene-Pathway netowork information
#######################
# wp.hs.gmt <- rWikiPathways::downloadPathwayArchive(organism="Homo sapiens", format = "gmt")
# wp2gene <- clusterProfiler::read.gmt(wp.hs.gmt)
# WikiPathways from .gpml folder
wp_files <- list.files(file.path(WP_DIR, "wikipathways-20210110-gpml-Homo_sapiens"))
wp_pathways <- sapply(wp_files, getSplit)
names(wp_pathways) <- NULL

# WikiPathways from .gmt
wp2gene <- clusterProfiler::read.gmt(file.path(WP_DIR, "wikipathways-20210110-gmt-Homo_sapiens.gmt"))
wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
# Subsetting using intersection  of .gmt and .gpml
wp2gene <- subset(wp2gene, wpid %in% intersect(unique(wp2gene$wpid), wp_pathways))
wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
wpid2name <- unique(wp2gene %>% dplyr::select(wpid,name)) #TERM2NAME
#wp_files <- list.files(file.path(WP_DIR, "wikipathways-20210110-gpml-Homo_sapiens"))
wp_pathways <- sapply(wp_files, getSplit)
names(wp_pathways) <- NULL
# check the number of pathways (2 Feb 2021, 628 pathways)
length(intersect(unique(wp2gene$wpid), wp_pathways))

source("getReactome.R")
wpid2gene <- rbind(wpid2gene, reactome2gene)
wpid2name <- rbind(wpid2name, reactome2name)

# final_wpid <- list(wpid2gene, wpid2name)
# save(final_wpid, file=file.path(RESULT_DIR, "final_wpid_info.rda"))


#######################
# Pathway Over-Representation Analysis (Begin)
#######################
## RUN WITH PATHVISIO ###
#######################
# Pathway Over-Representation Analysis (End)
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


# (OPTIONAL): Crude run of pathway interaction 
g <- neuro_g
source("crude_run_pathway_interaction.R")
neuro_full_paths_list <- full_paths_list

g <- nonneuro_g
source("crude_run_pathway_interaction.R")
nonneuro_full_paths_list <- full_paths_list

# res_rows <- sapply(full_paths_list, function(x) return(nrow(x)))
# full_paths_df <- bind_rows(full_paths_list)



# Calculate the total number of genes per pathway

#######################
# Calculate Pathway P-value (Begin)
#######################
# TODO: run with r+1/n+1 formula
# It is better to manually run the script, 
# for neuro and non-neuro, respectively (till line 82)
 

# g <- neuro_g
# merged_edge <- neuro_merged_edge
# merged_node <- neuro_merged_node
# # source("calculate_pval.R")
neuro_pval_df <- read.csv(file=file.path(RESULT_DIR, "neuro_pval_df_167_2021.csv"), 
                          stringsAsFactors=F)

# g <- nonneuro_g
# merged_edge <- nonneuro_merged_edge
# merged_node <- nonneuro_merged_node
# source("calculate_pval.R")
nonneuro_pval_df <- read.csv(file=file.path(RESULT_DIR, "nonneuro_pval_df_255_2021.csv"), 
                             stringsAsFactors=F)

# result
head(neuro_pval_df)
head(nonneuro_pval_df)

#######################
# Calculate Pathway P-value (End)
#######################


# Select pathways based on p-values
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
plot(seq(0.001, 0.1, by=0.001), neuro_pval_sensitivity)
plot(seq(0.001, 0.1, by=0.001), nonneuro_pval_sensitivity)

# Values are selected based on p-value
neuro_sel_val <- subset(neuro_pval_df, pvalue<0.05)$wpid
nonneuro_sel_val <- subset(nonneuro_pval_df, pvalue<0.05)$wpid
# Result of selection
# > length(neuro_sel_val)
# [1] 154
# > length(nonneuro_sel_val)
# [1] 246


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
sel_wp <- subset(neuro_pval_df, pvalue<0.05)$wpid
g <- neuro_g
source("collect_final_result.R")
neuro_sel_path_df <- sel_path_df
# checking if all selected pathways are indeed included
# > sum(neuro_pval_df$wpid %in% neuro_sel_path_df$target)
# [1] 154


sel_wp <- subset(nonneuro_pval_df, pvalue<0.05)$wpid
g <- nonneuro_g
source("collect_final_result.R")
nonneuro_sel_path_df <- sel_path_df
# checking if all selected pathways are indeed included
# > sum(nonneuro_pval_df$wpid %in% nonneuro_sel_path_df$target)
# [1] 246

sel_path_dfs <- list(neuro_sel_path_df, nonneuro_sel_path_df)
save(sel_path_dfs, file=file.path(RESULT_DIR, "sel_path_dfs_005.rda"))


getNumGenes <- function(input_df, wp_tobe_rm=wpid2name$wpid){
  # Calculate number of genes from
  # sel_path_df data frame
  input_nodes <- unique(c(input_df$source, input_df$target))
  filtered_genes <- setdiff(input_nodes, wp_tobe_rm)
  return(length(filtered_genes))
}
getNumGenes(neuro_sel_path_df)
# [1] 116
getNumGenes(nonneuro_sel_path_df)
# [1] 185
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
##
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

write.csv(neuro_gene_rank_df,
          file.path(SUPPL_DIR, "number_of_neighbors_Psychiatric.csv"),
          row.names = F)
write.csv(nonneuro_gene_rank_df,
          file.path(SUPPL_DIR, "number_of_neighbors_Nonpsychiatric.csv"),
          row.names = F)

#######################
# Check Hubs (End)
#######################



##########################################
#######
# Get Neuropsychiatric subgraph with pathway interaction (begin)
#######
neuro_path_nodes <- unique(c(neuro_sel_path_df$source, neuro_sel_path_df$target))
nonneuro_path_nodes <- unique(c(nonneuro_sel_path_df$source, nonneuro_sel_path_df$target))

source("get_neuro_genes.R")
neuro_path_genes <- intersect(total_neuro_genes, neuro_path_nodes)
nonneuro_path_genes <- intersect(total_neuro_genes, nonneuro_path_nodes)

neuro_path_entrez <- subset(neuro_limma, entrez_id %in% neuro_path_genes)$entrez_id
nonneuro_path_entrez <- subset(nonneuro_limma, entrez_id %in% nonneuro_path_genes)$entrez_id
neuro_path_hgnc <- subset(neuro_limma, entrez_id %in% neuro_path_genes)$hgnc_symbol
nonneuro_path_hgnc <- subset(nonneuro_limma, entrez_id %in% nonneuro_path_genes)$hgnc_symbol
excl_path_hgnc <- setdiff(neuro_path_hgnc, nonneuro_path_hgnc)
common_path_hgnc <- intersect(neuro_path_hgnc, nonneuro_path_hgnc)
out_path_hgnc <- setdiff(nonneuro_path_hgnc, neuro_path_hgnc)

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
                             collection = "Filtered Network with threshold 1.4 (for Paper), 14 Feb 2021")
toggleGraphicsDetails()
RCy3::loadTableData(neuro_limma, data.key.column = "ensembl_id", table.key.column = "ensembl")


pith_subgraph <- induced_subgraph(neuro_interaction_g, c(excl_path_hgnc, 'CRKL', 'GRAP2', 'WP4657'))
undirect_pith_subgraph <- as.undirected(pith_subgraph, mode="collapse")
RCy3::createNetworkFromIgraph(undirect_pith_subgraph,
                              title = "Neuro Excl. Subgraph (with CRKL/GRAP2/WP4657)",
                              collection = "Filtered Network with threshold 1.4 (for Paper), 14 Feb 2021")
toggleGraphicsDetails()
RCy3::loadTableData(neuro_limma, data.key.column = "ensembl_id", table.key.column = "ensembl")
#######
# Get Neuropsychiatric subgraph with pathway interaction (End)
#######
##########################################

##########################################
#######
# Get Neuropsychiatric subgraph with pathway ORA (begin)
#######

# run below script (manually, if possible)
# get_ora_network.R


#######
# Get Neuropsychiatric subgraph with pathway ORA (end)
#######


#######
# Combine PI and ORA (Begin)
#######

# pith_source means the current excl interaction subnetwork
pith_source <- V(undirect_pith_subgraph)$entrez
pith_source <- pith_source[pith_source != "WP4657"]
# pith_target: the list of 10 entrez represents exclusive DEGs by ORA
pith_target <- subset(neuro_limma, entrez_id %in% c("135",  "154",  "239",  "694", "2289", "3002", "3304", "3553", "5551", "9588"))

# pith, pith_edges, pith_nodes are current excl interaction subnetwork information
pith <- as_data_frame(pith_subgraph, what = "both")
pith_edges <- pith$edges
pith_nodes <- pith$vertices

one_from <- character()
one_to <- character()
targets <- character()

pith_all_edges <- subset(neuro_merged_edge, source %in% pith_source)
pith_new_edges <- subset(pith_all_edges, target %in% pith_target$entrez_id)
if ( nrow(pith_new_edges)>0 ){
  for (i in seq(1, nrow(pith_new_edges))){
    one_row <- pith_new_edges[i,]
    one_from <- c(one_from, subset(neuro_limma, entrez_id==one_row$source)$hgnc_symbol)
    one_to <- c(one_to, subset(neuro_limma, entrez_id==one_row$target)$hgnc_symbol)
    print(paste("target is", one_to))
    targets <- union(targets, one_to)
  }
  pith_source <- unique(pith_new_edges$target)
  print(paste("new source is", pith_source))
} 
# RUN lines (411-423) ThRICE  TO GET THE FOLLOWINGS:
# > one_from
# [1] "IL2"    "PIK3CA" "IL2"    "PRF1"  
# > one_to
# [1] "IL1B" "IL1B" "PRF1" "GZMB"
# > targets
# [1] "IL1B" "PRF1" "GZMB"

pith_edges_to_add <- data.frame(from=one_from, to=one_to, stringsAsFactors = F)
pith_edges <- rbind(pith_edges, pith_edges_to_add)

# remove duplicated nodes
pith_edges_sort <- as.data.frame(t(apply(pith_edges, 1, sort)))
pith_edges_uniq <- pith_edges_sort[!duplicated(pith_edges_sort),]
colnames(pith_edges_uniq) <- colnames(pith_edges)

pith_target_to_add <- subset(pith_target, hgnc_symbol %in% targets)[,c("hgnc_symbol", "entrez_id", "ensembl_id", "weight")]
colnames(pith_target_to_add) <- c("name", "entrez", "ensembl", "weight")
pith_target_to_add$type <- "gene_product"
rownames(pith_target_to_add) <- pith_target_to_add$name
pith_nodes <- rbind(pith_nodes, pith_target_to_add)

pi_ora_combined_g <- graph_from_data_frame(pith_edges_uniq, directed=FALSE, vertices=pith_nodes)
RCy3::createNetworkFromIgraph(pi_ora_combined_g, 
                              title = "PI_ORA_combined",
                              collection = "Filtered Network with threshold 1.4 (for Paper), 14 Feb 2021")
toggleGraphicsDetails()
RCy3::loadTableData(neuro_limma, data.key.column = "ensembl_id", table.key.column = "ensembl")


#######
# Combine PI and ORA (End)
#######


result_graphs <- list(neuro_g, nonneuro_g, all_g)
save(result_graphs, file=file.path(RESULT_DIR, "all_three_graphs.rda"))






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


