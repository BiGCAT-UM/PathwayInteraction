# vizPathInt2.R
# Polished version of visualizePathwayInteraction.R

#clean workspace
rm(list=ls())

library(dplyr)
library(ggplot2)
library(igraph)
library(rstudioapi)
library(RCy3)
library(rWikiPathways)
library(xlsx)

BASE_DIR <- "~/Desktop/BiGCAT/project_22q11ds/microarray_analysis/pathway_interaction"
DATA_DIR <- file.path(BASE_DIR, "data")
setwd(BASE_DIR)
source("functions.R")
# summary of selected candidate pathwyas
summary_df <- read.table(file=file.path(DATA_DIR, "sunmary_df_scores_p_values.csv"),
                         sep=",",
                         header=TRUE,
                         row.names=1)
limma_result <- read.csv(file=file.path(DATA_DIR, "limma_result_geomean.csv"))
limma_result$entrez_id <- as.character(limma_result$entrez_id)
# mapping for pathwya names
pathway_map <- read.csv(file=file.path(DATA_DIR,"pathway_map.csv"), row.names = 1, stringsAsFactors = F)
load(file=file.path(DATA_DIR, "merged_full_edges_for_kelder.rda"))
load(file=file.path(DATA_DIR, "merged_full_nodes_for_kelder.rda"))
merged_full_node_unique$hgnc_symbol <- as.character(merged_full_node_unique$hgnc_symbol)
g <- graph_from_data_frame(merged_full_edge, directed=TRUE, vertices=merged_full_node_unique)

one_name <- rownames(summary_df)[1]
one_res <- calculatePathwayScore(input_g=g,
                                 another_pathway_node=one_name,
                                 WEIGHT_THRESHOLD=1.4)
# check the length (number of edges) for each path
# only accept 0.9
checkLength <- function(one_p){
  return( length(one_p) - 1)
}
sapply(one_res[[2]], function(x) checkLength(x))

# sorting out paths
filterPath <- function(one_result){
  # List of paths, adjusted for length-threshold
  path_df <- data.frame(score=one_result[[1]],
                        plength=sapply(one_result[[2]],
                                       checkLength)
                        )
  path_df$thres <- 0.9 + (path_df$plength-3)*0.5
  filtered_paths <- one_result[[2]][path_df$score < path_df$thres]
  path_list <- lapply(filtered_paths, collectPaths)
  path_df <- bind_rows(path_list)
  return(path_df)
}
path_df <- filterPath(one_res)
# Crete a list of filtered paths
full_paths_list <- list()
for (one_name in rownames(summary_df)){
  one_res <- calculatePathwayScore(input_g=g,
                                   another_pathway_node=one_name,
                                   WEIGHT_THRESHOLD=1.4)
  one_path_df <- filterPath(one_res)
  full_paths_list[[(length(full_paths_list)+1)]] <- one_path_df
}

# combining resuls
all_paths <- bind_rows(full_paths_list)
res_paths <- all_paths[!duplicated(all_paths),]
rownames(res_paths) <- NULL
colnames(res_paths) <- c("source", "target")
save(res_paths, file=file.path(DATA_DIR, "res_paths.rds"))

load(file=file.path(DATA_DIR, "res_paths.rds"))

# if want, choose only unique edges
#
edge_df <- res_paths
orig_node_df=merged_full_node_unique
pathway_info=pathway_map
#ent_to_ens <- limma_result[,c("entrez_id", "ensembl_id", "hgnc_symbol")]



preparePathViz <- function(edge_df,
                           orig_node_df=merged_full_node_unique,
                           pathway_info=pathway_map){
  # Create a graph based on edge df
  #path_edge_df <- edge_df %>% left_join(merged_full_edge, by=c("from", "to"))
  #colnames(orig_node_df) <- c("", "ensembl_id", "name", "weight", "logFC", "P.Value", "type")
  path_edge_df <- edge_df[!duplicated(edge_df),]
  path_nodes <- unique(c(path_edge_df$source, path_edge_df$target))
  node_df <- subset(orig_node_df, name %in% path_nodes)
  rownames(node_df) <- NULL
  node_df$label <- node_df$hgnc_symbol
  node_df$desc <- node_df$hgnc_symbol
  node_df[node_df$type=="pathway",]$label <- node_df[node_df$type=="pathway",]$name
  node_df[node_df$type=="pathway",]$desc <- pathway_info[node_df[node_df$type=="pathway",]$name,]
  # create iGraph object
  res_g <- graph_from_data_frame(path_edge_df,
                                 directed=TRUE,
                                 vertices=node_df)
  return(res_g)
}

interaction_g <- preparePathViz(edge_df=res_paths)


# visualize the full results using RCy3
RCy3::cytoscapePing()
RCy3::createNetworkFromIgraph(interaction_g,
                              title = "One Interaction",
                              collection = "WP4657 and One Candidate Pathway")
toggleGraphicsDetails()
setNodeLabelMapping(table.column="label", network = "One Interaction")
RCy3::loadTableData(limma_result, data.key.column = "ensembl_id", table.key.column = "ensembl_id")

# visulize EXCLUDING pathway nodes (to focus on genes starting from WP4657)
s <- V(interaction_g)[[V(interaction_g)$type=="pathway"]]
names(s)
nodes_to_remove <- names(s)[names(s) != "WP4657"]
interaction_subg <- delete_vertices(interaction_g, nodes_to_remove)
RCy3::createNetworkFromIgraph(interaction_subg,
                              title = "Sub Interaction",
                              collection = "WP4657 and One Candidate Pathway")
toggleGraphicsDetails()
setNodeLabelMapping(table.column="label", network = "Sub Interaction")
RCy3::loadTableData(limma_result, data.key.column = "ensembl_id", table.key.column = "ensembl_id")


degree(interaction_subg, V(interaction_subg)[[1]])
V(interaction_subg)[[1]]$label

getDeg <- function(one_g){
  gene_name <- character()
  gene_deg <- numeric()
  for (one_v in V(one_g)){
    gene_name <- c(gene_name, V(one_g)[[one_v]]$label)
    gene_deg <- c(gene_deg, degree(one_g, one_v))
  }
  names(gene_deg) <- gene_name
  res_df <- data.frame(deg=gene_deg,
                       gene=gene_name,
                       stringsAsFactors = F)
  # print(res_df)
  res_df <- res_df[order(res_df$deg, decreasing = T),]
  # print(class(res_df))
  return(res_df)
}

g_deg <- getDeg(one_g=interaction_g)
subg_deg <- getDeg(one_g=interaction_subg)

getWPName(wp="WP4400")

getWPName <- function(wp, input_pathway_map=pathway_map){
  return(input_pathway_map[wp,])
}



genes <- subset(subg_deg, gene != "WP4657")$gene
sub_v <- V(interaction_subg)[[V(interaction_subg)$label %in% genes]]
one_g <- induced_subgraph(interaction_subg, sub_v)
V(one_g)$name <- V(one_g)$label
RCy3::cytoscapePing()
RCy3::createNetworkFromIgraph(one_g,
                              title = "Sub Node Graph for Extension",
                              collection = "WP4657 and One Candidate Pathway")


df_g <- read_graph(file=file.path(DATA_DIR, "drugtarget_g.graphml"),
                   format="graphml")

drug_df <- V(df_g)[V(df_g)$CTL.Type=="drug"]
drug_df$name
neighbors(df_g, drug_df$name[1])$name

s <- unlist(sapply(drug_df$name, function(x) neighbors(df_g, x)$name))

[1] "GAA"      "IGFBP7"   "NDUFC2"   "ITGAL"    "SERPIND1" "MGST1"    "CAT"      "SNAP25"   "TP53"     "FPR1"     "CPNE1"   
[12] "ACP3"     "CEACAM1"  "ATP8A1"   "KYAT3"    "LDLR"     "B4GALT1" 

## "CPNE1" "B4GALT1"

#### load(filtered_g_list, file=file.path(DATA_DIR, "filtered_graphs_for_viz.rda"))
# Do this for all candidates
# filtered_g_list <- list()
# count <- 0
# for (one_name in rownames(summary_df)){
#   count <- count + 1
#   print(paste("We are now on", one_name, ", count:", count))
#   one_res <- calculatePathwayScore(input_g=g,
#                                    another_pathway_node=one_name,
#                                    WEIGHT_THRESHOLD=1.4)
#   # filter each path
#   filtered_paths <- filterPath(one_res)
#   path_list <- lapply(filtered_paths, collectPaths)
#   path_df <- bind_rows(path_list)
#   # prepare data frames and creaate refined interaction graph
#   interaction_g <- preparePathViz(edge_df=path_df)
#   # Collect all data
#   filtered_g_list[[(length(filtered_g_list)+1)]] <- interaction_g
# }
#
#save(filtered_g_list, file=file.path(DATA_DIR, "filtered_graphs_for_viz.rda"))

fgl <- filtered_g_list
names(fgl) <- rownames(summary_df)
selected_p <- character()
count <- 0
for (one_g in fgl){
  count <- count + 1
  if ("TXNRD2" %in% V(one_g)$label){
    print("yes!")
    print(names(fgl)[count])
  }
}
selectGeneG <- function(gene_name, all_g=fgl, print_g=FALSE){
  # select graphs from a list of graphs
  # which contains a specific gene
  res_g <- list()
  count <- 0
  for (one_g in all_g){
    count <- count + 1
    if (gene_name %in% V(one_g)$label){
      res_g[[(length(res_g)+1)]] <- one_g
      if (print_g){
        pathway_code <- names(all_g)[count]
        RCy3::createNetworkFromIgraph(one_g,
                                      title = paste(pathway_code, V(one_g)[name==pathway_code]$label, sep=": "),
                                      collection = paste("Pathways with", gene_name))
      }
    }
  }
  return (res_g)
}

selectGeneG(gene_name="TXNRD2", print_g=TRUE)

