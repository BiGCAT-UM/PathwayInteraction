# functions.R
# Functions used for path

runLimma <- function(df, df_info, dep_var){
  # Run limma, summarize result by 
  # calculating median of logFC per each Entrez ID.
  # Next, add ensembl & hgnc ids and 
  # Finally, aadd weights
  dependent_variable <- df_info[,dep_var]
  design <- model.matrix(~0+dependent_variable)
  colnames(design) <- c("control", "disease")
  print("Calculating Correlation...")
  corfit <- duplicateCorrelation(df, design, block=df_info$family_group)
  fit <- lmFit(df, design, block=df_info$family_group, correlation=corfit$consensus)
  fit <- eBayes(fit)
  #compute the contrast fits
  cont.matrix <- makeContrasts(disease-control, levels=design)
  contrast.fit <- contrasts.fit(fit, cont.matrix)
  contrast.fit <- eBayes(contrast.fit)
  print("running Limma...")
  fit_ordered_by_p_val <- topTable(contrast.fit,
                                   adjust.method="BH",
                                   coef=1,
                                   number=dim(df)[1],
                                   resort.by="P")
  # Get illumina ID -> Entrez ID mapping
  probe_to_entrez <- data.frame(Gene=unlist(mget(x = rownames(df),
                                                 envir = illuminaHumanv4ENTREZID)))
  probe_to_entrez$probe_id <- rownames(probe_to_entrez)
  probe_to_entrez <- probe_to_entrez[rownames(fit_ordered_by_p_val),]
  fit_ordered_by_p_val$entrez_id <- probe_to_entrez$Gene
  complete_fit <- fit_ordered_by_p_val[complete.cases(fit_ordered_by_p_val),]
  # Median was obtained for eaach performed result
  getMedProb <- function(one_id){
    # retuen median of each probes
    sub_df <- subset(complete_fit, entrez_id==one_id)
    return(subset(sub_df, abs(logFC)==quantile(abs(sub_df$logFC), p = 0.5, type = 1)) )
  }
  print("Summarizing per each Entrez...")
  res_list <- lapply(unique(complete_fit$entrez_id), function(x) getMedProb(one_id=x))
  res <- bind_rows(res_list)
  rownames(res) <- res$entrez_id

  print("Add Ensembl and HGNC...")
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  gene_ids_map <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                        filters = "entrezgene_id",
                        values = res$entrez_id,
                        mart = mart
  )
  gene_ids_map$entrezgene_id <- as.character(gene_ids_map$entrezgene_id)
  gene_ids_map <- gene_ids_map[!duplicated(gene_ids_map$hgnc_symbol),]
  colnames(gene_ids_map) <- c("ensembl_id", "entrez_id", "hgnc_symbol")
  limma_result <- merge(res, gene_ids_map, by.x="entrez_id")
  ALPHA <- 2
  MU <- 3
  transformed_t <- 1 - (1 / (1 + exp(-ALPHA*(abs(limma_result$t) - MU))))
  limma_result$weight <- 1 - (1 / (1 + exp(-ALPHA*(abs(limma_result$t) - MU))))
  return(limma_result)
}


createModifiedGraph <- function(reference_edges,
                                reference_nodes,
                                pathway_node_to_remove,
                                to_match=TRUE){
  # Function to calculate P value
  # If to_match==TRUE, sample the same number of nodes from those below 0.8 and over 0.8
  # number of genes to choose
  nodes_to_replace <- subset(reference_edges,source==pathway_node_to_remove)
  num_nodes_to_choose <- nrow(nodes_to_replace)
  edges_to_add_to <- subset(reference_edges,
                            source!=pathway_node_to_remove & target!=pathway_node_to_remove)
  
  # a way to remove gene-pathway nodes
  # (for example, pathway nodes have weight==NA)
  nodes_to_choose <- reference_nodes[complete.cases(reference_nodes),]
  
  # to_match==TRUE means matching by range (threshold 0.8)
  if (to_match){
    nodes_small_to_choose <- subset(nodes_to_choose, weight<=0.8)
    nodes_big_to_choose <- subset(nodes_to_choose, weight>0.8)
    num_nodes_small <- nrow(nodes_to_replace[nodes_to_replace$weight<=0.8,])
    num_nodes_big <- nrow(nodes_to_replace[nodes_to_replace$weight>0.8,])
    
    sampled_nodes_small <- nodes_small_to_choose[sample(rownames(nodes_small_to_choose),
                                                        num_nodes_small),][,c("name", "weight")]
    sampled_nodes_big <- nodes_big_to_choose[sample(rownames(nodes_big_to_choose),
                                                    num_nodes_big),][,c("name", "weight")]
    sampled_nodes <- rbind(sampled_nodes_small, sampled_nodes_big)
  } else{
    sampled_nodes <- nodes_to_choose[sample(rownames(nodes_to_choose),
                                            num_nodes_to_choose),][,c("name", "weight")]
  }
  sampled_edge1 <- data.frame(source=rep(pathway_node_to_remove, nrow(sampled_nodes)),
                              target=sampled_nodes$name,
                              type=rep("gene_pathway", nrow(sampled_nodes)),
                              weight=sampled_nodes$weight,
                              stringsAsFactors = F)
  sampled_edge2 <- data.frame(source=sampled_nodes$name,
                              target=rep(pathway_node_to_remove, nrow(sampled_nodes)),
                              type=rep("gene_pathway", nrow(sampled_nodes)),
                              weight=sampled_nodes$weight,
                              stringsAsFactors = F)
  sampled_edges <- rbind(sampled_edge1, sampled_edge2)
  # not reactly random, but choose candidate pathways with similar weights
  modified_edge <- rbind(sampled_edges, edges_to_add_to)
  modified_g <- graph_from_data_frame(modified_edge,
                                      directed=TRUE,
                                      vertices=reference_nodes)
  return(modified_g)
}

calculatePathwayScore <- function(input_g,
                                  one_pathway_node="WP4657",
                                  another_pathway_node,
                                  print_path=FALSE,
                                  WEIGHT_THRESHOLD=0.9){
  # print(paste("New model", another_pathway_node, "is analyzed"))
  subg <- induced_subgraph(input_g,
                           c(one_pathway_node, another_pathway_node,
                             V(input_g)[type=="gene_product"]$name))
  # one shortest path
  paths <- shortest_paths(subg, one_pathway_node,
                          another_pathway_node, output = "epath")
  # get information of the first shortest path
  path_scores <- numeric()
  # path_lengths <- numeric()
  path_records <- list()
  one_path <- paths$epath[[1]]
  one_path_score <- sum(one_path$weight)
  one_path_length <- length(one_path)
  if (print_path) {
    print(one_path)
    print(one_path_score)
  }
  if (one_path_length==2 & 
      one_path_score < WEIGHT_THRESHOLD & 
      print_path) {
    print(paste(another_pathway_node, "has length 2 and <", WEIGHT_THRESHOLD))
  }
  while ((one_path_score < WEIGHT_THRESHOLD) & (one_path_length>1)){
    edges_to_remove <- numeric()
    if (one_path_length==2){
      edges_to_remove <- tail(one_path, 1)
      subg <- delete_edges(subg, edges_to_remove)
      # get the next path, wihtout updating path_score/path_length
      paths <- shortest_paths(subg, one_pathway_node, another_pathway_node, output = "epath")
      one_path <- paths$epath[[1]]
      one_path_score <- sum(one_path$weight)
      one_path_length <- length(one_path)
      if (print_path) {
        print(one_path) 
        print(paste("After removing it, new path score is:", one_path_score))
      }
      next
    }
    if (print_path){
      print("This is one stored path")
      print(one_path)
    }
    # Record path information 
    path_scores <- c(path_scores, one_path_score)
    path_records[[length(path_records)+1]] <- c(tail_of(subg, one_path)$name[1], head_of(subg, one_path)$name)
    # Remove edges and create an updated subgraph 
    for (one_edge_idx in one_path[2:(length(one_path)-1)]){
      edges_to_remove <- c(edges_to_remove, one_edge_idx)
    }
    subg <- delete_edges(subg, edges_to_remove)
    paths <- shortest_paths(subg, one_pathway_node, another_pathway_node, output = "epath")
    one_path <- paths$epath[[1]]
    # Update criteria (path score and path length)
    # to be checked by While loop
    one_path_score <- sum(one_path$weight)
    one_path_length <- length(one_path)
    if (print_path) print(one_path)
  }
  # After while loop, name the paths using the pathway name
  if (length(path_scores)>0){
    names(path_scores) <- rep(another_pathway_node, length(path_scores))
  }
  # return collected scores and and paths
  return(list(path_scores, path_records))
}

# Following code needs elements from functionPathwayInteraction.R
collectPaths <- function(one_res_path){
  # This function summarizes one string vector (path),
  # and returns one dataframe
  path_from <- head(one_res_path, -1)
  path_to <- tail(one_res_path, -1)
  path_df <- data.frame(source=path_from,
                        target=path_to,
                        stringsAsFactors = F)
  return(path_df)
}

checkLength <- function(one_p){
  return( length(one_p) - 1)
}

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

getDesc <- function(one_row, path_info){
  if (one_row$type=="gene_product"){
    return (one_row$hgnc_symbol)
  } else {
    return (path_info[path_info$wpid==one_row$entrez,]$name)
  }
}

getEdgeName <- function(one_name,
                        wpid_info=wpid2name,
                        node_info=orig_node_df){
  # get HGNC symbol for entrez and return
  # if it is WPID, return as it is
  if (one_name %in% wpid_info$wpid) return(one_name)
  else return(node_info[node_info$name==one_name,]$hgnc_symbol)
}


getOneRowName <- function(one_row){
  # Run getEdgeName for each row
  return(sapply(one_row, function(x) getEdgeName(one_name=x)))
}



preparePathViz <- function(input_edge_df,
                           orig_node_df,
                           pathway_info){
  # Create a graph based on input_edge_df
  path_edge_df <- input_edge_df[!duplicated(input_edge_df),]
  path_edge_mat <- apply(path_edge_df, 2, function(x) getOneRowName(x))
  edge_df <- as.data.frame(path_edge_mat)
  path_nodes <- unique(c(path_edge_df$source, path_edge_df$target))
  node_df <- subset(orig_node_df, name %in% path_nodes)
  rownames(node_df) <- NULL
  colnames(node_df) <- c("entrez", "ensembl", "hgnc_symbol", "weight", "type")
  node_df$desc <- sapply(node_df$entrez, function(x) getDesc(one_row=node_df[node_df$entrez==x,], path_info=pathway_info))
  node_df$name <- ifelse(node_df$type=="gene_product", node_df$hgnc_symbol, node_df$entrez)
  node_df <- node_df[,c("name", "entrez", "ensembl", "weight", "type")]
  # create iGraph object
  res_g <- graph_from_data_frame(edge_df,
                                 directed=TRUE,
                                 vertices=node_df)
  return(res_g)
}


s <- table(edge_df$source)
s[order(s, decreasing=TRUE)]
