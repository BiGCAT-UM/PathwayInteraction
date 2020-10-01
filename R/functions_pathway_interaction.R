# functions_pathway_interaction.R
# functions to calculate pathway interactions

calculatePathwayScore <- function(input_g=g,
                                  one_pathway_node="WP4657",
                                  another_pathway_node="WP2759"
                                  ){
  ## Calculate pathway scores.
  ## Disregard paths with length 2 or weighted sum > 0.9,
  ## and store all acceptable paths.
  ## Finally, return pathway score
  ## pathway score = sum(1/path_lengths)
  ########################################
  # param: iGraph object: g
  # param: string: one_pathway_node (22q11DS pathway by default)
  # param: string: another_pathway_node (target pathway)
  # param: bool: print_path
  # return: list of pathway scores and paths
  ########################################
  subg <- induced_subgraph(input_g,
                           c(one_pathway_node, another_pathway_node,
                             V(input_g)[type=="gene_product"]$name))
  
  # one shortest path
  paths <- shortest_paths(subg, one_pathway_node,
                          another_pathway_node, output = "epath")
  # get information of the first shortest path
  path_scores <- numeric()
  found_paths <- list()
  #path_lengths <- numeric()
  one_path <- paths$epath[[1]]
  one_path_score <- sum(one_path$weight)
  one_path_length <- length(one_path)

  # if length==2, remove the edge and move to the next while loop
  while ((one_path_score < 0.9) & (one_path_length>1)){
    edges_to_remove <- numeric()
    if (one_path_length==2){
      edges_to_remove <- tail(one_path, 1)
      subg <- delete_edges(subg, edges_to_remove)
      # get the next path, wihtout updating path_score/path_length
      paths <- shortest_paths(subg, one_pathway_node, another_pathway_node, output = "epath")
      one_path <- paths$epath[[1]]
      one_path_score <- sum(one_path$weight)
      one_path_length <- length(one_path)
      next
    }
    path_scores <- c(path_scores, one_path_score)
    found_paths[[length(found_paths)+1]] <- one_path
    #path_lengths <- c(path_lengths, one_path_length)
    # remove edges and create an updated subgraph
    for (one_edge_idx in one_path[2:(length(one_path)-1)]){
      edges_to_remove <- c(edges_to_remove, one_edge_idx)
    }
    subg <- delete_edges(subg, edges_to_remove)
    paths <- shortest_paths(subg, one_pathway_node, another_pathway_node, output = "epath")
    one_path <- paths$epath[[1]]
    one_path_score <- sum(one_path$weight)
    one_path_length <- length(one_path)
  }
  if (length(found_paths)>0){
    # names(path_lengths) <- rep(another_pathway_node, length(path_lengths))
    names(path_scores) <- rep(another_pathway_node, length(path_scores))
  }
  # return collected scores and paths, as long as lenght is < 0.9
  return(list(path_scores, found_paths))
}

#############################################################################################
createModifiedGraph <- function(reference_edges=merged_full_edge,
                                reference_nodes=merged_full_node_unique,
                                pathway_node_to_remove,
                                to_match=FALSE){
  ## Create a modified graph.
  ## Remove pathway edges collected to the node,
  ## and replace with randomly sampled edges
  ## Can be either match by range ((0.8)
  ## or completely random.
  ########################################
  # param: dataframe: reference_edges 
  # param: daraframe: reference_nodes
  # param: string: pathway_node_to_remove
  # param: bool: to_match
  # return: iGraph g:
  ########################################
  # number of genes to choose
  nodes_to_replace <- subset(reference_edges,from==pathway_node_to_remove)
  num_nodes_to_choose <- nrow(nodes_to_replace)
  edges_to_add_to <- subset(reference_edges,
                            from!=pathway_node_to_remove & to!=pathway_node_to_remove)
  
  # a way to remove gene-pathway nodes (for example, weight==NA)
  nodes_to_choose <- reference_nodes[complete.cases(reference_nodes),]
  
  # to_match==TRUE means matching by range (threshold 0.8)
  if (to_match){
    nodes_small_to_choose <- subset(nodes_to_choose, weight<=PATHWAY_SAMPLING_THRESHOLD)
    nodes_big_to_choose <- subset(nodes_to_choose, weight>PATHWAY_SAMPLING_THRESHOLD)
    num_nodes_small <- nrow(nodes_to_replace[nodes_to_replace$weight<=PATHWAY_SAMPLING_THRESHOLD,])
    num_nodes_big <- nrow(nodes_to_replace[nodes_to_replace$weight>PATHWAY_SAMPLING_THRESHOLD,])
    
    sampled_nodes_small <- nodes_small_to_choose[sample(nrow(nodes_small_to_choose),
                                                        num_nodes_small),][,c("name", "weight")]
    sampled_nodes_big <- nodes_big_to_choose[sample(nrow(nodes_big_to_choose),
                                                    num_nodes_big),][,c("name", "weight")]
    sampled_nodes <- rbind(sampled_nodes_small, sampled_nodes_big)
  } else{
    sampled_nodes <- nodes_to_choose[sample(nrow(nodes_to_choose),
                                            num_nodes_to_choose),][,c("name", "weight")]
  }
  sampled_edge1 <- data.frame(from=rep(pathway_node_to_remove, nrow(sampled_nodes)),
                              to=sampled_nodes$name,
                              type=rep("gene_pathway", nrow(sampled_nodes)),
                              weight=sampled_nodes$weight,
                              stringsAsFactors = F)
  sampled_edge2 <- data.frame(from=sampled_nodes$name,
                              to=rep(pathway_node_to_remove, nrow(sampled_nodes)),
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

#############################################################################################


