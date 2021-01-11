# visualize_path.R
# Visualize path on Cytoscape

# All nodes network
interaction_g <- preparePathViz(input_edge_df=sel_path_df,
                                orig_node_df=merged_node,
                                pathway_info=wpid2name)
RCy3::cytoscapePing()
RCy3::createNetworkFromIgraph(interaction_g,
                              title = paste("All", one_prefix, "Network (with pathways)"),
                              collection = "Filtered Network with threshold 1.4 (for Paper)")
toggleGraphicsDetails()
# load data (limma_result needs to be specified)
RCy3::loadTableData(limma_result, data.key.column = "ensembl_id", table.key.column = "ensembl")

# Gene - only network
wp_to_rm <- V(interaction_g)[[V(interaction_g)$type=="pathway"]]
names(wp_to_rm)
nodes_to_remove <- names(wp_to_rm)[names(wp_to_rm) != "WP4657"]
interaction_subg <- delete_vertices(interaction_g, nodes_to_remove)
RCy3::createNetworkFromIgraph(interaction_subg,
                              title = paste("Gene-only", one_prefix, "Network (without pathways)"),
                              collection = "Filtered Network with threshold 1.4 (for Paper)")
toggleGraphicsDetails()
setNodeLabelMapping(table.column="name", network = paste("Gene-only", one_prefix, "Network (without pathways)"))
RCy3::loadTableData(limma_result, data.key.column = "ensembl_id", table.key.column = "ensembl")