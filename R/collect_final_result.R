# collect_final_result.R
# Collect final results (path and filtered pathway list)
sel_paths_list <- list()
for (one_name in sel_val){
  one_res <- calculatePathwayScore(input_g=g,
                                   another_pathway_node=one_name,
                                   WEIGHT_THRESHOLD=1.4,
                                   print_path = TRUE)
  one_path_df <- filterPath(one_res)
  sel_paths_list[[(length(sel_paths_list)+1)]] <- one_path_df
}
sel_path_df <- bind_rows(sel_paths_list)