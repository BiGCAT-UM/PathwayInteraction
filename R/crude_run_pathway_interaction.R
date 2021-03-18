# crude_run_pathway_interaction.R
# Run all pathways the list of wpids, 
# and without filtering the network

full_paths_list <- list()
full_res_score <- list()
candidate_wp <- wpid2name[wpid2name$wpid != "WP4657",]$wpid
count <- 0
for (one_name in candidate_wp){
  if (count %% 100 == 0) print(paste("We are at:", count))
  count <- count + 1
  one_res <- calculatePathwayScore(input_g=g,
                                   another_pathway_node=one_name,
                                   WEIGHT_THRESHOLD=0.9,
                                   print_path = FALSE)
  # full_res_score is used to calculate overall distribution of path scores
  # print(paste("pathway score of", one_name))
  # print(one_res[[2]])
  if (length(one_res[[1]])>0){
    full_res_score[[(length(full_res_score)+1)]] <- one_res[[1]]
    one_path_df <- filterPath(one_res)
    full_paths_list[[(length(full_paths_list)+1)]] <- one_path_df
  }
}