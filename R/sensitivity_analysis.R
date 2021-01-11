# sensitivity_analysis.R
# Run sensitivity analysis, i.e. 
# check what happens if one changes the parameters

# OVERALL SENSITIVITY ANALYSIS IS BASED ON ALL_LIMMA (ALL_G)

#############################
## Change pathway threshold
############################
# Change threshold from 0.8 -> 1.8, to see what happens

thresh_score_sum <- numeric()
all_res_scores <- list()
threshold_range <- seq(0.9, 1.8, 0.1)
candidate_wp <- wpid2name[wpid2name$wpid != "WP4657",]$wpid
wp_list <- candidate_wp
for (thresh in threshold_range){
  res_score <- list()
  # thresh==1.0 means, the result for thresh==0.9 was obtained.
  # we use this as the criterion
  if (thresh == 1.0){
    wp_list <- unique(unlist(sapply(all_res_scores[[1]], names)))
  }
  count <- 0
  for (one_name in wp_list){
    count <- count + 1
    print(paste("We are at threshold", thresh, " and pathway count", count))
    one_res <- calculatePathwayScore(input_g=nonneuro_g,
                                     another_pathway_node=one_name,
                                     WEIGHT_THRESHOLD=thresh)
    # res_score is used to calculate overall distribution of path scores
    if (length(one_res[[1]])>0){
      res_score[[(length(res_score)+1)]] <- one_res[[1]]
    }
  }
  one_sum <- sum(sapply(res_score, function(x) sum(1/x)))
  all_res_scores[[(length(all_res_scores)+1)]] <- res_score
  print(paste("The score sum is:", one_sum))
  thresh_score_sum <- c(thresh_score_sum, one_sum)
}

# nonneuro_g_scores <- all_res_scores
# nonneuro_thresh_score_sum <- thresh_score_sum
# save(nonneuro_g_scores, file=file.path(DATA_DIR, "nonneuro_g_scores.rda"))
# neuro_g_scores <- all_res_scores
# neuro_thresh_score_sum <- thresh_score_sum
# save(neuro_g_scores, file=file.path(DATA_DIR, "neuro_g_scores.rda"))
#all_g_scores <- all_res_scores

compare_sensitivity <- data.frame(val=threshold_range,
                                  neuro_score=neuro_thresh_score_sum,
                                  nonneuro_score=nonneuro_thresh_score_sum)
write.csv(compare_sensitivity, file.path(RESULT_DIR, "compare_sensitivity.csv"), row.names = F)
