# sensitivity_analysis.R
# Run sensitivity analysis, i.e. 
# check what happens if one changes the parameters

# OVERALL SENSITIVITY ANALYSIS IS BASED ON ALL_LIMMA (ALL_G)

#############################
## Change pathway threshold
############################
library(foreach)
library(doParallel)
library(parallel)
numCores <- detectCores()-2
numCores
registerDoParallel(numCores)  # use multicore, set to the number of our cores


threshold_range <- seq(0.9, 1.8, 0.1)
g <- neuro_g
thresh_score_sum <- numeric()
candidate_wp <- wpid2name[wpid2name$wpid != "WP4657",]$wpid
wp_list <- candidate_wp
for (thresh in threshold_range){
  #res_score <- list()
  # We'll only use pathways having a significant path for thresh==0.9.
  # Update and fix wp_list accordingly, starting from thresh==1.0
  if (thresh == 1.0){
    wp_list <- one_names
  }
  reses <- foreach (one_name=wp_list) %dopar% {
    one_res <- calculatePathwayScore(input_g=g,
                                     another_pathway_node=one_name,
                                     WEIGHT_THRESHOLD=thresh,
                                     print_path = FALSE)
  }
  one_sums <- unlist(lapply(reses, function(x) if (length(x[[1]])>0) return(sum(1/x[[1]]))))
  one_names <- unlist(lapply(reses, function(x) if (length(x[[1]])>0) return(names(x[[1]])[1])))
  print(paste("threshold is", thresh))
  print(paste("The score sum is:", sum(one_sums)))
  thresh_score_sum <- c(thresh_score_sum, sum(one_sums))
}
neuro_thresh_score_sum <- thresh_score_sum

g <- nonneuro_g
thresh_score_sum <- numeric()
candidate_wp <- wpid2name[wpid2name$wpid != "WP4657",]$wpid
wp_list <- candidate_wp
for (thresh in threshold_range){
  #res_score <- list()
  # We'll only use pathways having a significant path for thresh==0.9.
  # Update and fix wp_list accordingly, starting from thresh==1.0
  if (thresh == 1.0){
    wp_list <- one_names
  }
  reses <- foreach (one_name=wp_list) %dopar% {
    one_res <- calculatePathwayScore(input_g=g,
                                     another_pathway_node=one_name,
                                     WEIGHT_THRESHOLD=thresh,
                                     print_path = FALSE)
  }
  one_sums <- unlist(lapply(reses, function(x) if (length(x[[1]])>0) return(sum(1/x[[1]]))))
  one_names <- unlist(lapply(reses, function(x) if (length(x[[1]])>0) return(names(x[[1]])[1])))
  print(paste("threshold is", thresh))
  print(paste("The score sum is:", sum(one_sums)))
  thresh_score_sum <- c(thresh_score_sum, sum(one_sums))
}
nonneuro_thresh_score_sum <- thresh_score_sum

# sensitivity_scores <- list(neuro_thresh_score_sum, nonneuro_thresh_score_sum)
# save(sensitivity_scores, file=file.path(RESULT_DIR, "sensitivity_scores.rda"))

splitVecSum <- function(inp_vec, n){
  # Split inp_vec into
  # equl lengths of n
  num_subvec <- length(inp_vec)/n
  vec_list <- lapply(seq(0, num_subvec-1), function(x) return(inp_vec[(x*n+1):((x+1)*n)]))
  return (vec_list)
}

neuro_thresh_score_list <- splitVecSum(neuro_thresh_score_sum, 167)
neuro_sensitivity_sum <- unlist(lapply(neuro_thresh_score_list, sum))

nonneuro_thresh_score_list <- splitVecSum(nonneuro_thresh_score_sum, 255)
nonneuro_sensitivity_sum <- unlist(lapply(nonneuro_thresh_score_list, sum))

compare_sensitivity <- data.frame(val=threshold_range,
                                  neuro_score=neuro_sensitivity_sum,
                                  nonneuro_score=nonneuro_sensitivity_sum)
write.csv(compare_sensitivity, file.path(RESULT_DIR, "compare_sensitivity.csv"), row.names = F)

stopImplicitCluster()
