# calculate_pval.R
# Calculate p value of the graph,
# By checking empirical test

library(foreach)
library(doParallel)
library(parallel)
numCores <- detectCores()-1
numCores
registerDoParallel(numCores)  # use multicore, set to the number of our cores


# Just test if it runs
foreach(icount(1000), .combine='+') %dopar% rnorm(4)
#??
#cl <- makeCluster(numCores) 
# registerDoParallel(cl)

# We match the number of genes, those with weight less than 0.8 and greater than 0.8
# 'g', "merged_edge", "merged_node" will be specified before running this script
# g <- nonneuro_g
# merged_edge <- nonneuro_merged_edge
# merged_node <- nonneuro_merged_node

############################################################
# First, using parallel computing, run with threshold 0.9, in order to obtain the list of 
# eligible pathways
reses <- foreach (one_name=candidate_wp) %dopar% {
  one_res <- calculatePathwayScore(input_g=g,
                                   another_pathway_node=one_name,
                                   WEIGHT_THRESHOLD=0.9,
                                   print_path = FALSE)
}
eligible_wp <- unlist(lapply(reses, function(x) if (length(x[[1]])>0) return(unique(names(x[[1]])))))
pathway_scores <- unlist(lapply(reses, function(x) if (length(x[[1]])>0) return(sum(1/x[[1]]))))
names(pathway_scores) <- eligible_wp
print(paste("Total number of eligible pathways:", length(eligible_wp)))
##############################################################

# Calculate p-values
one_candidate_node <- eligible_wp[1]

calculateOneRandScore <- function(one_cand){
  new_g <- createModifiedGraph(reference_edges=merged_edge,
                               reference_nodes=merged_node,
                               pathway_node_to_remove=one_cand,
                               to_match=TRUE)
  one_random_path_score <- calculatePathwayScore(input_g=new_g, another_pathway_node=one_cand)[[1]]
  one_random_score <- sum(1/one_random_path_score)
  return(one_random_score)
}

NUM_REPEAT <- 1000
sorted_wpid2name <- wpid2name
rownames(sorted_wpid2name) <- sorted_wpid2name$wpid
collected_pvals <- numeric()
count <- 0

# run next..
count <- length(collected_pvals)
for (one_candidate_node in eligible_wp[201:length(eligible_wp)]){
  count <- count + 1
  print(paste("We are pathway:", one_candidate_node, "with count:", count, "from total", length(eligible_wp), "nodes"))
  one_rand_scores_set <- foreach (one_name=1:NUM_REPEAT) %dopar% {
    calculateOneRandScore(one_cand=one_candidate_node)
  }
  one_pathway_score <- pathway_scores[one_candidate_node]
  one_pval <- sum(unlist(one_rand_scores_set) >= one_pathway_score)/NUM_REPEAT
  print(paste("p-value is:", one_pval))
  collected_pvals <- c(collected_pvals, one_pval)
}
res_pvals <- collected_pvals
names(res_pvals) <- eligible_wp[1:length(res_pvals)]
pval_df <- data.frame(wpid=names(res_pvals),
                      pvalue=res_pvals,
                      description=sorted_wpid2name[names(res_pvals),]$name,
                      stringsAsFactors = F)
# write.csv(pval_df, file.path(RESULT_DIR, paste("nonneuro_pval_df_", length(res_pvals), ".csv", sep="")), row.names = FALSE)
# nonneuro_pval_df <- pval_df

# Stop cluster
#stopCluster(cl)
stopImplicitCluster()

# input_g <- g
# one_pathway_node="WP4657"
# one_candidate_node="WP3972"
# another_pathway_node <- one_candidate_node
# print_path=TRUE
# WEIGHT_THRESHOLD=0.9
# 
# 
# # Try one pathway between two nodes
# 
# one_pathway_res <- calculatePathwayScore(input_g=g, another_pathway_node="WP4541")
# one_pathway_score <- sum(1/one_pathway_res[[1]])



############################################################
# First, run with threshold 0.9, in order to obtain the list of 
# eligible pathways
# eligible_wp <- character()
# candidate_wp <- wpid2name[wpid2name$wpid != "WP4657",]$wpid
# pathway_scores <- numeric()
# count <- 0
# for (one_name in candidate_wp){
#   if (count %% 100 == 0) print(paste("We are at:", count))
#   count <- count + 1
#   one_res <- calculatePathwayScore(input_g=g,
#                                    another_pathway_node=one_name,
#                                    WEIGHT_THRESHOLD=0.9,
#                                    print_path = FALSE)
#   # full_res_score is used to calculate overall distribution of path scores
#   print(paste("pathway score of", one_name))
#   print(one_res[[2]])
#   if (length(one_res[[1]])>0){
#     eligible_wp <- c(eligible_wp, unique(names(one_res[[1]])))
#     pathway_scores <- c(pathway_scores, sum(1/one_res[[1]]))
#   }
# }
# names(pathway_scores) <- eligible_wp
# print(paste("Total number of eligible pathways", length(eligible_wp)))
##############################################################

########################################
# one_candidate_node <- "WP3972"
# Using the above; eligible_wp', calculate the p-values
# NUM_REPEAT <- 1000
# collected_pvals <- numeric()
# all_random_scores <- list()
# count <- 0
# for (one_candidate_node in eligible_wp){
#   count <- count + 1
#   print(paste("We are pathway count", count, "from total", length(eligible_wp), "nodes"))
#   random_scores <- numeric()
#   for (idx in seq(1:NUM_REPEAT)){
#     if (idx %% 100 == 0) {
#       print(paste("Pathway:", one_candidate_node, "We are repetition at", idx))
#     }
#     new_g <- createModifiedGraph(reference_edges=merged_edge,
#                                  reference_nodes=merged_node,
#                                  pathway_node_to_remove=one_candidate_node,
#                                  to_match=TRUE)
#     one_random_path_score <- calculatePathwayScore(input_g=new_g, another_pathway_node=one_candidate_node)[[1]]
#     one_random_score <- sum(1/one_random_path_score)
#     random_scores <- c(random_scores, one_random_score)
#   }
#   one_pathway_score <- pathway_scores[one_candidate_node]
#   one_pval <- sum(random_scores >= one_pathway_score)/NUM_REPEAT
#   print(paste("p-value is:", one_pval))
#   collected_pvals <- c(collected_pvals, one_pval)
#   all_random_scores[[(length(all_random_scores)+1)]] <- random_scores
#   
#   if ((count %% 50 == 0) | (count==length(eligible_wp)) ){
#     sorted_wpid2name <- wpid2name
#     rownames(sorted_wpid2name) <- sorted_wpid2name$wpid
#     names(collected_pvals) <- eligible_wp
#     pval_df <- data.frame(wpid=names(collected_pvals),
#                           pvalue=collected_pvals,
#                           description=sorted_wpid2name[names(collected_pvals),]$name,
#                           stringsAsFactors = F)
#     write.csv(pval_df, file.path(RESULT_DIR, paste("nonneuro_pval_df_", count, ".csv", sep="")), row.names = FALSE)
#   }
# }


# neuro_collected_pvals <- collected_pvals
# sorted_wpid2name <- wpid2name
# rownames(sorted_wpid2name) <- sorted_wpid2name$wpid
# names(neuro_collected_pvals) <- eligible_wp
# neuro_pval_df <- data.frame(wpid=names(neuro_collected_pvals),
#                             pvalue=neuro_collected_pvals,
#                             description=sorted_wpid2name[names(neuro_collected_pvals),]$name,
#                             stringsAsFactors = F)
# write.csv(neuro_pval_df, file.path(RESULT_DIR, "neuro_pval_df.csv"), row.names = FALSE)

#sorted_wpid2name[sorted_wpid2name$wpid=="WP2583",]

########################################

# reference_edges = merged_edge
# reference_nodes = merged_node
# pathway_node_to_remove = another_pathway_node
# to_match=TRUE

