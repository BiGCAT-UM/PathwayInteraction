# calculate_pval.R
# Calculate p value of the graph,
# By using empirical test

library(foreach)
library(doParallel)
library(parallel)
numCores <- detectCores()-1
numCores
registerDoParallel(numCores)  # use multicore, set to the number of our cores

# Just test if it runs
foreach(icount(1000), .combine='+') %dopar% rnorm(4)

############################################################
# First, using parallel computing, run with threshold 0.9, in order to obtain the list of 
# eligible pathways
candidate_wp <- wpid2name$wpid
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
  ### Feb1, 2021, pval formula CHANGED to  (r+1)/(n+1). For comparison, '>=' is correct
  ### source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC379178/
  print(paste("r = ", sum(unlist(one_rand_scores_set) >= one_pathway_score)))
  one_pval <- (sum(unlist(one_rand_scores_set) >= one_pathway_score)+1)/(NUM_REPEAT+1)
  print(paste("p-value is:", one_pval))
  collected_pvals <- c(collected_pvals, one_pval)
}

#collected_pvals <- c(collected_pvals[1:180], collected_pvals[182:length(collected_pvals)])
res_pvals <- collected_pvals
names(res_pvals) <- eligible_wp[1:length(res_pvals)]
pval_df <- data.frame(wpid=names(res_pvals),
                      pvalue=res_pvals,
                      description=sorted_wpid2name[names(res_pvals),]$name,
                      stringsAsFactors = F)
# write.csv(pval_df, file.path(RESULT_DIR, paste("neuro_pval_df_", length(res_pvals), "_2021.csv", sep="")), row.names = FALSE)
# write.csv(pval_df, file.path(RESULT_DIR, paste("nonneuro_pval_df_", length(res_pvals), "_2021.csv", sep="")), row.names = FALSE)

# Stop cluster
stopImplicitCluster()
