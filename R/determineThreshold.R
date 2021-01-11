# determineThreshold.R
# by changing the change of scores, of selected 55 wikipathways,
# choose the optimal threshold value other than 0.9

# First, we need final_df from visualizePathwayInteraction.R

wiki_map <- read.csv(file.path(DATA_DIR, "wiki_map.csv"), stringsAsFactors = F)

all_wiki_score_list <- list()
for (one_thresh in seq(0.9, 3.0, by=0.1)){
  one_wiki_score <- numeric()
  print("****************************")
  print(paste("We are checking threshold=", one_thresh))
  print("****************************")
  for (one_name in final_df$name){
    print(paste("We are now with ", one_name))
    one_res <- calculatePathwayScore(input_g=g,
                                     another_pathway_node=one_name,
                                     WEIGHT_THRESHOLD = one_thresh)
    wiki_score <- c(wiki_score, sum(1/one_res[[1]]))
  }
  all_wiki_score_list[[(length(all_wiki_score_list)+1)]] <- wiki_score
}

ss <- sapply(all_wiki_score_list, sum)
thresh_ss <- seq(0.9, 2.1, by=0.1)
plot(x=thresh_ss, y=ss, xlab="Score threshold", ylab="Sum of total path scores")
# [1]  1324.536  1647.767  2000.666  2384.981  2867.830  3419.383  4139.040  5629.248  7828.699 11205.279 17837.989 27101.223
# [13] 36451.299
slopes = numeric()
for (idx in 2:length(ss)){
  one_slope = (ss[idx]-ss[(idx-1)])/0.1
  slopes <- c(slopes, one_slope)
}

xbreaks <- numeric()
for (val in head(thresh_ss,-1)){
  xbreaks <- c(xbreaks, paste(val, val+0.1, sep="-"))
}

plot(slopes)
# > slopes
# [1]  3232.310  3528.988  3843.149  4828.489  5515.527  7196.578 14902.077 21994.505 33765.802 66327.099 92632.343 93500.755


