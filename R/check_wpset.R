# check_wpset.R
# check the list of wp's
# in the wikipathway and reactome folders

# RUN THIS SCRIPT AFTER LINE 110 OF main.R
ALL_DIR <- file.path(WP_DIR, "wp_matched")
all_wp_files <- list.files(ALL_DIR)

# getSplit <- function(x){
#   val <- unlist(strsplit(x, "_"))
#   return (tail(val, 2)[1])
# }
all_wpid <- sapply(all_wp_files, getSplit)
names(all_wpid) <- NULL

common_wpid <- unique(wpid2name$wpid)
wps_to_remove <- setdiff(all_wpid, common_wpid)

for (one_wp_file in all_wp_files) {
  one_wp_file_id <- getSplit(one_wp_file)
  if (one_wp_file_id %in% wps_to_remove){
    file.remove(file.path(ALL_DIR, one_wp_file))
  }
}


# Removed!
length(all_wp_files)
# [1] 1150
length(list.files(ALL_DIR))
# [1] 1136

