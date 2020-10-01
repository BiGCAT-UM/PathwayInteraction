# main.R
# Running the functions

library(dplyr)
library(igraph)
library(limma)
library(tidyr)

BASE_DIR <- "~/Desktop/BiGCAT/github/PathwayInteraction/"
DATA_DIR <- file.path(BASE_DIR, "data")
setwd(BASE_DIR)

# load limma result
limma_result <- read.csv(file=file.path(DATA_DIR, "limma_result_geomean.csv"))
