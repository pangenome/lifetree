#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)
tsv_file <- args[1]
output_base <- gsub(".similarity.tsv", "", tsv_file)

tryCatch({
  # Read and process distance matrix
  sparse_matrix_df <- read_tsv(tsv_file, show_col_types = FALSE)
  
  # Check if we have enough data
  n_groups <- length(unique(c(sparse_matrix_df$group.a, sparse_matrix_df$group.b)))
  if (n_groups < 2) {
    cat("Skipping", tsv_file, ": only", n_groups, "group(s) found\n")
    quit(status = 0)
  }
  
  jaccard_dist_df <- sparse_matrix_df %>%
    arrange(group.a, group.b) %>%
    select(group.a, group.b, jaccard.distance) %>%
    pivot_wider(names_from = group.b, values_from = jaccard.distance) %>%
    column_to_rownames(var = "group.a")
  
  # Create distance object
  jaccard_dist <- as.dist(jaccard_dist_df)
  
  # Check for NA/NaN/Inf values
  if (any(!is.finite(jaccard_dist))) {
    cat("Skipping", tsv_file, ": contains NA/NaN/Inf values\n")
    quit(status = 0)
  }
  
  # Create tree
  tree <- as.phylo(hclust(jaccard_dist))
  
  # Write Newick file
  write.tree(tree, file = paste0(output_base, ".nwk"))
  cat("Created tree for", tsv_file, "\n")
  
}, error = function(e) {
  cat("Error processing", tsv_file, ":", e$message, "\n")
})
