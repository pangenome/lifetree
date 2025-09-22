#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(ape)
})

# Function to determine topology for 3-taxon tree
get_topology <- function(tree_file) {
  tree <- read.tree(tree_file)
  
  # Get the cophenetic distance matrix
  dist_mat <- cophenetic(tree)
  
  # Ensure we have the expected taxa
  taxa <- rownames(dist_mat)
  human <- grep("chm13", taxa, value = TRUE)[1]
  chimp <- grep("mPanTro3", taxa, value = TRUE)[1]
  bonobo <- grep("mPanPan1", taxa, value = TRUE)[1]
  
  if(is.na(human) || is.na(chimp) || is.na(bonobo)) {
    return(NA)
  }
  
  # Get pairwise distances
  d_hc <- dist_mat[human, chimp]
  d_hb <- dist_mat[human, bonobo]
  d_cb <- dist_mat[chimp, bonobo]
  
  # Determine topology based on which pair has smallest distance
  if(d_cb < d_hc && d_cb < d_hb) {
    return("species_tree")  # ((chimp,bonobo),human)
  } else if(d_hb < d_hc && d_hb < d_cb) {
    return("ILS_HB")  # ((human,bonobo),chimp)
  } else if(d_hc < d_hb && d_hc < d_cb) {
    return("ILS_HC")  # ((human,chimp),bonobo)
  } else {
    return("equal")  # All distances equal
  }
}

# Get all newick files
tree_files <- list.files("apes-chm13", pattern = "\\.nwk$", full.names = TRUE)

# Process each tree with corrected parsing
results <- tibble(file = tree_files) %>%
  mutate(
    basename = basename(file),
    # Extract the full region string after "apes."
    region_full = gsub("^apes\\.", "", gsub("\\.nwk$", "", basename)),
    # Extract chromosome (everything before the last two underscores)
    chr = gsub("_[0-9]+_[0-9]+$", "", region_full),
    # Extract start position
    start = as.numeric(gsub(".*_([0-9]+)_[0-9]+$", "\\1", region_full)),
    # Extract end position  
    end = as.numeric(gsub(".*_([0-9]+)$", "\\1", region_full)),
    # Get topology
    topology = map_chr(file, get_topology)
  ) %>%
  arrange(chr, start)

# Summary statistics
cat("\n=== Topology Summary ===\n")
topology_counts <- table(results$topology)
print(topology_counts)

cat("\n=== Percentage of each topology ===\n")
print(round(prop.table(topology_counts) * 100, 2))

# Write detailed results
write_tsv(results %>% select(chr, start, end, topology, file), 
          "tree_topology_comparison.tsv")

# Find ILS regions - create proper BED format
ils_regions <- results %>%
  filter(topology %in% c("ILS_HB", "ILS_HC")) %>%
  mutate(
    # Create clean chromosome name for BED
    chr_clean = gsub(".*#", "", chr)  # Remove "chm13#1#" prefix to get just "chr1"
  ) %>%
  select(chr_clean, start, end, topology) %>%
  rename(chr = chr_clean)

write_tsv(ils_regions, "ILS_regions.bed", col_names = FALSE)

cat("\nResults written to:\n")
cat("- tree_topology_comparison.tsv (all regions with headers)\n")
cat("- ILS_regions.bed (regions showing ILS, BED format)\n")

# Create visualization if we have valid coordinates
if(all(!is.na(results$start))) {
  library(ggplot2)
  
  # Clean chromosome names for plotting
  plot_data <- results %>%
    mutate(chr_clean = gsub(".*#", "", chr))
  
  p <- ggplot(plot_data, aes(x = start, y = 1, fill = topology)) +
    geom_tile(aes(width = end - start, height = 0.5)) +
    facet_wrap(~chr_clean, scales = "free_x", ncol = 1) +
    scale_fill_manual(values = c("species_tree" = "green", 
                                "ILS_HB" = "red", 
                                "ILS_HC" = "blue",
                                "equal" = "gray",
                                "NA" = "white")) +
    theme_minimal() +
    labs(title = "Tree topology along chromosomes",
         x = "Position", y = "", fill = "Topology") +
    theme(axis.text.y = element_blank())
  
  ggsave("topology_changes.pdf", p, width = 12, height = 8)
  cat("- topology_changes.pdf (visualization)\n")
}

# Print first few ILS regions to verify
cat("\n=== First few ILS regions ===\n")
print(head(ils_regions))