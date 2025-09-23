#!/usr/bin/env Rscript
options(scipen = 10000)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ape)
  library(ggplot2)
  library(scales)
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

# Process each tree
results <- tibble(file = tree_files) %>%
  mutate(
    basename = basename(file),
    region_full = gsub("^apes\\.", "", gsub("\\.nwk$", "", basename)),
    chr = gsub("_[0-9]+_[0-9]+$", "", region_full),
    start = as.numeric(gsub(".*_([0-9]+)_[0-9]+$", "\\1", region_full)),
    end = as.numeric(gsub(".*_([0-9]+)$", "\\1", region_full)),
    topology = map_chr(file, get_topology)
  ) %>%
  arrange(chr, start)

# Clean chromosome names and convert positions to Mbp
results <- results %>%
  mutate(
    chr_clean = gsub(".*#", "", chr),
    start_mbp = start / 1e6,
    end_mbp = end / 1e6,
    midpoint_mbp = (start_mbp + end_mbp) / 2
  )

# Create factor levels for chromosomes
chrom_levels <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
results$chr_clean <- factor(results$chr_clean, levels = chrom_levels)

# Summary statistics
cat("\n=== Topology Summary ===\n")
topology_counts <- table(results$topology)
print(topology_counts)

cat("\n=== Percentage of each topology ===\n")
print(round(prop.table(topology_counts) * 100, 2))

# Create the wide format plot
p_ils_wide <- ggplot(results) +
  # Add rectangles for each topology region
  geom_rect(aes(xmin = start_mbp, xmax = end_mbp, 
                ymin = 0, ymax = 1, 
                fill = topology),
            alpha = 0.9) +
  
  # Facet by chromosome with free x scales
  facet_grid(rows = vars(chr_clean), scales = "free_x", switch = "y") +
  
  # Color scheme matching your original
  scale_fill_manual(
    values = c(
      "species_tree" = "#00BA38",  # Green
      "ILS_HB" = "#F8766D",        # Red
      "ILS_HC" = "#619CFF",        # Blue
      "equal" = "gray80",
      "NA" = "white"
    ),
    name = "Topology"
  ) +
  
  # Set x-axis 
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 5),
    expand = c(0.01, 0)
  ) +
  
  # Fixed y-axis since we're just showing presence/absence
  scale_y_continuous(
    limits = c(0, 1),
    breaks = NULL
  ) +
  
  labs(
    title = "Incomplete Lineage Sorting (ILS) Patterns Across Primate Genomes",
    subtitle = "Human-Chimpanzee-Bonobo topology analysis using IMPG",
    x = "Position (Mbp)",
    y = NULL
  ) +
  
  theme_minimal() +
  theme(
    # Facet strip styling
    strip.text.y.left = element_text(size = 9, angle = 0),
    strip.background = element_rect(fill = "gray95", color = NA),
    strip.placement = "outside",
    
    # Axis styling
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 8),
    axis.title = element_text(size = 11),
    axis.ticks.y = element_blank(),
    
    # Title styling
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    
    # Panel and grid styling
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.major.y = element_blank(),
    panel.spacing = unit(0.1, "lines"),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.3),
    
    # Legend styling
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(1.5, "cm")
  )

# Save the plot
ggsave("ILS_topology_by_chromosome.pdf", p_ils_wide, 
       width = 14, height = 10, dpi = 300)
ggsave("ILS_topology_by_chromosome.png", p_ils_wide, 
       width = 14, height = 10, dpi = 300)

cat("\n=== Plots saved ===\n")
cat("- ILS_topology_by_chromosome.pdf\n")
cat("- ILS_topology_by_chromosome.png\n")

# Create summary by chromosome
chr_summary <- results %>%
  group_by(chr_clean) %>%
  summarise(
    total_regions = n(),
    species_tree = sum(topology == "species_tree", na.rm = TRUE),
    ILS_HB = sum(topology == "ILS_HB", na.rm = TRUE),
    ILS_HC = sum(topology == "ILS_HC", na.rm = TRUE),
    pct_ILS = round((ILS_HB + ILS_HC) / total_regions * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(chr_clean)

cat("\n=== ILS Summary by Chromosome ===\n")
print(as.data.frame(chr_summary))

# Write results
write_tsv(results %>% select(chr_clean, start, end, start_mbp, end_mbp, topology), 
          "ILS_topology_results.tsv")
write_tsv(chr_summary, "ILS_summary_by_chromosome.tsv")

cat("\n=== Files written ===\n")
cat("- ILS_topology_results.tsv\n")
cat("- ILS_summary_by_chromosome.tsv\n")
