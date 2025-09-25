#!/usr/bin/env Rscript
options(scipen = 10000)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ape)
  library(ggplot2)
  library(scales)
  library(phangorn)
  library(parallel)  # For parallel processing
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript topology_stability_analysis_improved.R <path_to_nwk_files> [reference_tree.nwk]
       path_to_nwk_files: Directory containing .nwk files to process
       reference_tree.nwk: Optional reference tree file (default: trees.nwk in current directory)")
}

target_path <- args[1]
reference_tree_file <- if (length(args) >= 2) args[2] else "trees.nwk"

# Check if target path exists
if (!dir.exists(target_path)) {
  stop(paste("Error: Directory", target_path, "does not exist"))
}

message(paste("Processing .nwk files in:", target_path))
message(paste("Using reference tree:", reference_tree_file))

# Function to calculate topological distance between two trees
calculate_topological_distance <- function(tree1, tree2) {
  tryCatch({
    # Convert to phylo objects if needed
    if (is.character(tree1)) tree1 <- read.tree(text = tree1)
    if (is.character(tree2)) tree2 <- read.tree(text = tree2)
    
    # Ensure both trees have the same taxa
    taxa1 <- tree1$tip.label
    taxa2 <- tree2$tip.label
    
    # Find common taxa
    common_taxa <- intersect(taxa1, taxa2)
    if (length(common_taxa) < 3) {
      return(NA)  # Need at least 3 taxa for meaningful comparison
    }
    
    # Drop uncommon taxa
    if (length(taxa1) > length(common_taxa)) {
      tree1 <- drop.tip(tree1, setdiff(taxa1, common_taxa))
    }
    if (length(taxa2) > length(common_taxa)) {
      tree2 <- drop.tip(tree2, setdiff(taxa2, common_taxa))
    }
    
    # Calculate Robinson-Foulds distance (normalized)
    rf_dist <- RF.dist(tree1, tree2, normalize = TRUE)
    
    return(rf_dist)
    
  }, error = function(e) {
    message(paste("Error calculating distance:", e$message))
    return(NA)
  })
}

# Function to calculate tree stability (lower distance = higher stability)
calculate_stability_score <- function(rf_distance) {
  if (is.na(rf_distance)) return(NA)
  # Convert RF distance to stability score (0 = identical, 1 = maximally different)
  # Stability score: 1 - normalized_rf_distance (1 = most stable, 0 = least stable)
  stability <- 1 - rf_distance
  return(stability)
}

# Improved function to extract genomic coordinates from filename
extract_coordinates <- function(filename) {
  basename <- basename(filename)
  
  # Try multiple patterns for chromosome/contig extraction
  chr_patterns <- c(
    "CP[0-9]+\\.[0-9]+",     # CP068277.2
    "chr[0-9XYM]+",          # chr1, chrX, etc
    "Chr[0-9XYM]+",          # Chr1, ChrX, etc
    "scaffold_[0-9]+",       # scaffold_1
    "#([^#]+)_[0-9]+_[0-9]+" # Extract between last # and first _
  )
  
  chr_contig <- NA
  for (pattern in chr_patterns) {
    chr_match <- str_extract(basename, pattern)
    if (!is.na(chr_match)) {
      chr_contig <- chr_match
      break
    }
  }
  
  # If still no match, try to extract between # symbols
  if (is.na(chr_contig)) {
    # Pattern: anything#number#contig_start_end.nwk
    if (grepl("#[^#]+#([^_]+)", basename)) {
      chr_contig <- str_extract(basename, "#[^#]+#([^_]+)") %>%
        str_extract("[^#]+$")
    }
  }
  
  # Extract positions - look for pattern _number_number at the end
  position_match <- str_extract(basename, "_([0-9]+)_([0-9]+)\\.[^.]+$")
  if (!is.na(position_match)) {
    numbers <- str_extract_all(position_match, "[0-9]+")[[1]]
    if (length(numbers) >= 2) {
      start <- as.numeric(numbers[1])
      end <- as.numeric(numbers[2])
    } else {
      start <- NA
      end <- NA
    }
  } else {
    start <- NA
    end <- NA
  }
  
  return(list(chr_contig = chr_contig, start = start, end = end))
}

# Process a single file
process_single_file <- function(file_path, reference_tree) {
  coords <- extract_coordinates(file_path)
  
  segment_tree <- tryCatch({
    read.tree(file_path)
  }, error = function(e) {
    return(NULL)
  })
  
  if (!is.null(segment_tree)) {
    rf_dist <- calculate_topological_distance(reference_tree, segment_tree)
    stability <- calculate_stability_score(rf_dist)
  } else {
    rf_dist <- NA
    stability <- NA
  }
  
  return(data.frame(
    file = file_path,
    chr_contig = coords$chr_contig,
    start = coords$start,
    end = coords$end,
    rf_distance = rf_dist,
    stability_score = stability,
    stringsAsFactors = FALSE
  ))
}

# Main execution
message("Reading reference tree...")
reference_tree <- read.tree(reference_tree_file)

# Check reference tree
message(paste("Reference tree has", length(reference_tree$tip.label), "taxa"))
message("Sample taxa names:")
print(head(reference_tree$tip.label, 10))

# Get all newick files in the target directory
message("Scanning for segment tree files...")
tree_files <- list.files(target_path, pattern = "\\.nwk$", recursive = TRUE, full.names = TRUE)
# Exclude the reference tree if it's in the same directory
tree_files <- tree_files[!grepl(basename(reference_tree_file), tree_files)]

message(paste("Found", length(tree_files), "segment tree files"))

# Test coordinate extraction on a few files
message("\nTesting coordinate extraction on first 5 files:")
for (i in 1:min(5, length(tree_files))) {
  coords <- extract_coordinates(tree_files[i])
  message(sprintf("  File: %s", basename(tree_files[i])))
  message(sprintf("    Chr: %s, Start: %s, End: %s", 
                  coords$chr_contig, coords$start, coords$end))
}

# Determine number of cores to use
n_cores <- min(detectCores() - 1, 8)  # Leave one core free, max 8 cores
message(paste("\nUsing", n_cores, "cores for parallel processing"))

# Process files in parallel batches
batch_size <- 1000
all_results <- list()

for (batch_start in seq(1, length(tree_files), by = batch_size)) {
  batch_end <- min(batch_start + batch_size - 1, length(tree_files))
  batch_files <- tree_files[batch_start:batch_end]
  
  message(paste("\nProcessing batch", ceiling(batch_start/batch_size), 
                "of", ceiling(length(tree_files)/batch_size),
                "- files", batch_start, "to", batch_end))
  
  # Process batch in parallel
  if (n_cores > 1) {
    cl <- makeCluster(n_cores)
    clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(ape)
        library(phangorn)
        library(stringr)
      })
    })
    clusterExport(cl, c("calculate_topological_distance", 
                        "calculate_stability_score",
                        "extract_coordinates",
                        "process_single_file",
                        "reference_tree"))
    
    batch_results <- parLapply(cl, batch_files, process_single_file, reference_tree)
    stopCluster(cl)
  } else {
    batch_results <- lapply(batch_files, process_single_file, reference_tree)
  }
  
  all_results <- append(all_results, batch_results)
}

# Combine results
message("\nCombining results...")
results <- bind_rows(all_results)

# Filter and process
results <- results %>%
  filter(!is.na(start) & !is.na(end) & !is.na(stability_score)) %>%
  mutate(
    start_mbp = start / 1e6,
    end_mbp = end / 1e6,
    midpoint_mbp = (start_mbp + end_mbp) / 2,
    # Create stability categories for visualization
    stability_category = case_when(
      stability_score >= 0.8 ~ "Highly Stable",
      stability_score >= 0.6 ~ "Stable", 
      stability_score >= 0.4 ~ "Moderately Stable",
      stability_score >= 0.2 ~ "Unstable",
      TRUE ~ "Highly Unstable"
    )
  ) %>%
  arrange(chr_contig, start)

message(paste("Successfully analyzed", nrow(results), "segments"))

# Summary statistics
cat("\n=== Topology Stability Summary ===\n")
stability_summary <- results %>%
  summarise(
    total_segments = n(),
    mean_stability = round(mean(stability_score, na.rm = TRUE), 3),
    median_stability = round(median(stability_score, na.rm = TRUE), 3),
    min_stability = round(min(stability_score, na.rm = TRUE), 3),
    max_stability = round(max(stability_score, na.rm = TRUE), 3),
    highly_stable = sum(stability_score >= 0.8, na.rm = TRUE),
    highly_unstable = sum(stability_score < 0.2, na.rm = TRUE)
  )

print(stability_summary)

# Category distribution
cat("\n=== Stability Category Distribution ===\n")
category_counts <- table(results$stability_category)
print(category_counts)
print(round(prop.table(category_counts) * 100, 1))

# Create chromosome-level visualization
message("\nCreating visualization...")
p_stability <- results %>%
  filter(!is.na(chr_contig)) %>%
  mutate(chr_contig = factor(chr_contig)) %>%
  ggplot() +
  geom_rect(aes(xmin = start_mbp, xmax = end_mbp, 
                ymin = as.numeric(chr_contig) - 0.4, ymax = as.numeric(chr_contig) + 0.4, 
                fill = stability_score),
            alpha = 0.9) +
  
  # Color scale: deep blue for stable (high score), deep red for unstable (low score)
  scale_fill_gradient2(
    low = "#8B0000",      # Deep red for unstable
    mid = "#FFFF99",      # Light yellow for moderate
    high = "#00008B",     # Deep blue for stable
    midpoint = 0.5,
    name = "Stability\nScore",
    labels = function(x) sprintf("%.2f", x),
    limits = c(0, 1)
  ) +
  
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 6),
    expand = c(0.01, 0)
  ) +
  
  scale_y_continuous(
    breaks = function(x) {
      chr_levels <- unique(results$chr_contig[!is.na(results$chr_contig)])
      chr_levels <- sort(chr_levels)
      seq_along(chr_levels)
    },
    labels = function(x) {
      chr_levels <- unique(results$chr_contig[!is.na(results$chr_contig)])
      chr_levels <- sort(chr_levels)
      chr_levels[x]
    },
    expand = c(0.02, 0)
  ) +
  
  labs(
    title = "Chromosomal Segments Topological Stability Analysis",
    subtitle = paste("Analysis of", nrow(results), "genomic segments compared to reference phylogeny"),
    x = "Position (Mbp)",
    y = "Chromosome",
    caption = "Deep Blue = Conserved topology | Deep Red = Variable topology"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.title = element_text(size = 11),
    axis.title.y = element_text(angle = 90, vjust = 0.5),
    
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    plot.caption = element_text(size = 9, color = "gray60"),
    
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.2),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.3),
    
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.width = unit(0.8, "cm"),
    legend.key.height = unit(1.5, "cm")
  )

# Save the plot
message("Saving plots...")
# Calculate dynamic height based on number of chromosomes
n_chr <- length(unique(results$chr_contig[!is.na(results$chr_contig)]))
plot_height <- max(6, n_chr * 0.5 + 3)  # 0.5 inches per chromosome plus margins
ggsave("topology_stability_analysis.pdf", p_stability, 
       width = 12, height = plot_height, dpi = 300)
ggsave("topology_stability_analysis.png", p_stability, 
       width = 12, height = plot_height, dpi = 300)

# Create a summary heatmap by chromosome
if (length(unique(results$chr_contig)) > 1) {
  chr_summary <- results %>%
    group_by(chr_contig) %>%
    summarise(
      n_segments = n(),
      mean_stability = mean(stability_score, na.rm = TRUE),
      median_stability = median(stability_score, na.rm = TRUE),
      min_stability = min(stability_score, na.rm = TRUE),
      max_stability = max(stability_score, na.rm = TRUE),
      stable_proportion = sum(stability_score >= 0.6, na.rm = TRUE) / n(),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_stability))
  
  p_chr_summary <- chr_summary %>%
    ggplot(aes(x = 1, y = reorder(chr_contig, mean_stability))) +
    geom_tile(aes(fill = mean_stability), color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.3f", mean_stability)), 
              color = "white", fontface = "bold") +
    scale_fill_gradient2(
      low = "#8B0000", mid = "#FFFF99", high = "#00008B",
      midpoint = 0.5,
      name = "Mean\nStability",
      limits = c(0, 1)
    ) +
    labs(
      title = "Chromosome-level Topology Stability",
      subtitle = "Mean stability score per chromosome/contig",
      x = NULL,
      y = "Chromosome/Contig"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      legend.position = "right"
    )
  
  ggsave("chromosome_stability_summary.pdf", p_chr_summary, 
         width = 10, height = max(6, length(unique(results$chr_contig)) * 0.3), 
         dpi = 300)
  ggsave("chromosome_stability_summary.png", p_chr_summary, 
         width = 10, height = max(6, length(unique(results$chr_contig)) * 0.3), 
         dpi = 300)
  
  cat("\n=== Chromosome Summary ===\n")
  print(as.data.frame(chr_summary), row.names = FALSE)
}

# Write results to files
message("\nWriting output files...")
write_tsv(results %>% select(file, chr_contig, start, end, start_mbp, end_mbp, 
                            rf_distance, stability_score, stability_category), 
          "topology_stability_results.tsv")

if (exists("chr_summary")) {
  write_tsv(chr_summary, "chromosome_stability_summary.tsv")
}

# Report problematic files
n_failed <- length(tree_files) - nrow(results)
if (n_failed > 0) {
  cat("\n=== Files with issues ===\n")
  cat(paste(n_failed, "files could not be processed (missing coordinates or read errors)\n"))
}

cat("\n=== Files Generated ===\n")
cat("- topology_stability_analysis.pdf/png\n")
cat("- topology_stability_results.tsv\n")
if (exists("chr_summary")) {
  cat("- chromosome_stability_summary.pdf/png\n")
  cat("- chromosome_stability_summary.tsv\n")
}

message("\nAnalysis completed successfully!")
