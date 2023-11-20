#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
# min length of mappings on a query chromosome to be included in plot
min_query_length <- as.numeric(args[2])

data <- read.delim(path)

# Function to filter data for a given chromosome
filter_data_for_chromosome <- function(data, chromosome) {
    data %>%
        filter(target.name == chromosome)
}

# Filter data for query chromosomes whose sum of mapping lengths in the data is greater than min_query_length
filter_data_for_query <- function(data, min_query_length) {
    data %>%
        group_by(query.name) %>%
        summarise(query_length = sum(query.end - query.start)) %>%
        filter(query_length > min_query_length) %>%
        select(query.name) %>%
        inner_join(data, by = "query.name")
}


# Unique list of chromosomes
chromosomes <- unique(data$target.name)

# Find the maximum length of target.name
max_target_length <- max(data$target.end)

# List to store individual ggplot objects
plots <- list()

for (chrom in chromosomes) {
  filtered_data <- filter_data_for_chromosome(data, chrom)
  filtered_data <- filter_data_for_query(filtered_data, min_query_length)
  p <- ggplot(filtered_data, aes(x = target.start, y = query.name, fill = identity)) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_minimal() +
    labs(title = paste("all vs.", chrom),
         x = "queries",
         y = "target",
         fill = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
        #scale_x_continuous(limits = c(0, max_target_length)) + # Adjust x-axis scale
        #scale_fill_gradient(low = "#B3B3B3", high = "#000000", limits = c(70, 100)) +
  #plots[[chrom]] <- p
  # replace any hashes in names with _
  name <- gsub("#", "_", chrom)
  # scale height by number of queries in dataset
  height <- length(unique(filtered_data$query.name)) * 0.15
  ggsave(paste0(path, "_", name, ".heat.pdf"), p, width = 10, height = height, limitsize = FALSE)
}

# Combine plots
#combined_plot <- wrap_plots(plots, ncol = 1)

# Save the combined plot
#ggsave(paste0(path, ".heat.svg"), combined_plot, width = 10, height = 500, limitsize = FALSE)
