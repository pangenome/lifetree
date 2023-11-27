# Alignment Heatmap Scripts

This directory contains two scripts to generate a heatmap visualization of sequence alignment coverage across chromosomes.

## overlaps.sh

This Bash script takes a sequence alignment file, reference genome, prefix, and window size as input. It does the following:

- Extracts and reorders relevant fields from the alignment file into BED format
- Generates fixed-size windows across the genome using `bedtools makewindows` 
- Intersects the alignments with the windows using `bedtools intersect`
- Sorts the output by chromosome and position
- Adds a header and converts to TSV format

The output is a file named `${prefix}.win.bed.tsv` containing the intersection of the alignments and windows, ready for downstream analysis and plotting.

## plot.R

This R script takes the output TSV file from `overlaps.sh` and generates a heatmap plot for each chromosome. It does the following:

- Reads in the TSV data
- Filters the data by chromosome 
- Plots a tile heatmap for each chromosome with position on x-axis, query sequence on y-axis, and identity as fill
- Combines the chromosome plots into one tall plot
- Saves the combined plot as a PDF  

The output is a PDF file named `${prefix}.heat.pdf` containing the heatmap visualization.

## Usage

The typical usage would be:

```
# Generate windows
./overlaps.sh z.4.paf ~/yeast/cerevisiae.fa.gz S288C 10000

# Plot
./plot.R z.4.paf.win.bed.tsv
```

This will produce a heatmap PDF visualization of the provided alignments file across 10kb windows of the genome.
The scripts are customizable for different data sources and parameters.
