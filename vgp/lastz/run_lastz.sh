#!/bin/bash

# Arguments from submission script
query_fna="$1"      # path to query FASTA file
reference_fna="$2"  # path to reference FASTA file
dir_output="$3"     # output directory
LASTZ="$4"          # lastz binary

# Get the basename for this query_fna file
name=$(basename "$query_fna" .fna.gz)
echo "Processing $name"
echo "Query FNA file: $query_fna"
echo "Output directory: $dir_output"
echo "Reference FNA file: $reference_fna"

# Create temporary directory on scratch
temp_dir="/scratch/${SLURM_JOB_ID}_${name}"
mkdir -p "$temp_dir"

# Create output directories if they don't exist
mkdir -p "$dir_output/logs/$name"
mkdir -p "$dir_output"


# Export variables for parallel
export dir_output query_fna name temp_dir LASTZ reference_fna

# Generate all combinations and run in parallel
cut -f 1 "${reference_fna}.fai" | while read target; do
    cut -f 1 "${query_fna}.fai" | while read query; do
        echo "$target $query"
    done
done | parallel -j 48 --colsep ' ' '
    # Prepare sequences
    target={1}
    query={2}
    
    # Use unique temp files with PID and job slot number
    temp_target="${temp_dir}/${target}.${PARALLEL_PID}.${PARALLEL_JOBSLOT}.fna"
    temp_query="${temp_dir}/${query}.${PARALLEL_PID}.${PARALLEL_JOBSLOT}.fna"
    
    # Extract sequences
    samtools faidx "$reference_fna" "$target" > "$temp_target"
    samtools faidx "$query_fna" "$query" > "$temp_query"
    
    # Run alignment
    \time -v "$LASTZ" "$temp_target" "$temp_query" --format=PAF:wfmash \
        > "${temp_dir}/${name}.aln.${query}-vs-${target}.paf" \
        2> "$dir_output/logs/${name}/${name}.${query}-vs-${target}.log"
    
    # Clean up temp files immediately after each job
    rm -f "$temp_target" "$temp_query"
'

# Merge all PAF files
echo "Merging PAF files for $name"
find "${temp_dir}" -name "${name}.aln.*.paf" -print0 | \
    xargs -0 cat > "$dir_output/${name}.aln.paf"

# Final cleanup
echo "Cleaning up temporary directory"
rm -rf "$temp_dir"

echo "Completed processing $name"
