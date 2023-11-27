#!/bin/bash

alignments=$1
genome=$2
prefix=$3
window_size=$4

# Extract relevant fields
# reordering fields to match bedtools intersect such that the target is the first 3 columns
# and the orientation and query is the next 4 columns
# and the score is the last column
# and removing the id:f: prefix from the query name
# order should be 6,8,9,1,3,4,5,13
awk -v OFS='\t' '{print $6,$8,$9,$1,$3,$4,$5,$13}' $alignments \
    | sed 's/id:f://g' | sed 's/gi:f://' > $alignments.bed

# Generate windows
bedtools makewindows -g <(grep ^$prefix $genome.fai) -w $window_size > windows.bed 

# Intersect 
bedtools intersect -a windows.bed -b $alignments.bed -wa -wb > $alignments.win.bed.unsrt

# Sort
sort -k1,1 -k2,2n $alignments.win.bed.unsrt > $alignments.win.bed
rm -f $alignments.win.bed.unsrt

# Add header and convert to tsv
# rewrite reference_ to target. and query_ to query.

( echo -e "target.name\ttarget.start\ttarget.end\tquery.name\tquery.start\tquery.end\tstrand\tidentity";
  cut -f 1-3,7- $alignments.win.bed ) \
    | sed 's/reference_/target./g' | sed 's/query_/query./g' \
    > $alignments.win.bed.tsv

