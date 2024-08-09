#!/bin/bash

target_genome=$1
panfasta=$2
joblist=$3
outbase=$4
script=$5
map_threads=$6
aln_threads=$7

a=$(cat $joblist | awk '$1 ~ /'$target_genome'/ { print NR, $0 }' | head -1 | cut -f 1 -d\ )
b=$(cat $joblist | awk '$1 ~ /'$target_genome'/ { print NR, $0 }' | tail -1 | cut -f 1 -d\ )

echo "writing commands.txt for job ids $a to $b"

for id in $(seq $a $b);
do
    echo $script $panfasta $joblist $id $outbase $map_threads $aln_threads
done >commands.txt
