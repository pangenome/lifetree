#!/bin/bash

refs=$1
panfasta=$2
outbase=$3

echo "writing commands for indexing to commands.txt"

for ref in $(cat $refs);
do
    echo wfmash-v0.21.0-0-g4521c10 -t 48 -T $ref --mm-index $outbase/$ref.mm3 -n 1 -k 15 -p 70 -s 1k -l 3k -c 30k --create-index-only $panfasta '2>'$outbase/$ref.mm3.log
done >commands.txt
