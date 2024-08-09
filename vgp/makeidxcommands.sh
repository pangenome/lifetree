#!/bin/bash

refs=$1
panfasta=$2
outbase=$3

echo "writing commands for indexing to commands.txt"

for ref in $(cat $refs);
do
    echo wfmash -t 48 -T $ref --mm-index $outbase/$ref.mm3 -n 1 -k 19 -p 70 -s 5k -c 30k --create-index-only $panfasta '2>'$outbase/$ref.mm3.log
done >commands.txt
