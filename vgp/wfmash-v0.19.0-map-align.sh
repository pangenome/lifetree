#!/bin/bash

# wfmash.sh: run wfmash on a pair of sequences
# for version v0.19.0

# break on error
set -euox pipefail

seqs=$1
jobs=$2
id=$3
out=$4

cpus=$SLURM_CPUS_PER_TASK
# set to max of parameter defined and available cpus
cpus_for_mapping=$(($cpus > $5 ? $5 : $cpus))
cpus_for_aligning=$(($cpus > $6 ? $6 : $cpus))
echo "cpus for mapping: $cpus_for_mapping"
echo "cpus for aligning: $cpus_for_aligning"

time=/usr/bin/time
fmt="%C\n%Us user %Ss system %P cpu %es total %MKb max memory"
timer=$(which time)

ls -l $seqs
mkdir -p $out

spec=$(head -$id $jobs | tail -1)
target=$(echo $spec | awk '{print $1}')
query=$(echo $spec | awk '{print $2}')
mm3idx=$SCRATCH/idxs/$target.mm3

echo "wfmash mapping $id = $target against $query on $(hostname) with $cpus CPUs at $(date +%s) / $(date) with $(which wfmash)"
wfmash-v0.19.0-0-gc95b6a3 --version

map_out=$out/$id.$target.map

( $timer -f "$fmt" \
   wfmash-v0.19.0-0-gc95b6a3 \
        -t $cpus_for_mapping \
        -m \
        --mm-index $mm3idx \
        -T $target \
        -Q $query \
        -Y '#' \
        -n 1 \
        -k 19 \
        -p 70 \
        -s 5k \
        -c 30k \
        -P 1M \
        $seqs \
        >$map_out.paf ) 2>$map_out.log
touch $map_out.ok

echo "wfmash aligning $id = $target against $query on $(hostname) with $cpus CPUs at $(date +%s) / $(date)"

aln_out=$out/$id.$target.aln
( $timer -f "$fmt" \
    wfmash-v0.19.0-0-gc95b6a3 \
        -t $cpus_for_aligning \
        -i $map_out.paf \
        $seqs \
        >$aln_out.paf ) 2>$aln_out.log
touch $aln_out.ok

echo "zipping $map_out.paf and $aln_out.paf"
pigz -f $map_out.paf $aln_out.paf

touch $id.$target.done
