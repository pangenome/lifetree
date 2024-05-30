#!/bin/bash
#SBATCH --job-name=wfmash
#SBATCH --nodes=1
#SBATCH --time=5-00:00:00
#SBATCH --partition=skx

# map vgp+477 to itself at 70% identity 5kb segments
# skipping self alignments
# using the one-to-one best mapping filter

eval "$(micromamba shell hook --shell bash)"
micromamba activate env

wfmash=$(which wfmash)
time=/usr/bin/time

# Determine number of CPUs to use
CPUS=$(($SLURM_NTASKS_PER_NODE))
base=$WORK
scratch=$SCRATCH
seqs=$WORK/vgp+477.fa.gz
out=$SCRATCH/vgp+477_wfmash-v0.13.0
mkdir -p $out

id=$SLURM_ARRAY_TASK_ID
target=$(head -$id $base/refs.txt | tail -1)
prefix=$target'#'

echo "mapping against $id = $target on $(hostname)"

$time $wfmash -t $CPUS \
      -m \
      -P $prefix \
      --one-to-one \
      -Y '#' \
      -n 1 \
      -k 19 \
      -p 70 \
      -s 5k \
      -c 20k \   
      $seqs \
      >$out/$target.map.paf \
      2>$out/$target.map.log \
      && touch $out/$target.map.ok

echo "aligning against $id = $target on $(hostname)"

$time $wfmash -t $CPUS \
      -i $out/$target.map.paf \
      $seqs \
      >$out/$target.aln.paf \
      2>$out/$target.aln.log \
      && touch $out/$target.aln.ok      
