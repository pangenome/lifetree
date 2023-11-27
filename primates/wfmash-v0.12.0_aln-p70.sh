#!/bin/bash
#SBATCH --job-name=alignment
#SBATCH --account=TG-MCB140147
#SBATCH --constraint="lustre"
#SBATCH --output=map_%j.out
#SBATCH --error=map_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --time=2-00:00:00
#SBATCH --partition=compute

# align t2t-primates to itself at 70% identity 5kb segments
# skipping self alignments

wfmash=~/wfmash/wfmash-01f812e5.sif
time=/usr/bin/time

module load singularitypro

# Determine number of CPUs to use
CPUS=$(($SLURM_NTASKS_PER_NODE))
base=/expanse/lustre/scratch/egarrison/temp_project/t2t-primates
seqs=$base/primates13.fa.gz
in=$base/map/p70_v1
out=$base/aln/p70_v1
mkdir -p $out

id=$SLURM_ARRAY_TASK_ID
target=$(head -$id $base/refs.txt | tail -1)
prefix=$target'#'

echo "mapping against $id = $target on $(hostname)"

$time singularity run --bind /scratch,/expanse \
      $wfmash -t $CPUS \
      -i $in/$target.paf \
      $seqs \
      >$out/$target.paf \
      2>$out/$target.log \
      && touch $out/$target.ok      
