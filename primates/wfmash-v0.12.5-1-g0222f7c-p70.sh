#!/bin/bash
#SBATCH --job-name=wfmash
#SBATCH --account=TG-MCB140147
#SBATCH --constraint="lustre"
#SBATCH --output=wfmash_%j.out
#SBATCH --error=wfmash_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=249208M
#SBATCH --time=2-00:00:00
#SBATCH --partition=compute

# map t2t-primates to itself at 70% identity 5kb segments
# skipping self alignments

wfmash=~/wfmash/wfmash-v0.12.5-1-g0222f7c.sif
time=/usr/bin/time

module load singularitypro

# Determine number of CPUs to use
CPUS=$(($SLURM_NTASKS_PER_NODE))
base=/expanse/lustre/scratch/egarrison/temp_project/t2t-primates
seqs=$base/primates16.20231205.fa.gz
out=$base/primates16.20231205_wfmash-v0.12.5
mkdir -p $out

id=$SLURM_ARRAY_TASK_ID
target=$(head -$id $base/refs.txt | tail -1)
prefix=$target'#'

echo "mapping against $id = $target on $(hostname)"

$time singularity run --bind /scratch,/expanse \
      $wfmash -t $CPUS \
      -m \
      -P $prefix \
      --one-to-one \
      -Y '#' \
      -n 1 \
      -p 70 \
      -s 5k \
      -c 20k \
      $seqs \
      >$out/$target.map.paf \
      2>$out/$target.map.log \
      && touch $out/$target.map.ok

echo "aligning against $id = $target on $(hostname)"

$time singularity run --bind /scratch,/expanse \
      $wfmash -t $CPUS \
      -i $out/$target.map.paf \
      $seqs \
      >$out/$target.aln.paf \
      2>$out/$target.aln.log \
      && touch $out/$target.aln.ok      
