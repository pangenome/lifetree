#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=thin
#SBATCH --time=48:00:00

# map the vgp290 to itself at 90% identity 5kb segments

nix_run=/home/egarrison/bin/nix-it
wfmash=/nix/store/gpzlqsy0b5sjywwy1nhvfmrivwbqss7i-wfmash-0.10.7/bin/wfmash
time=/usr/bin/time

base=/gpfs/nvme1/0/egarrison/lifetree
seqs=$base/asm/vgp290.fa.gz
out=$base/map/p90
mkdir -p $out

id=$SLURM_ARRAY_TASK_ID
target=$(head -$id $base/asm/vgp290.names.txt | tail -1)
prefix=$target'#'

echo "mapping against $id = $target"

$time -v $nix_run \
      $wfmash \
      -t 32 \
      -m \
      -Y# \
      -P $prefix \
      -n 1 \
      -p 90 \
      -s 5k \
      -c 20k \
      $seqs \
      >$out/$target.paf 2>$out/$target.log \
      && touch $out/$target.ok
