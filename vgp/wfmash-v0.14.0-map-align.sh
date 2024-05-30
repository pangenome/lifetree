#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2-00:00:00

# break on error
set -euo pipefail

# our wfmash binary
wfmash=$WORK/micromamba/bin/wfmash

seqs=$1
jobs=$2
id=$3
out=$4

cpus=$SLURM_CPUS_PER_TASK
time=/usr/bin/time
ls -l $seqs
mkdir -p $out

spec=$(head -$id $jobs | tail -1)
target=$(echo $spec | awk '{print $1}')
query=$(echo $spec | awk '{print $2}')

echo "wfmash mapping $id = $target against $query on $(hostname) with $cpus CPUs at $(date +%s) / $(date) with $wfmash "
$wfmash --version

map_out=$out/$id.$target.map
$time -vv $wfmash -t $cpus \
      -m \
      -T $target \
      -Q $query \
      -Y '#' \
      -n 1 \
      -k 19 \
      -p 70 \
      -s 5k \
      -c 20k \
      $seqs \
      >$map_out.paf \
      2>$map_out.log \
    && touch $map_out.ok

echo "wfmash aligning $id = $target against $query on $(hostname) with $cpus CPUs at $(date +%s) / $(date)"

aln_out=$out/$id.$target.aln
$time -vv $wfmash -t $cpus \
      -i $out/$target.id_$id.map.paf \
      $seqs \
      >$aln_out.paf \
      2>$aln_out.log \
    && touch $aln_out.ok

echo "zipping $map_out.paf and $aln_out.paf"
pigz -f $map_out.paf $aln_out.paf

touch $id.$target.done
