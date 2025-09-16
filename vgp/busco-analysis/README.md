# BUSCO analysis

## Paths

The following paths are used in the analysis:

```shell
qbed=/lizardfs/guarracino/vgp/GCA_023851605.1.busco_vertebrata.bed
tbed=/lizardfs/guarracino/vgp/GCA_009914755.4.busco_vertebrata.bed
paf=/lizardfs/guarracino/vgp/alignment-vs-human_lz/GCA_023851605.1.aln.paf.gz
```

`impg` has to be in the `PATH`.

## Analysis

Full analysis:

```shell
python analyze_busco_alignments.py \
    --query-bed $qbed \
    --target-bed $tbed \
    --paf $paf \
    --threads 16 \
    --output .gene_details.tsv --summary-output summary_stats.tsv --coverage-output coverage_stats.tsv
```

Single gene analysis:

```shell
python analyze_busco_alignments.py \
    --query-bed $qbed \
    --target-bed $tbed \
    --paf $paf \
    --threads 16 \
    --output gene_details.11584at7742.tsv --summary-output summary_stats.11584at7742.tsv \
    --gene 11584at7742
```