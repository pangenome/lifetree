# ILS

## Full(ish) IMPG way

Limited by SPOA memory usage:

```shell
mkdir apes-chm13
awk '{print $1"\t0\t"$2}' /lizardfs/guarracino/vgp/apes/chm13v2.pansn.fa.gz.fai | bedtools makewindows -b - -w 25000 | \
parallel -j 48 --colsep '\t' \
    'impg similarity -p /lizardfs/guarracino/vgp/apes/chm13v2+mPanPan1P+mPanTro3h2.paf.gz \
     -r {1}:{2}-{3} \
     --sequence-files /lizardfs/guarracino/vgp/apes/*.pansn.fa.gz \
     --force-large-region --delim "#" --delim-pos 1 --distances --all \
     > apes-chm13/apes.{1}_{2}_{3}.similarity.tsv && \
     Rscript /lizardfs/guarracino/vgp/apes/make_tree.R apes-chm13/apes.{1}_{2}_{3}.similarity.tsv'

Rscript /lizardfs/guarracino/vgp/apes/compare_trees.R
```

## IMPG + <SEQUENCE-ALIGNER> + SEQWISH + ODGI way

It allows larger windows:

```shell
mkdir apes-chm13
awk '{print $1"\t0\t"$2}' /lizardfs/guarracino/vgp/apes/chm13v2.pansn.fa.gz.fai | bedtools makewindows -b - -w 500000 | \
parallel -j 48 --colsep '\t' \
    'impg query -p /lizardfs/guarracino/vgp/apes/chm13v2+mPanPan1P+mPanTro3h2.paf.gz \
     -r {1}:{2}-{3} \
     --transitive \
     --sequence-files /lizardfs/guarracino/vgp/apes/*.pansn.fa.gz \
     -o fasta | seqkit rmdup -s - \
     > apes-chm13/apes.{1}_{2}_{3}.fa && samtools faidx apes-chm13/apes.{1}_{2}_{3}.fa && \
     wfmash apes-chm13/apes.{1}_{2}_{3}.fa > apes-chm13/apes.{1}_{2}_{3}.paf && \
     seqwish -s apes-chm13/apes.{1}_{2}_{3}.fa -p apes-chm13/apes.{1}_{2}_{3}.paf -g apes-chm13/apes.{1}_{2}_{3}.gfa --temp-dir /scratch && \
     odgi similarity -i apes-chm13/apes.{1}_{2}_{3}.gfa --delim "#" --delim-pos 1 --distances --all > apes-chm13/apes.{1}_{2}_{3}.similarity.tsv && \
     Rscript /lizardfs/guarracino/vgp/apes/make_tree.R apes-chm13/apes.{1}_{2}_{3}.similarity.tsv'
     
# "seqkit rmdup -s -" to patch a buggy behavior of impg query that outputs duplicated sequences in FASTA formats with transitive queries

Rscript /lizardfs/guarracino/vgp/apes/compare_trees.R
```
