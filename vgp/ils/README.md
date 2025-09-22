# ILS

```shell
mkdir apes-chm13
awk '{print $1"\t0\t"$2}' /lizardfs/guarracino/vgp/apes/chm13v2.pansn.fa.gz.fai  | bedtools makewindows -b - -w 25000 | \
parallel -j 20 --colsep '\t' \
    'impg similarity -p /scratch/chm13v2+mPanPan1P+mPanTro3h2.paf.gz \
     -r {1}:{2}-{3} \
     --sequence-files /lizardfs/guarracino/vgp/apes/*.pansn.fa.gz \
     --force-large-region --delim "#" --delim-pos 1 --distances --all \
     > apes-chm13/apes.{1}_{2}_{3}.similarity.tsv && \
     Rscript scripts/make_tree.R apes-chm13/apes.{1}_{2}_{3}.similarity.tsv'

Rscripts scripts/compare_trees.R
```
