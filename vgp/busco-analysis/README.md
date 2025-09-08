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
    --output gene_details.tsv --summary-output summary_stats.tsv

cat summary_stats.tsv | column -t
total_analyzed  missing  incomplete  wrong_location  good
3378            2        2669        592             115

cat summary_stats.tsv | column -t
gene_id        target_name                   target_start  target_end  target_length  gene_present_in_query  gene_present_in_alignments  query_n_total_alignments  query_n_busco_alignments  query_n_non_busco_alignments  query_busco_completeness
1000205at7742  GCA_009914755.4#0#CP068262.2  74097059      74168116    71057          Yes                    Yes                         185                       4                         181                           1.0000
1001576at7742  GCA_009914755.4#0#CP068273.2  141222171     141231147   8976           Yes                    Yes                         64                        3                         61                            1.0000
1003583at7742  GCA_009914755.4#0#CP068266.2  113172200     113210057   37857          Yes                    Yes                         165                       5                         160                           0.9004
1005433at7742  GCA_009914755.4#0#CP068261.2  47128637      47130701    2064           Yes                    Yes                         5                         1                         4                             1.0000
1006116at7742  GCA_009914755.4#0#CP068274.2  17793208      17825769    32561          Yes                    Yes                         522                       13                        509                           0.1244
1006359at7742  GCA_009914755.4#0#CP068267.2  108149951     108154935   4984           Yes                    Yes                         5                         3                         2                             0.1671
1007560at7742  GCA_009914755.4#0#CP068277.2  204339358     204348699   9341           Yes                    Yes                         20                        4                         16                            0.5543
1007920at7742  GCA_009914755.4#0#CP068277.2  223416549     223519008   102459         Yes                    Yes                         999                       32                        967                           0.0917
1008884at7742  GCA_009914755.4#0#CP068273.2  141173195     141187722   14527          Yes                    Yes                         46                        11                        35                            0.8586
1009797at7742  GCA_009914755.4#0#CP068276.2  85288218      85392937    104719         Yes                    Yes                         2117                      13                        2104                          0.5487
1010709at7742  GCA_009914755.4#0#CP068275.2  45370497      45563886    193389         Yes                    Yes                         2569                      56                        2513                          0.5375
1014237at7742  GCA_009914755.4#0#CP068275.2  118383874     118441093   57219          Yes                    Yes                         1396                      17                        1379                          0.4652
...
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