awk '{print $1"\t0\t"$2}' Mammals/GCA_009914755.4.fna.gz.fai | bedtools makewindows -b - -w 100000 | \
parallel -j 30 --colsep '\t' \
 'impg query -p ./vgp-chm13.paf -r {1}:{2}-{3} \
 --transitive \
 --sequence-files genome/*.fna.gz \
 -o fasta > VGP/tmp/{2}-{3}.fa && samtools faidx VGP/tmp/{2}-{3}.fa
 FastGA -pafx VGP/tmp/{2}-{3}.fa VGP/tmp/{2}-{3}.fa> VGP/tmp/{2}-{3}.paf 
 #wfmash VGP/tmp/{2}-{3}.fa > VGP/tmp/{2}-{3}.paf 
 seqwish -s VGP/tmp/{2}-{3}.fa -p VGP/tmp/{2}-{3}.paf -g VGP/tmp/{2}-{3}.gfa
 odgi similarity -i VGP/tmp/{2}-{3}.gfa --delim "#" --delim-pos 1 --distances --all > VGP/out/{1}_{2}_{3}.similarity.tsv
 Rscript make_trees.R VGP/out/{1}_{2}_{3}.similarity.tsv'
