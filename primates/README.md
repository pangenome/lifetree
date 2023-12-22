# wfmash 16 primate assemblies

## Preparation (release 2023/12/05)

Download assemblies:

```shell
# Primates from https://genomeark.s3.amazonaws.com/index.html?prefix=species/
wget -c https://s3.amazonaws.com/genomeark/species/Gorilla_gorilla/mGorGor1/assembly_curated/mGorGor1.dip.cur.20231122.fasta.gz
wget -c https://s3.amazonaws.com/genomeark/species/Pan_paniscus/mPanPan1/assembly_curated/mPanPan1.dip.cur.20231122.fasta.gz
wget -c https://s3.amazonaws.com/genomeark/species/Pan_troglodytes/mPanTro3/assembly_curated/mPanTro3.dip.cur.20231122.fasta.gz
wget -c https://s3.amazonaws.com/genomeark/species/Pongo_abelii/mPonAbe1/assembly_curated/mPonAbe1.dip.cur.20231205.fasta.gz
wget -c https://s3.amazonaws.com/genomeark/species/Pongo_pygmaeus/mPonPyg2/assembly_curated/mPonPyg2.dip.cur.20231122.fasta.gz
wget -c https://genomeark.s3.amazonaws.com/species/Symphalangus_syndactylus/mSymSyn1/assembly_curated/mSymSyn1.dip.cur.20231205.fasta.gz

# HG002v1.0.1
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.0.1.fasta.gz
```

Apply [PanSN-spec](https://github.com/pangenome/PanSN-spec):

```shell
for species in mGorGor1 mPanPan1; do
    zcat $species.dip.cur.20231122.fasta.gz | \
        sed "/^>/ s/chr\(.*\)_\(mat\|pat\)_*/$species##\2##\0/; s/#mat#/M/; s/#pat#/P/" | \
        bgzip -@ 48 -l 9 > $species.fa.gz && samtools faidx $species.fa.gz
done
for species in mPanTro3 mPonPyg2; do
    zcat $species.dip.cur.20231122.fasta.gz  | \
        sed "/^>/ s/chr\(.*\)_\(hap1\|hap2\)_*/$species##\2##\0/; s/#hap1#/1/; s/#hap2#/2/" | \
        bgzip -@ 48 -l 9 > $species.fa.gz && samtools faidx $species.fa.gz
done
for species in mPonAbe1 mSymSyn1; do
    zcat $species.dip.cur.20231205.fasta.gz  | \
        sed "/^>/ s/chr\(.*\)_\(hap1\|hap2\)_*/$species##\2##\0/; s/#hap1#/1/; s/#hap2#/2/" | \
        bgzip -@ 48 -l 9 > $species.fa.gz && samtools faidx $species.fa.gz
done
rm *fasta.gz

zcat hg002v1.0.1.fasta | sed -e 's/^>chr\(.*\)_MATERNAL/>hg002#M#chr\1_MATERNAL/' \
    -e 's/^>chr\(.*\)_PATERNAL/>hg002#P#chr\1_PATERNAL/' \
    -e 's/^>chrEBV/>hg002#P#chrEBV/' \
    -e 's/^>chrM/>hg002#M#chrM/' | bgzip -@ 48 -l 9 > hg002v101.fa.gz && samtools faidx hg002v101.fa.gz
rm hg002v1.0.1.fasta.gz
```

Put all haplotypes together (assumes you already have PanSN-ed CHM13v2 and PanSN-ed GRCh38):

```shell
zcat chm13v2.0.fa.gz grch38.fa.gz hg002v101.fa.gz mGorGor1.fa.gz mPanPan1.fa.gz mPanTro3.fa.gz mPonAbe1.fa.gz mPonPyg2.fa.gz mSymSyn1.fa.gz | bgzip -@ 48 -l 9 > primates16.20231205.fa.gz
samtools faidx primates16.20231205.fa.gz
```

## Mapping

The `wfmash-v0.12.0_map-p70.sh` script is used to generate mappings between all genomes and each target genome.
By running this for all genomes, we can evoke an all-to-all homology mapping and base-level alignment in parallel.

Mapping is accomplished by the `mashmap3` algorithm as parameterized in `wfmash -m`.
It's achieved by splitting each query genome into overlapping 5kb segments and mapping each segment to the chosen target.
Mappings are filtered to keep only those with >70% identity and a cumulative length >20kb. 

The script runs `wfmash` with these key parameters:

- `-m` - generate mappings (instead of alignments)  
- `-p 70` - minimum 70% identity
- `-s 5k` - segment size 5kb
- `-c 20k` - chains with cumulative 20kb mapping

This results in PAF format mapping files for each target sequence in the genome.

## Alignment

The `wfmash-v0.12.0_aln-p70.sh` script can then be used to convert the mappings to alignments.

It takes the PAF mapping files generated in the previous step as input. 

It runs `wfmash` again for each target, this time in alignment mode, using the mappings to guide the alignment process.

This results in aligned PAF files for each target sequence which include `--eqx` compatible cigars usable in rustybam, seqwish, and other tools.

## Visualization

The `plot.R` script is then used to visualize the mappings. In this example we use the homology mappings only, for simplicity.

It loads the mapping PAF file and bins the mappings into windows along each target chromosome using `bedtools intersect`.

Then for each target chromosome, it generates a heatmap showing which query sequences align to that target, filtering for queries above a minimum length.

The heatmaps are colored by query sequence identity and scaled based on the number of unique query sequences. 

This results in PDF and PNG heatmap images for each target chromosome showing the mapping coverage and origins.

## Usage

1. Generate mappings with `wfmash-v0.12.0_map-p70.sh` or alignments with `wfmash-v0.12.0_aln-p70.sh`.

2. Bin and process mappings into windows with `overlaps.sh`

3. Visualize each chromosome with `plot.R` 

For example:

```
# Generate mappings with an array job
sbatch -a 1-13 wfmash-v0.12.0_map-p70.sh 

# Bin mappings 
./overlaps.sh mappings.paf genome.fa.gz prefix 100000

# Visualize
./plot.R binned_mappings.bed.tsv 1000000
```

This will generate heatmap images for each chromosome showing the mapping results.

## Data

* Mappings only (data shown here) https://garrisonlab.s3.amazonaws.com/t2t-primates/primates13-v0.1-wfmash-01f812e5-map.tar.gz 
* Alignments plus mappings (full cigars for all homologies shown here) https://garrisonlab.s3.amazonaws.com/t2t-primates/primates13-v0.1-wfmash-01f812e5.tar.gz
* Assemblies with names formatted as shown here: https://garrisonlab.s3.amazonaws.com/t2t-primates/primates13.fa.gz 
