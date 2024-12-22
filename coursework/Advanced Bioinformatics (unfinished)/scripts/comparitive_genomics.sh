#!/bin/bash

conda activate syri
## make sure you've done the setting environment for SyRI


mkdir compare_genomics
cd ./compare_genomics

## files download
# wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_Genome/SK1.genome.fa.gz
# wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_Genome/Y12.genome.fa.gz
# wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_Genome/N44.genome.fa.gz
# wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_GFF/SK1.all_feature.gff.gz

## configuration files for plotsr (https://github.com/schneebergerlab/plotsr)
# touch genomes.txt
# touch chrord.txt
# touch track.txt

ls -l

## files required in current working directory:
## SK1.genome.fa.gz
## Y12.genome.fa.gz
## N44.genome.fa.gz
## SK1.all_feature.gff.gz
## genomes.txt
## chrord.txt
## track.txt

gunzip *.gz

## Genomes Alignment
nucmer --maxmatch -c 200 -b 500 -l 100 -p SK1_vs_Y12_nucmer SK1.genome.fa Y12.genome.fa
nucmer --maxmatch -c 50 -l 20 -p Y12_vs_N44_nucmer Y12.genome.fa N44.genome.fa

delta-filter -m -i 90 -l 100 SK1_vs_Y12_nucmer.delta > SK1_vs_Y12_nucmer.filter
delta-filter -m -i 90 -l 100 Y12_vs_N44_nucmer.delta > Y12_vs_N44_nucmer.filter

show-coords -THrd SK1_vs_Y12_nucmer.filter > SK1_vs_Y12_nucmer.filter.coords
show-coords -THrd Y12_vs_N44_nucmer.filter > Y12_vs_N44_nucmer.filter.coords


mkdir SK1_vs_Y12
# mkdir SK1_vs_N44
mkdir Y12_vs_N44


## SV detection
syri -c SK1_vs_Y12_nucmer.filter.coords -d SK1_vs_Y12_nucmer.filter -r SK1.genome.fa -q Y12.genome.fa --nc 5 --dir SK1_vs_Y12/ --prefix SK1_vs_Y12_nucmer.filter. --lf SK1_vs_Y12_nucmer.filter.syri.log
syri -c Y12_vs_N44_nucmer.filter.coords -d Y12_vs_N44_nucmer.filter -r Y12.genome.fa -q N44.genome.fa --nc 5 --dir Y12_vs_N44/ --prefix Y12_vs_N44_nucmer.filter. --lf Y12_vs_N44_nucmer.filter.syri.log


## Plot syri
plotsr --sr ./SK1_vs_Y12/*.out --sr ./Y12_vs_N44/*.out --genomes genomes.txt  --chrord chrord.txt --tracks track.txt -d 600 -W 10 -H 14 -f 8 -o SK1_Y12_N44.plotsr.pdf
