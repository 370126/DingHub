#!/bin/bash

## open the well-prepared environment
conda active pan_genomics
mkdir pan_genomics
cd pan_genomics

## check the raw data
ls -l

## extract chrI in each genome and index it.
for file in *_genome.tidy.fa; do
idx=$(basename "$file" .nuclear_genome.tidy.fa)
pan_sn="${idx}#1#chrI"
seqkit faidx --quiet $file chrI | sed "s/>chrI/>${pan_sn}/" > ${pan_sn}.fa
    ## if data is not well-prepared (sorted by chromosomes), delete/move away.
if [ ! -s "${pan_sn}.fa" ]; then
rm ${pan_sn}.fa
mv $file ../
rm ${file}.fai
echo "${file} do not match requirement! Moved out :("
fi
done

## combine the *chrI.fa and index it.
mkdir raw_chrI
mkdir out_chrI
mv *chrI.fa ./raw_chrI
cat ./raw_chrI/*.fa > ./out_chrI/total_chrI.fa
samtools faidx ./out_chrI/total_chrI.fa
num=$( grep ">.*" ./out_chrI/total_chrI.fa | wc -l )
echo "total genomes number: ${num}"

## run pggb
cd out_chrI
mkdir out_pggb
pggb -i total_chrI.fa -o out_pggb -n ${num} -t 5 -p 90

## dotplot
wget https://raw.githubusercontent.com/waveygang/wfmash/refs/heads/main/scripts/paf2dotplot
./paf2dotplot png large ./out_pggb/*.wfmash.paf

## pangenome graph statistics
odgi stats -i ./out_pggb/*.smooth.final.og -S > ./out_pggb/out_pggb_chrI.fa.smooth.final.og.stats

## vg dive into pangenome structure
mkdir out_vg
    # 450:500
vg chunk -x ./out_pggb/out_pggb_chrI.smooth.final.vg -p AIC.asm01.HP0#1#chrI#0:450-500 -c 1 | vg view -pd - | dot -Tpng > ./out_vg/out_pggb_chrI.smooth.final.AIC_chrI_450-500_real.png
    # 34500:35500
vg chunk -x ./out_pggb/out_pggb_chrI.smooth.final.vg -p AIC.asm01.HP0#1#chrI#0:35450-35500 -c 1 | vg view -pd - | dot -Tpng > ./out_vg/out_pggb_chrI.smooth.final.AIC_chrI_35450-35500_real.png

## run panacus histgrowth to calculate coverage and pangenome growth for nodes (default) with coverage/quorum thresholds 1/0, 2/0, 1/1, 1/0.5, and 1/0.1
RUST_LOG=info panacus histgrowth -t4 -l 1,2,1,1,1 -q 0,0,1,0.5,0.1 -S -a \
    ./out_pggb/total_chrI.fa.bf3285f.11fba48.78baeea.smooth.final.gfa \
    > ./out_pggb/total_chrI.fa.bf3285f.11fba48.78baeea.smooth.final.histgrowth.node.tsv

## visualize the growth statistics
panacus-visualize -e ./out_pggb/total_chrI.fa.bf3285f.11fba48.78baeea.smooth.final.histgrowth.node.tsv --split_subfigures -f png --split_prefix chrI_panacus_ -s 20 12

cd ..
