#!/bin/bash

conda activate syri
#mkdir salmon_q
cd salmon_q/
# wget http://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
# salmon index -p 4 -t ./Mus_musculus.GRCm38.cdna.all.fa.gz -i ./index
# rm -rf Mus_musculus.GRCm38.cdna.all.fa.gz
index="./index/"
cat ../SRR_Acc_List.txt | while read id ; do
echo "$id start salmon"

salmon quant -i $index -l A -1 ../clean_fastq/${id}*val*_1.fq* -2 ../clean_fastq/${id}*val*_2.fq* -p 6 --output ${id}_quant

echo "$id salmon done"


if [ -e ./${id}_quant/quant.sf ]; then
echo "$id salmon quant done"
sed 's/\.[0-9]*//' ./${id}_quant/quant.sf > ./${id}_quant/quant_clean.sf
else
rm -rf ./${id}_quant/
echo "$id salmon quant failed!"
fi
done

