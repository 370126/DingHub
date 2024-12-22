#!/bin/bash
mkdir RNA-seq
cd RNA-seq

mkdir raw_fastq
mkdir SRR

list="./SRR_Acc_List.txt"

for SRR  in $(cat "$list"); do
echo "$SRR prefetch begin"
prefetch -p $SRR -O SRR
echo "$SRR prefetch done"
done


list="./SRR_Acc_List.txt"
for SRR  in $(cat "$list"); do
echo "$SRR fastq-dump begin"
fastq-dump --gzip --skip-technical --readids --dumpbase --split-3 --clip ./SRR/$SRR -O ./raw_fastq/
echo "$SRR fastq-dump maybe done"
done

## if disk is limited
# rm -rf ./SRR


mkdir fastQC

ls raw_fastq | xargs fastqc -o fastQC -t 12

multiqc ./fastQC/ -o ./fastQC/


conda activate syri
mkdir clean_fastq

ls ./raw_fastq/*_1* > _1_id.txt
ls ./raw_fastq/*_2* > _2_id.txt

paste _1_id.txt _2_id.txt > IDs.txt
rm *_id.txt

cat IDs.txt | while read trial 
do
temp=($trial)
fq1=${temp[0]}
fq2=${temp[1]}

# echo "fq1=$fq1 fq2=$fq2"

trim_galore -j 4 -q 25 --phred33 --length 30 --stringency 3 \
            --paired --gzip -o ./clean_fastq/ $fq1 $fq2

done


mkdir fastQC_trimmed
ls ./clean_fastq/*val*fq* | xargs fastqc -o fastQC_trimmed -t 12
multiqc ./fastQC_trimmed/ -o ./fastQC_trimmed

## if disk is limited
# rm -rf ./raw_fastq
