#!/bin/bash

mkdir genome

cd genome

trim_galore --paired *.fastq.gz

bwa index GCF_000182925.2_NC12_genomic.fna

bwa mem GCF_000182925.2_NC12_genomic.fna SRR5177532_2_val_2.fq.gz SRR5177532_1_val_1.fq.gz > SRR5177532.sam

samtools view -bS SRR5177532.sam | samtools sort -o SRR5177532_sorted.bam && samtools index SRR5177532_sorted.bam

bcftools mpileup -f GCF_000182925.2_NC12_genomic.fna SRR5177532_sorted.bam | bcftools call -mv -Oz -o SRR5177532_variants.vcf.gz
(Note: none of --samples-file, --ploidy or --ploidy-file given, assuming all sites are diploid)

bcftools filter -s LOWQUAL -e 'QUAL<20 || INFO/DP<10' -Oz -o SRR5177532_variants_filtered.vcf.gz SRR5177532_variants.vcf.gz

gunzip *.vcf.gz

igv
