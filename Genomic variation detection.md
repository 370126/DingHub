## 第一次作业 基因组变异检测

**221505023 张牧原**

---

### 1. 数据准备

#### 1.1 选择物种及参考基因组

- 选择物种：粗糙脉胞菌  *Neurospora Crassa* OR74A

- 参考基因组来源：[NCBI-genome-NC12](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000182925.2/)
- GCF_000182925.2
- 基本信息：
    - 7条染色体
    - length: 41 Mb
    - GC% : 48.5
- 参考序列文件：GCF_000182925.2_NC12_genomic.fna (41.6 Mb)

#### 1.2 重测序数据集获取

- 数据来源：[NCBI-SRA](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR5177532&display=download)
- 搜索方法: ("Neurospora crassa OR74A"[Organism] OR Neurospora crassa OR74A[All Fields]) AND ("strategy wgs"[Properties] AND "library layout paired"[Properties] AND "filetype fastq"[Properties])
- 基本信息：
    - Name: FGSC2489
    - Instrument: Illumina HiSeq 2500
    - Read Count: 18924722
    - Base Count: 3721527023
    - Strategy: WGS
    - Source: GENOMIC
    - Selection: RANDOM
    - Layout: PAIRED
- 数据文件：
    - SRR5177532_1.fastq.gz (1.6 G)
    - SRR5177532_2.fastq.gz (1.6 G)

### 2. 数据分析


1. 按下面的脚本，对数据进行处理：

```{shell}
mkdir genome
cd ./genome

trim_galore --paired *.fastq.gz

bwa index GCF_000182925.2_NC12_genomic.fna

bwa mem GCF_000182925.2_NC12_genomic.fna SRR5177532_2_val_2.fq.gz SRR5177532_1_val_1.fq.gz > SRR5177532.sam

samtools view -bS SRR5177532.sam | samtools sort -o SRR5177532_sorted.bam && samtools index SRR5177532_sorted.bam

bcftools mpileup -f GCF_000182925.2_NC12_genomic.fna SRR5177532_sorted.bam | bcftools call -mv -Oz -o SRR5177532_variants.vcf.gz
(Note: none of --samples-file, --ploidy or --ploidy-file given, assuming all sites are diploid)

bcftools filter -s LOWQUAL -e 'QUAL<20 || INFO/DP<10' -Oz -o SRR5177532_variants_filtered.vcf.gz SRR5177532_variants.vcf.gz

gunzip *.vcf.gz

igv
```

2. 包括最初的参考序列以及重测序文件，最终该目录下的文件有：
```
GCF_000182925.2_NC12_genomic.fna
GCF_000182925.2_NC12_genomic.fna.amb
GCF_000182925.2_NC12_genomic.fna.ann
GCF_000182925.2_NC12_genomic.fna.bwt
GCF_000182925.2_NC12_genomic.fna.fai
GCF_000182925.2_NC12_genomic.fna.pac
GCF_000182925.2_NC12_genomic.fna.sa
SRR5177532_1.fastq.gz
SRR5177532_1.fastq.gz_trimming_report.txt
SRR5177532_1_val_1.fq.gz
SRR5177532_2.fastq.gz
SRR5177532_2.fastq.gz_trimming_report.txt
SRR5177532_2_val_2.fq.gz
SRR5177532.bam
SRR5177532.sam
SRR5177532_sorted.bam
SRR5177532_sorted.bam.bai
SRR5177532_variants_filtered.vcf.gz
SRR5177532_variants.vcf.gz
```

### 3. 使用igv可视化结果

1. 加载参考基因组
> 左上角'Genomes' -> 'Load Genome from File...' -> 'GCF_000182925.2_NC12_genomic.fna'

2. 载入.bam文件显示覆盖度

> 左上角'File' -> 'Load from File...' -> 'SRR5177532_sorted.bam'

点击某一条染色体，然后放大，则会看到可视化的效果。

有两个道：上面的一道名字为'SRR5177532_sorted.bam Coverage', 即每个位点的覆盖度；下面一道名字为'SRR5177532_sorted.bam', 所有的read的覆盖情况。

可以点击鼠标右键对道的样貌进行设置。对于道的大小有'Collapsed', 'Expanded', 'Squished'三种模式，其中'Expanded'展示信息最多，也更宽; 'Squished'最紧凑; 'Collapsed'适中。

3. 载入.vcf文件显示基因组的SNP/INDEL位点

```
gunzip SRR5177532_variants_filtered.vcf.gz
```
> 左上角'File' -> 'Load from File...' -> 'SRR5177532_variants_filtered.vcf'

添加了两条新的道，一条为'SRR5177532_variants_filtered.vcf', 另一条为'SRR5177532_sorted.bam'
注意此时如果视野较广(即所查看范围较大)，能看见视野中存在许多细蓝条(但此时因视野过大可能无法查看2.中的reads覆盖)，这些即是SNP/INDEL位点，点击相应的位置可以查看相关突变信息。

若想要测序reads覆盖度和SNP/INDEL位点同时呈现在画面中，应找到对应的SNP/INDEL位点，然后逐渐放大，确保该变异位点始终在视野当中直至显示reads覆盖信息。

4. 实例展示

- SNP
<img src="https://github.com/370126/DingHub/blob/main/FIGS/genomics/SNP.png" width="100000px"/>

- INDEL
<img src="https://github.com/370126/DingHub/blob/main/FIGS/genomics/INDEL.png">
