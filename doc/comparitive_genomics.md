## 第二次生信作业<!--no toc--> 

***221505023 张牧原***

---


- Exercise 02: Compare and plot three genomes


---

### 目录<!--no toc--> 


- [**ONE. Compare and plot three genomes**](#one-compare-and-plot-three-genomes)
  - [I. 数据准备](#i-数据准备)
  - [II.基因组比对](#ii基因组比对)
  - [III.结构变异检测统计](#iii结构变异检测统计)
  - [IV.plotsr可视化结果](#ivplotsr可视化结果)


---

### **ONE. Compare and plot three genomes**

### I. 数据准备

#### 1.1. 用于比较的三个酵母基因组和来源

| species | strain | data source                                                                     |
| ------- | ------ | ------------------------------------------------------------------------------- |
| *S.c.*  | SK1    | http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_Genome/SK1.genome.fa.gz |
| *S.c.*  | Y12    | http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_Genome/Y12.genome.fa.gz |
| *S.p.*  | N44    | http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_Genome/N44.genome.fa.gz |

#### 1.2. 用于plotsr的配置文件

根据plotsr的[使用手册](https://github.com/schneebergerlab/plotsr)，可以添加输入文件用于**1.说明基因组信息；2.设置染色体排列顺序；3.添加额外的轨道**。因此，我使用三个文件来完成相关设置：

1. genomes.txt
```
/home/ubuntu/Desktop/advanced_bioinfo/compare_genomics/SK1.genome.fa	SK1	lw:1.5
/home/ubuntu/Desktop/advanced_bioinfo/compare_genomics/Y12.genome.fa	Y12	lw:1.5
/home/ubuntu/Desktop/advanced_bioinfo/compare_genomics/N44.genome.fa	N44	lw:1.5
```
2. chrord.txt
```
chrI
chrII
chrIII
chrIV
chrV
chrVI
chrVII
chrVIII
chrIX
chrX
chrXI
chrXII
chrXIII
chrXIV
chrXV
chrXVI
```
> 注意：该文件的总行数与染色体总数应保持一致，否则plotsr会将空行视作一条染色体而报错！

3. track.txt

在这里根据练习的附加要求，选择添加SK1菌株的基因密度作为新的轨道。

输入数据要求为.gff格式，数据来源：http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_GFF/SK1.all_feature.gff.gz

```
SK1.all_feature.gff	Genes on SK1	ft:gff;bw:100000;lc:blue;lw:4;bc:lightblue;ba:0.5
```
### II.基因组比对

当需要三个以上基因之间进行比较基因组学分析时，不必所有基因组之间两两比对，只需“**首尾相接**”即可。

e.g. 有三个基因组 **A.fa, B.fa, C.fa**

那么如果使用plotsr整合三个基因组的比较结果，只需保证:
- Alignment 1: -r A.fa -q B.fa
- Alignment 2: -r B.fa -q C.fa

具体操作步骤见下。

1. 使用mummer对基因组进行比对

对于同种菌株之间(SK1/Y12)，比对参数可设置得稍严格些；对于不同菌株之间(Y12/N44)，比对参数可稍放松。
相关的参数含义是：

- *-c*    Sets the minimum length of a cluster of matches (default 65)
- *-b*    Set the distance an alignment extension will attempt to extend poor scoring regions before giving up (default 200)
- *-l*    Set the minimum length of a single match (default 20)
- *-p*    Set the prefix of the output files (default "out")


```{shell}
nucmer --maxmatch -c 200 -b 500 -l 100 -p SK1_vs_Y12_nucmer SK1.genome.fa Y12.genome.fa
nucmer --maxmatch -c 50 -l 20 -p Y12_vs_N44_nucmer Y12.genome.fa N44.genome.fa
```

其输出文件为 **.delta**文件。

2. 过滤.delta文件，去除不符合要求的比对片段

相关参数的含义是：

- *-m*    Many-to-many alignment allowing for rearrangements (union of -r and -q alignments)
- *-i*    Set the minimum alignment identity [0, 100], default 0
- *-l*    Set the minimum alignment length, default 0
- (*-r*)  Maps each position of each reference to its best hit in the query, allowing for query overlaps
- (*-q*)  Maps each position of each query to its best hit in the reference, allowing for reference overlaps


```{shell}
delta-filter -m -i 90 -l 100 SK1_vs_Y12_nucmer.delta > SK1_vs_Y12_nucmer.filter
delta-filter -m -i 90 -l 100 Y12_vs_N44_nucmer.delta > Y12_vs_N44_nucmer.filter
```

3. 将.delta文件转化为更可读的.coords文件

```{shell}
show-coords -THrd SK1_vs_Y12_nucmer.filter > SK1_vs_Y12_nucmer.filter.coords
show-coords -THrd Y12_vs_N44_nucmer.filter > Y12_vs_N44_nucmer.filter.coords
```

### III.结构变异检测统计

相关参数的含义是：

- *-c*  File containing alignment coordinates (default: None)
- *-d*  .delta file from mummer. Required for short variation (SNPs/indels) identification when CIGAR string is not available (default: None)
- *-r*  Genome A (which is considered as reference for the alignments). Required for local variation (large indels, CNVs) identification. (default: None)
- *-q*  Genome B (which is considered as query for the alignments). Required for local variation (large indels, CNVs) identification. (default: None)
- *--nc*    number of cores to use in parallel (max is number of chromosomes) (default: 1)
- *--dir*   path to working directory (if not current directory). All files must be in this directory. (default: None)
- *--prefix*    Prefix to add before the output file Names (default: )
- *--lf*    Name of log file (default: syri.log)


```
syri -c SK1_vs_Y12_nucmer.filter.coords -d SK1_vs_Y12_nucmer.filter -r SK1.genome.fa -q Y12.genome.fa --nc 5 --dir SK1_vs_Y12/ --prefix SK1_vs_Y12_nucmer.filter. --lf SK1_vs_Y12_nucmer.filter.syri.log
syri -c Y12_vs_N44_nucmer.filter.coords -d Y12_vs_N44_nucmer.filter -r Y12.genome.fa -q N44.genome.fa --nc 5 --dir Y12_vs_N44/ --prefix Y12_vs_N44_nucmer.filter. --lf Y12_vs_N44_nucmer.filter.syri.log
```


### IV.plotsr可视化结果

相关参数的含义是：

- *--sr*    Structural annotation mappings (syri.out) identified by SyRI (default: None)
- *--genomes*   File containing path to genomes (default: None)
- *--chrord*    File containing reference (first genome) chromosome IDs in the order in which they are to be plotted. File requires one chromosome ID per line. Not compatible with --chr (default: None)
- *--tracks*     File listing paths and details for all tracks to be plotted (default: None)
- *-o*  Output file name. Acceptable format: pdf, png, svg (default: plotsr.pdf)

```{shell}
## make sure you've prepared the configuration files.
plotsr --sr ./SK1_vs_Y12/*.out --sr ./Y12_vs_N44/*.out --genomes genomes.txt  --chrord chrord.txt --tracks track.txt -d 600 -W 10 -H 14 -f 8 -o SK1_Y12_N44.plotsr.pdf
```
可视化的结果如下。

![FIG.1](/home/ubuntu/Desktop/advanced_bioinfo/compare_genomics/SK1_Y12_N44.plotsr.png "Fig.1 Compare 3 genomes of S.p. strains")

---

下面是实现整个分析流程的 **shell script**:

```{shell}
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

```
