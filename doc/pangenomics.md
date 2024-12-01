## 第二次生信作业<!--no toc--> 

***221505023 张牧原***

---

所选练习：

- Exercise 06: Create pangenome graph for more chromosomes or more strains

---

### 目录<!--no toc--> 


- [**TWO. Create pangenome graph for more chromosomes or more strains**](#two-create-pangenome-graph-for-more-chromosomes-or-more-strains)
  - [0.环境配置](#0环境配置)
  - [I.数据/文件准备](#i数据文件准备)
  - [II.提取所有基因组特定染色体](#ii提取所有基因组特定染色体)
  - [III. Run *pggb* and related visualization](#iii-run-pggb-and-related-visualization)
  - [IV. Calculate pangenome openness with *panacus*](#iv-calculate-pangenome-openness-with-panacus)


---


### **TWO. Create pangenome graph for more chromosomes or more strains**

### 0.环境配置

由于我在配pggb环境时遇到了不少困难，有必要在此说明一下配置环境时的注意事项。

```{shell}
## 更新conda及所有的包
conda update conda
conda update --all

## 创建一个干净的环境，避免冲突
conda create -n pangenome
conda activate pangenome

## (optional) 更新gcc/g++
conda install -c conda-forge gcc_linux-64 gxx_linux-64

## 安装seqkit
conda install -c bioconda seqkit

## 安装pggb
conda install -c bioconda -c conda-forge pggb

## 如果出现'HTTP error'，考虑换源，参考 https://mirrors4.tuna.tsinghua.edu.cn/help/anaconda/

## 安装graphvz
sudo apt-get install graphviz  

## 安装panacus
conda install -c bioconda -c conda-forge panacus

## (optional) 安装pggb后wfmash好像并没有配置好，总是报错，可以单独安装一下
conda insall -c bioconda wfmash
```
上面这些命令并不能解决所有可能遇到的问题，比如清华的镜像会莫名其妙地拒绝路由器的网络访问（换成校园网/手机热点就没问题，但手机热点奇慢无比），有时候新建一个环境或者换一个网络就好了、、、

### I.数据/文件准备

1. 数据来源

从[ScRAP数据库](https://www.evomicslab.org/db/ScRAPdb/)中获取属于物种*Saccharomyces cerevisiae*的基因组序列文件。我随机选择了29个数据进行下载，如下所示：


```
273614N.asm01.HP0.nuclear_genome.tidy.fa
273614N.asm02.HP0.nuclear_genome.tidy.fa
ABL.asm01.HP0.nuclear_genome.tidy.fa
AFH.asm01.HP0.nuclear_genome.tidy.fa
AGK.asm01.HP0.nuclear_genome.tidy.fa
AGY731.asm01.HP0.nuclear_genome.tidy.fa
AHG.asm01.HP0.nuclear_genome.tidy.fa
AIC.asm01.HP0.nuclear_genome.tidy.fa
AIE.asm01.HP0.nuclear_genome.tidy.fa
ALS_1a.asm01.HP0.nuclear_genome.tidy.fa
ANM.asm01.HP.nuclear_genome.tidy.fa
BAI_1a.asm01.HP0.nuclear_genome.tidy.fa
BAK_1a.asm01.HP0.nuclear_genome.tidy.fa
BAW6.asm01.HP0.nuclear_genome.tidy.fa
BBM_1a.asm01.HP0.nuclear_genome.tidy.fa
BBT.asm01.HP.nuclear_genome.tidy.fa
BC187.asm01.HP0.nuclear_genome.tidy.fa
BDC_1c.asm01.HP0.nuclear_genome.tidy.fa
biocodex.asm01.HP0.nuclear_genome.tidy.fa
BJ4.asm01.HP0.nuclear_genome.tidy.fa
BY.asm01.HP0.nuclear_genome.tidy.fa
CAS_1a.asm01.HP0.nuclear_genome.tidy.fa
CBK.asm01.HP0.nuclear_genome.tidy.fa
CFC.asm01.HP.nuclear_genome.tidy.fa
CPA_1a.asm01.HP0.nuclear_genome.tidy.fa
CPG_1a.asm01.HP0.nuclear_genome.tidy.fa
CRL.asm01.HP0.nuclear_genome.tidy.fa
unique28.asm01.HP0.nuclear_genome.tidy.fa
Y12.asm01.HP0.nuclear_genome.tidy.fa
```

其中有些基因组文件不符合后续分析要求（按染色体排列），将在后续被剔除。

2. paf2dot文件

这个脚本后续分析要用，但我没在系统中找到它。故可以从以下方式下载。

```{shell}
## 下载paf2dot
wget https://raw.githubusercontent.com/waveygang/wfmash/refs/heads/main/scripts/paf2dotplot
```


### II.提取所有基因组特定染色体

1.  提取特定染色体序列

在此，我选择提取所有基因组的第一条染色体，即chrI。若基因组序列文件中没有找到chrI的序列，则该基因组不用。

```{shell}
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
```

2. 合并所有的染色体序列

由于在上一步已经将不同基因组来源的染色体标签命名好了，只需将上一步提取的所有染色体序列合并即可。

```{shell}
## combine the *chrI.fa and index it.
mkdir raw_chrI
mkdir out_chrI
mv *chrI.fa ./raw_chrI
cat ./raw_chrI/*.fa > ./out_chrI/total_chrI.fa
samtools faidx ./out_chrI/total_chrI.fa

num=$( grep ">.*" ./out_chrI/total_chrI.fa | wc -l )
echo "total genomes number: ${num}"
```
终端输出可用的基因组数目为**26** 。

### III. Run *pggb* and related visualization

相关参赛的含义：

- *-i*  input FASTA/FASTQ file
- *-o*  output directory
- *-n*  number of haplotypes
- *-t*  number of compute threads to use in parallel steps
- *-p*  percent identity for mapping/alignment

```{shell}
## run pggb
mkdir out_pggb
pggb -i total_chrI.fa -o out_pggb -n ${num} -t 5 -p 90

## dotplot
./paf2dotplot png large ./out_pggb/*.wfmash.paf

## pangenome graph statistics
odgi stats -i ./out_pggb/*.smooth.final.og -S > ./out_pggb/out_pggb_chrI.fa.smooth.final.og.stats
```


![FIG.2](/home/ubuntu/Desktop/advanced_bioinfo/raw_material_pang(副本)/out_chrI/out_pggb/total_chrI.fa.bf3285f.11fba48.78baeea.smooth.final.og.lay.draw_multiqc.png "Fig.2 2D visualization of Pangenome")


![FIG.3](/home/ubuntu/Desktop/advanced_bioinfo/raw_material_pang(副本)/out_chrI/out_pggb/total_chrI.fa.bf3285f.11fba48.78baeea.smooth.final.og.viz_multiqc.png "Fig.3 1D visualization of Pangenome")


![FIG.4](/home/ubuntu/Desktop/advanced_bioinfo/raw_material_pang(副本)/out_chrI/out_pggb/total_chrI.fa.bf3285f.11fba48.78baeea.smooth.final.og.viz_O_multiqc.png "Fig.4 0D visualization of Pangenome")


从上图可以看出，基因组一号染色体的**前$\frac{1}{2}$区段结构变异较大，形成的泛基因组图结构较为复杂**，分叉较多；而在一号染色体**后半部分（除了尾部），不同菌株之间差异较小，形成的泛基因组图结构较简单**。

为了更好地反应这一特点，利用vg分别提取染色体前后50bp的序列（相对于AIC.asm01.HP0序列450:500，和35450:35500）的区域，进行可视化：

```{shell}
## vg dive into pangenome structure
# 450:500
vg chunk -x ./out_pggb/out_pggb_chrI.smooth.final.vg -p AIC.asm01.HP0#1#chrI#0:450-500 -c 1 | vg view -pd - | dot -Tpng > ./out_vg/out_pggb_chrI.smooth.final.AIC_chrI_450-500_real.png
# 34500:35500
vg chunk -x ./out_pggb/out_pggb_chrI.smooth.final.vg -p AIC.asm01.HP0#1#chrI#0:35450-35500 -c 1 | vg view -pd - | dot -Tpng > ./out_vg/out_pggb_chrI.smooth.final.AIC_chrI_35450-35500_real.png
```
> vg chunk -p path究竟是什么我找了很久没找到:( 报错提示是在xg index里，但我也没有找到相关的xg index <br>不过可以先通过其他参数得到一张图(e.g.通过提取结点 -r 50:55)，在图中的最右边即是各个path的名字。


![FIG.5](/home/ubuntu/Desktop/advanced_bioinfo/raw_material_pang(副本)/out_chrI/out_vg/out_pggb_chrI.smooth.final.AIC_chrI_450-500_real.png "Fig.5 450:500")

![FIG.6](/home/ubuntu/Desktop/advanced_bioinfo/raw_material_pang(副本)/out_chrI/out_vg/out_pggb_chrI.smooth.final.AIC_chrI_35450-35500_real.png "Fig.6 35450:35500")

从上两图可以看出，450：500bp的泛基因组图结构明显比35450：35000bp复杂许多，即不同菌株基因组之间差异较大。

### IV. Calculate pangenome openness with *panacus*

1. panacus histgrowth

- *-l*  Ignore all countables with a coverage lower than the specified
 threshold. The coverage of a countable corresponds to the number of path/walk that contain it. Repeated appearances of a countable in the same path/walk are counted as one.

- *-q*  Unlike the -l parameter, which specifies a minimum constant number of paths for all growth point m (1 <= m <= num_paths), --quorum adjust the threshold based on m. At each m, a countable is counted in the average growth if the countable is contained in at least *floor(m\*quorum)* paths.

- *-S*  Merge counts from paths belonging to same sample

- *-a*  Also include histogram in output

2. panacus-visualize

- *-e*  Estimate growth parameters based on least-squares fit (default: False)
- *--split_subfigures*  Split output into multiple files (default: False)



```{shell}
## run panacus histgrowth to calculate coverage and pangenome growth for nodes (default) with coverage/quorum thresholds 1/0, 2/0, 1/1, 1/0.5, and 1/0.1
RUST_LOG=info panacus histgrowth -t4 -l 1,2,1,1,1 -q 0,0,1,0.5,0.1 -S -a \
    ./out_pggb/total_chrI.fa.bf3285f.11fba48.78baeea.smooth.final.gfa \
    > ./out_pggb/total_chrI.fa.bf3285f.11fba48.78baeea.smooth.final.histgrowth.node.tsv

## visualize the growth statistics
panacus-visualize \
    -e ./out_pggb/total_chrI.fa.bf3285f.11fba48.78baeea.smooth.final.histgrowth.node.tsv --split_subfigures -f png --split_prefix chrI_panacus_ -s 20 12

```

![FIG.7](/home/ubuntu/Desktop/advanced_bioinfo/raw_material_pang(副本)/out_chrI/chrI_panacus_0_0.png "Fig.7 panacus_1")

![FIG.8](/home/ubuntu/Desktop/advanced_bioinfo/raw_material_pang(副本)/out_chrI/chrI_panacus_0_1.png "Fig.8 panacus_2")

![FIG.9](/home/ubuntu/Desktop/advanced_bioinfo/raw_material_pang(副本)/out_chrI/chrI_panacus_0_2.png "Fig.9 panacus_3")

从**Fig.9**可以看出，随着基因组的添加，泛基因组的图结点增加数越来越少。然而，从**Fig.8**曲线拟合得到的$\gamma=0.216$并不是一个很小的数值；这可能是因为不同菌株间差异较大，使得该泛基因组闭合程度尚不很高。

---

下面是实现整个分析流程的 **shell script**:

```{shell}

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
```