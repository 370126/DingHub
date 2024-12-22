rm(list = ls())
options(stringAsFactors=F)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("tximport")
setwd("./salmon_q");getwd()
library(tximport)
library(tidyverse) # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
library(data.table) #多核读取文件

## data import
t2s <- fread("./transcrip_id_2_gene_id.txt", data.table = F, header = F); head(t2s)
files<-list.files(path="./",pattern="*quant_clean.sf",recursive=T,full.names=T);files
# 删掉SRR29422225,SRR29422229两个不好的样本
files <- files[!grepl("SRR29422225", files)]
files <- files[!grepl("SRR29422229", files)];files
txi_salmon <- tximport(files, type = "salmon", tx2gene = t2s)

cn <- sapply(strsplit(files,'\\/'), function(x) x[length(x)-1]); cn
colnames(txi_salmon$counts) <- gsub('_quant','',cn); colnames(txi_salmon$counts)

counts <- as.data.frame(apply(txi_salmon$counts,2,as.integer)) #将counts数取整
rownames(counts) <- rownames(txi_salmon$counts) 
tpm <- as.data.frame(txi_salmon$abundance)  ###abundance为基因的Tpm值
colnames(tpm) <- colnames(txi_salmon$counts)

#筛选出至少在重复样本数量内的表达量counts大于1的行（基因）
keep_feature <- rowSums(counts > 1) >= 2               #ncol(counts)/length(table(group_list)) 
table(keep_feature)  #查看筛选情况
counts_filt <- counts[keep_feature, ] #替换counts为筛选后的基因矩阵（保留较高表达量的基因）
tpm_filt <- tpm[keep_feature, ]


## 数据检查
library(FactoMineR)
library(factoextra)  
library(tidyverse) # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
library(pheatmap)
library(DESeq2)
library(RColorBrewer)

# 归一化
#dat <- as.data.frame(log2(edgeR::cpm(counts)+1))
dat <- log2(tpm_filt+1)

sample_table <- data.frame(
  treatment = as.factor(c(rep("HOX5aKO_somites", 2), rep("WT_somites", 2)))
);sample_table
rownames(sample_table)=gsub('_quant','',cn);rownames(sample_table)

# 各样本整体表达量
boxplot(dat, ylab="dat", main=" normalized data ",
        outline = F, notch = F)
# hClust
sampleDists <- dist(t(dat))   #dist默认计算矩阵行与行的距离， 因此需要转置
sampleDistMatrix <- as.matrix(sampleDists)  
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)  #选取热图的颜色
pheatmap::pheatmap(sampleDistMatrix,
                         fontsize=7,
                         clustering_distance_rows=sampleDists,
                         clustering_distance_cols=sampleDists,
                         angle_col=45,
                         col=colors)
# PCA
dat.pca <- PCA(t(dat),graph = F) 
fviz_pca_ind(dat.pca,
                    title = "Principal Component Analysis",
                    legend.title = "Groups",
                    geom.ind = c("point", "text"), 
                    pointsize = 1.5,
                    labelsize = 3,
                    col.ind = sample_table$treatment, # 分组上色
                    # axes.linetype=NA,  # remove axeslines
                    mean.point=F#去除分组中心点
) 

cg <- names(tail(sort(apply(dat,1,sd)),200)) #取每一行的方差，从小到大排序，取最大的500个
pheatmap::pheatmap(dat[cg, ],show_colnames =T,show_rownames = F,
                         fontsize=7,
                         # legend_breaks = -3:3,
                         scale = "row",
                         angle_col=45,
                        ) 

dat_500 <- dat[names(sort(apply(dat,1,mad),decreasing = T)[1:500]),]#取高表达量前500基因
M <- cor(dat_500)

pheatmap::pheatmap(M,
                        show_rownames = T,
                        angle_col=45,
                        fontsize=7)



## DESeq

library(DESeq2)
# browseVignettes("DESeq2")


# dds <- DESeqDataSetFromTximport(
#   txi=txi_salmon,
#   colData = sample_table,
#   design = ~treatment
# )


dds <- DESeqDataSetFromTximport(txi=txi_salmon, 
                              colData = sample_table, 
                              design = ~ treatment)
dds$treatment <- relevel(dds$treatment,ref="WT_somites")
dds <- DESeq(dds)
res <- results(dds)
#View(as.data.frame(res))
res <- as.data.frame(res)
res$pvalue[is.na(res$pvalue)]<- 1
res <- res[order(res$pvalue),]
gene_list=rownames(res)[res$pvalue<0.05]

## volcano plot
# cut-off: p_value<0.05, log2FoldChange>5
log2FC_cutoff=log2(5)
valc<- res[,c(2,5)]
valc$significance  <- as.factor(ifelse(valc$pvalue < 0.05 & abs(valc$log2FoldChange) > log2FC_cutoff,
                                           ifelse(valc$log2FoldChange > log2FC_cutoff ,'UP','DOWN'),'NOT'))
valc=as.data.frame(valc)
ggplot(data=valc, 
       aes(x=log2FoldChange, y=-log10(pvalue), 
           color=significance)) +
  #点和背景
  geom_point(alpha=0.2, size=1) +
  theme_classic()+ #无网格线
  #坐标轴
  xlab("log2 ( FoldChange )") + 
  ylab("-log10 ( P.value )") +
  #标题文本
  ggtitle("p-value v.s. log2FoldChange") +
  #分区颜色                  
  scale_colour_manual(values = c('blue','grey','red'))+ 
  #辅助线
  geom_vline(xintercept = c(-log2FC_cutoff,log2FC_cutoff),lty=4,col="grey",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.8) +
  coord_cartesian(xlim = c(-15, 15), ylim = c(0, 10))+
  #图例标题间距等设置
  theme(plot.title = element_text(hjust = 0.5), 
        plot.margin=unit(c(2,2,2,2),'lines'), #上右下左
        legend.title = element_blank(), #不显示图例标题
        legend.position="right")  #图例位置





## gene cluster analysis

keep <- rowSums(counts(dds)) >= 1.5*ncol(counts)  #Pre-filtering ，过滤低表达基因
dds <- dds[keep,] 
dds <- DESeq(dds,quiet = F) 
# res <- results(dds,contrast=c("sample_table", "HOX5aKO_somites", "WT_somites"))  #指定提取为exp/ctr结果
res<-results(dds)
resOrdered <- res[order(res$padj),]  #order根据padj从小到大排序结果
tempDEG <- as.data.frame(resOrdered)
DEG_DEseq2 <- na.omit(tempDEG)

need_DEG <- DEG_DEseq2[,c(2,5,6)]  
head(need_DEG)
colnames(need_DEG) <- c('log2FoldChange','pvalue','padj')
# 筛选条件设置
log2FC_cutoff = log2(5)
pvalue_cutoff = 0.05
gene_up=rownames(need_DEG[with(need_DEG,log2FoldChange > log2FC_cutoff & pvalue<pvalue_cutoff),])
gene_down=rownames(need_DEG[with(need_DEG,log2FoldChange < -log2FC_cutoff & pvalue<pvalue_cutoff),])



library(clusterProfiler)
# BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

gene_up_entrz<-bitr(gene_up,
           fromType="ENSEMBL",
           toType="ENTREZID",
           OrgDb = "org.Mm.eg.db")
gene_down_entrz<-bitr(gene_down,
                      fromType="ENSEMBL",
                      toType="ENTREZID",
                      OrgDb = "org.Mm.eg.db")


# KEGG clustering
kegg_enrich_down <- enrichKEGG(gene  = gene_down_entrz$ENTREZID,
                                  organism  = "mmu", #物种人类 hsa 小鼠mmu
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.2)
kegg_enrich_up <- enrichKEGG(gene  = gene_up_entrz$ENTREZID,
                               organism  = "mmu", #物种人类 hsa 小鼠mmu
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.2)
browseKEGG(kegg_enrich_down, 'mmu') #网页查看通路
browseKEGG(kegg_enrich_up, 'mmu') #网页查看通路
# GO clustering
ego_enrich_down=enrichGO(
  gene=gene_down,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont="ALL",
  readable = TRUE,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
  
)
ego_enrich_up=enrichGO(
  gene=gene_up,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont="ALL",
  readable = TRUE,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
  
)
# 绘制气泡图
dotplot(ego_enrich_down, showCategory = 10) +
  ggtitle("GO Enrichment Analysis") +
  theme_minimal()
gop <- goplot(go_enrich_down, showCategory = 10)

# GSEA 
library(ReactomePA)    
