setwd("C:/Users/86189/Desktop/neo_project/RNA")

library(readr)
rna0 <- read_csv("../OmicsExpressionProteinCodingGenesTPMLogp1.csv")
rna1 <- read_csv("CCLE_depMap_19Q1_TPM.csv")
rna2 <- read_csv("OmicsExpressionProteinCodingGenesTPMLogp1 (3).csv")
rna3 <- read_csv("CCLE_expression.csv")
rna4 <- read_csv("CCLE_expression (1).csv")
rna5 <- read_csv("CCLE_expression (2).csv")

meta=read.csv("../Model.csv", header = T, sep = ",")
ID2NAME=meta[,c(1,4)];colnames(ID2NAME) <- c("...1", "cell_line")

rna_list=list(rna0,rna1,rna3,rna4,rna5)
for (i in 1:length(rna_list)){
  rna_list[[i]] <- as.data.frame((rna_list[[i]]))
  rna_list[[i]]=merge(rna_list[[i]], ID2NAME, by.x="...1")
  rownames(rna_list[[i]]) <- rna_list[[i]]$cell_line
  rna_list[[i]]=subset(rna_list[[i]], select = -c(cell_line, ...1))
  print(i)
}


LOI=c("A549","K562","MCF7","UO31")
LOI_list=list()
for (i in 1:length(rna_list)) {
  LOI_list[[i]] <- rna_list[[i]][rownames(rna_list[[i]]) %in% LOI,]
  colnames(LOI_list[[i]]) <- gsub("\\([0-9]*\\)", "", colnames(LOI_list[[i]]))
  LOI_list[[i]] <- t(LOI_list[[i]])
  rownames(LOI_list[[i]])=gsub(" ","",rownames(LOI_list[[i]]))
  rownames(LOI_list[[i]])=gsub("or.*","",rownames(LOI_list[[i]]))
  rownames(LOI_list[[i]])=gsub("\\(ENSG.*","",rownames(LOI_list[[i]]))
}
gc()

A549_flux_list <- list()
K562_flux_list <- list()
MCF7_flux_list <- list()
UO31_flux_list <- list()
library(METAFlux)
for (i in 1:length(LOI_list)) {
  A549_scores <- calculate_reaction_score(as.matrix(LOI_list[[i]][,1]))
  A549_flux <- compute_flux(mras = A549_scores, medium = human_blood)
  A549_flux_list[[i]] <- A549_flux
  gc()
  K562_scores <- calculate_reaction_score(as.matrix(LOI_list[[i]][,2]))
  K562_flux <- compute_flux(mras = K562_scores, medium = human_blood)
  K562_flux_list[[i]] <- K562_flux
  gc()
  MCF7_scores <- calculate_reaction_score(as.matrix(LOI_list[[i]][,3]))
  MCF7_flux <- compute_flux(mras = MCF7_scores, medium = human_blood)
  MCF7_flux_list[[i]] <- MCF7_flux
  gc()
  UO31_scores <- calculate_reaction_score(as.matrix(LOI_list[[i]][,4]))
  UO31_flux <- compute_flux(mras = UO31_scores, medium = human_blood)
  UO31_flux_list[[i]] <- UO31_flux
  gc()
}


A549_flux_total <- as.data.frame(do.call(cbind, A549_flux_list));
colnames(A549_flux_total) <- c(1,2,3,4,5)
MCF7_flux_total <- as.data.frame(do.call(cbind, MCF7_flux_list));
colnames(MCF7_flux_total) <- c(1,2,3,4,5)
K562_flux_total <- as.data.frame(do.call(cbind, K562_flux_list));
colnames(K562_flux_total) <- c(1,2,3,4,5)
UO31_flux_total <- as.data.frame(do.call(cbind, UO31_flux_list));
colnames(UO31_flux_total) <- c(1,2,3,4,5)
rownames(A549_flux_total) <- rownames(A549_flux_list[[1]]);rm(A549_flux_list)
rownames(K562_flux_total) <- rownames(K562_flux_list[[1]]);rm(K562_flux_list)
rownames(MCF7_flux_total) <- rownames(MCF7_flux_list[[1]]);rm(MCF7_flux_list)
rownames(UO31_flux_total) <- rownames(UO31_flux_list[[1]]);rm(UO31_flux_list)
gc()
save(A549_flux_total,K562_flux_total,MCF7_flux_total,UO31_flux_total,file="flux_total.RData")

data("human_gem")
pathway_show(human_gem,A549_flux_total,"A549")
pathway_show(human_gem,K562_flux_total,"K562")
pathway_show(human_gem,MCF7_flux_total,"MCF7")
pathway_show(human_gem,UO31_flux_total,"UO31")

library(readxl)
rna_A549=read_xlsx("GSE206606_read_count_data.xlsx")
rna_A549=rna_A549[,c(1,14,15,16)];
rna_A549=as.data.frame(rna_A549)
rna_A549 <- rna_A549[!duplicated(rna_A549$Gene.names), ]
rownames(rna_A549)=rna_A549$Gene.names;
rna_A549=calculate_log2cpm(rna_A549)
rna_A549=rna_A549[,-1]

rna_A549_list <- list()
for (i in 1:ncol(rna_A549)){
  temp=as.matrix(rna_A549[,i]);rownames(temp)=rownames(rna_A549)
  rna_A549_scores <- calculate_reaction_score(temp)
  rna_A549_flux <- compute_flux(mras = rna_A549_scores, medium = human_blood)
  rna_A549_list[[i]] <- rna_A549_flux
  gc()
}

A549_flux <- as.data.frame(do.call(cbind, rna_A549_list));
pathway_show(human_gem,A549_flux,"A549")                        

A549_flux_total_total=cbind(A549_flux_total,A549_flux)
pathway_show(human_gem,A549_flux_total_total,"A549 more")





################## A549 more!
setwd("C:/Users/86189/Desktop/neo_project/RNA/new")
GENE_LENGTH=read.delim("C:/Users/86189/Desktop/neo_project/mart_export.txt")
GENE_LENGTH$length=abs(GENE_LENGTH$Gene.end..bp.-GENE_LENGTH$Gene.start..bp.)
GENE_LENGTH=subset(GENE_LENGTH,Gene.name!="")
gene_length=data.frame(gene_id=GENE_LENGTH$Gene.name,length=GENE_LENGTH$length)

library(readxl)
rna_A549_1=read_excel("GSE193930_RAW/GSM5823771_Control-24h-1_FPKM.xlsx")
rna_A549_1=subset(rna_A549_1,select=-c(Gene_ID))
rna_A549_1$`log2TPM1`=log2(fpkm_to_tpm(rna_A549_1$`Control-24h-1_FPKM`)+1)

rna_A549_1=as.data.frame(rna_A549_1)
rna_A549_1 <- rna_A549_1[!duplicated(rna_A549_1$`Gene Symbol`), ]
rownames(rna_A549_1)=rna_A549_1$`Gene Symbol`;rna_A549_1=subset(rna_A549_1,select=-c(`Gene Symbol`,`Control-24h-1_FPKM`))

rna_A549_2=read.csv("GSE290536_RNAseq_A549_raw_counts.csv/GSE290536_RNAseq_A549_raw_counts.csv")
rna_A549_2.1=rna_A549_2[,c("X","A549_1")];rna_A549_2.2=rna_A549_2[,c("X","A549_2")]
colnames(rna_A549_2.1)=c("gene_name","counts");colnames(rna_A549_2.2)=c("gene_name","counts");
rna_A549_2.1=calculate_log2tpm(rna_A549_2.1,gene_length);
rna_A549_2.2=calculate_log2tpm(rna_A549_2.2,gene_length);
rna_A549_2.1 <- rna_A549_2.1[!duplicated(rna_A549_2.1$`gene_name`), ]
rownames(rna_A549_2.1)=rna_A549_2.1$`gene_name`;rna_A549_2.1=subset(rna_A549_2.1,select=-c(`gene_name`))
rna_A549_2.2 <- rna_A549_2.2[!duplicated(rna_A549_2.2$`gene_name`), ]
rownames(rna_A549_2.2)=rna_A549_2.2$`gene_name`;rna_A549_2.2=subset(rna_A549_2.2,select=-c(`gene_name`))


rna_A549_3 <- read.delim("GSE211711_gene_expression.txt/GSE211711_gene_expression.txt", 
                   header = TRUE, 
                   stringsAsFactors = FALSE)
rna_A549_3=rna_A549_3[,c("gene_name","mock_count")];colnames(rna_A549_3)=c("gene_name","counts")
rna_A549_3=calculate_log2tpm(rna_A549_3,gene_length)
rna_A549_3 <- rna_A549_3[!duplicated(rna_A549_3$`gene_name`), ]
rownames(rna_A549_3)=rna_A549_3$`gene_name`;rna_A549_3=subset(rna_A549_3,select=-c(`gene_name`))


rna_A549_4 <- read.delim("GSE221796_RAW/GSM6896615_A549cntrl1_count.txt",header=TRUE)
colnames(rna_A549_4)=c("ENSEMBL","counts")
library(clusterProfiler)
library(org.Hs.eg.db)
gene_symbols <- bitr(
  rna_A549_4$Geneid,                 # 输入的 ENSEMBL ID
  fromType = "ENSEMBL",        # 输入 ID 的类型
  toType = "SYMBOL",           # 输出的字段
  OrgDb = org.Hs.eg.db         # 使用的数据库
)
rna_A549_4=merge(rna_A549_4,gene_symbols,by="ENSEMBL")
rna_A549_4=subset(rna_A549_4,select=-c(ENSEMBL));colnames(rna_A549_4)=c("counts","gene_name")
rna_A549_4=calculate_log2tpm(rna_A549_4,gene_length)
rna_A549_4 <- rna_A549_4[!duplicated(rna_A549_4$`gene_name`), ]
rownames(rna_A549_4)=rna_A549_4$`gene_name`;rna_A549_4=subset(rna_A549_4,select=-c(`gene_name`))


rna_A549_5 <- read.delim("GSE229282_RAW/GSM7157943_A549_Mock.txt")
rna_A549_5=rna_A549_5[,c("symbol","A549_Mock")];
rna_A549_5=subset(rna_A549_5,symbol!="")
rna_A549_5 <- rna_A549_5[!duplicated(rna_A549_5$`symbol`), ]
rownames(rna_A549_5)=rna_A549_5$`symbol`;rna_A549_5=subset(rna_A549_5,select=-c(`symbol`))


rna_A549_6 <- read.delim("GSE261915_A549_sgTNFAIP3_vs_sgCtrl_counts.txt/GSE261915_A549_sgTNFAIP3_vs_sgCtrl_counts.txt")
rna_A549_6.1=rna_A549_6[,c("Gene","A549.sgCtrl.1")];colnames(rna_A549_6.1)=c("gene_name","counts")
rna_A549_6.2=rna_A549_6[,c("Gene","A549.sgCtrl.2")];colnames(rna_A549_6.2)=c("gene_name","counts")
rna_A549_6.3=rna_A549_6[,c("Gene","A549.sgCtrl.3")];colnames(rna_A549_6.3)=c("gene_name","counts")
rna_A549_6.1=calculate_log2tpm(rna_A549_6.1,gene_length);
rna_A549_6.2=calculate_log2tpm(rna_A549_6.2,gene_length);
rna_A549_6.3=calculate_log2tpm(rna_A549_6.3,gene_length);
rna_A549_6.1 <- rna_A549_6.1[!duplicated(rna_A549_6.1$`gene_name`), ]
rownames(rna_A549_6.1)=rna_A549_6.1$`gene_name`;rna_A549_6.1=subset(rna_A549_6.1,select=-c(`gene_name`))
rna_A549_6.2 <- rna_A549_6.2[!duplicated(rna_A549_6.2$`gene_name`), ]
rownames(rna_A549_6.2)=rna_A549_6.2$`gene_name`;rna_A549_6.2=subset(rna_A549_6.2,select=-c(`gene_name`))
rna_A549_6.3 <- rna_A549_6.3[!duplicated(rna_A549_6.3$`gene_name`), ]
rownames(rna_A549_6.3)=rna_A549_6.3$`gene_name`;rna_A549_6.3=subset(rna_A549_6.3,select=-c(`gene_name`))

rna_A549_7 <- read.delim("GSE265781_RAW/GSM8229009_A549normoxic1_normalized.txt")
rna_A549_7=rna_A549_7[,c("Gene.Symbol","A5_N1")];colnames(rna_A549_7)=c("gene_name","counts")
rna_A549_7=calculate_log2tpm(rna_A549_7,gene_length)
rna_A549_7 <- rna_A549_7[!duplicated(rna_A549_7$`gene_name`), ]
rownames(rna_A549_7)=rna_A549_7$`gene_name`;rna_A549_7=subset(rna_A549_7,select=-c(`gene_name`))


rna_A549_more_list=list(rna_A549_1,rna_A549_2.1,rna_A549_2.2,rna_A549_3,rna_A549_4,rna_A549_5,rna_A549_6.1,rna_A549_6.2,rna_A549_6.3,rna_A549_7)
rna_A549_flux_list <- list()
library(METAFlux)
for (i in 1:length(rna_A549_more_list)){
  rna_A549_scores <- calculate_reaction_score(rna_A549_more_list[[i]])
  rna_A549_flux <- compute_flux(mras = rna_A549_scores, medium = human_blood)
  rna_A549_flux_list[[i]] <- rna_A549_flux
  gc()
}

A549_flux_more <- as.data.frame(do.call(cbind, rna_A549_flux_list));
pathway_show(human_gem,A549_flux_more,"A549")                        

A549_flux_total_total_total=cbind(A549_flux_total_total,A549_flux_more)
pathway_show(human_gem,A549_flux_total_total_total,"A549")



################## K562 more!
setwd("C:/Users/86189/Desktop/neo_project/RNA/new_new")
rna_K562_1 <- read.delim("GSE188415_K562.Ctrl_and_LEO1_KO.fragment_count_table.txt", 
                         header = TRUE, 
                         stringsAsFactors = FALSE)
rna_K562_1.1=rna_K562_1[,c("Geneid","Ctrl_4")];colnames(rna_K562_1.1)=c("gene_name","counts")
rna_K562_1.2=rna_K562_1[,c("Geneid","Ctrl_5")];colnames(rna_K562_1.2)=c("gene_name","counts")
rna_K562_1.3=rna_K562_1[,c("Geneid","Ctrl_6")];colnames(rna_K562_1.3)=c("gene_name","counts")
rna_K562_1.1=calculate_log2tpm(rna_K562_1.1,gene_length);
rna_K562_1.1 <- rna_K562_1.1[!duplicated(rna_K562_1.1$`gene_name`), ];rownames(rna_K562_1.1)=rna_K562_1.1$`gene_name`;rna_K562_1.1=subset(rna_K562_1.1,select=-c(`gene_name`))
rna_K562_1.2=calculate_log2tpm(rna_K562_1.2,gene_length);
rna_K562_1.2 <- rna_K562_1.2[!duplicated(rna_K562_1.2$`gene_name`), ];rownames(rna_K562_1.2)=rna_K562_1.2$`gene_name`;rna_K562_1.2=subset(rna_K562_1.2,select=-c(`gene_name`))
rna_K562_1.3=calculate_log2tpm(rna_K562_1.3,gene_length);
rna_K562_1.3 <- rna_K562_1.3[!duplicated(rna_K562_1.3$`gene_name`), ];rownames(rna_K562_1.3)=rna_K562_1.3$`gene_name`;rna_K562_1.3=subset(rna_K562_1.3,select=-c(`gene_name`))





rna_K562_2 <- read.delim("GSE188415_K562.Ctrl_and_PNUTS_KO.fragment_count_table.txt", 
                         header = TRUE, 
                         stringsAsFactors = FALSE)
rna_K562_2.1=rna_K562_2[,c("Geneid","Ctrl_1")];colnames(rna_K562_2.1)=c("gene_name","counts")
rna_K562_2.2=rna_K562_2[,c("Geneid","Ctrl_2")];colnames(rna_K562_2.2)=c("gene_name","counts")
rna_K562_2.3=rna_K562_2[,c("Geneid","Ctrl_3")];colnames(rna_K562_2.3)=c("gene_name","counts")
rna_K562_2.1=calculate_log2tpm(rna_K562_2.1,gene_length);
rna_K562_2.1 <- rna_K562_2.1[!duplicated(rna_K562_2.1$`gene_name`), ];rownames(rna_K562_2.1)=rna_K562_2.1$`gene_name`;rna_K562_2.1=subset(rna_K562_2.1,select=-c(`gene_name`))
rna_K562_2.2=calculate_log2tpm(rna_K562_2.2,gene_length);
rna_K562_2.2 <- rna_K562_2.2[!duplicated(rna_K562_2.2$`gene_name`), ];rownames(rna_K562_2.2)=rna_K562_2.2$`gene_name`;rna_K562_2.2=subset(rna_K562_2.2,select=-c(`gene_name`))
rna_K562_2.3=calculate_log2tpm(rna_K562_2.3,gene_length);
rna_K562_2.3 <- rna_K562_2.3[!duplicated(rna_K562_2.3$`gene_name`), ];rownames(rna_K562_2.3)=rna_K562_2.3$`gene_name`;rna_K562_2.3=subset(rna_K562_2.3,select=-c(`gene_name`))


rna_K562_3 <- read.delim("GSE278906_raw_counts.txt", 
                         header = TRUE, 
                         stringsAsFactors = FALSE)
rna_K562_3.1=rna_K562_3[,c("Gene","ScrambleIMA1")];colnames(rna_K562_3.1)=c("gene_name","counts")
rna_K562_3.2=rna_K562_3[,c("Gene","ScrambleIMA3")];colnames(rna_K562_3.2)=c("gene_name","counts")
rna_K562_3.3=rna_K562_3[,c("Gene","ScrambleIMA2")];colnames(rna_K562_3.3)=c("gene_name","counts")
rna_K562_3.1=calculate_log2tpm(rna_K562_3.1,gene_length);
rna_K562_3.1 <- rna_K562_3.1[!duplicated(rna_K562_3.1$`gene_name`), ];rownames(rna_K562_3.1)=rna_K562_3.1$`gene_name`;rna_K562_3.1=subset(rna_K562_3.1,select=-c(`gene_name`))
rna_K562_3.2=calculate_log2tpm(rna_K562_3.2,gene_length);
rna_K562_3.2 <- rna_K562_3.2[!duplicated(rna_K562_3.2$`gene_name`), ];rownames(rna_K562_3.2)=rna_K562_3.2$`gene_name`;rna_K562_3.2=subset(rna_K562_3.2,select=-c(`gene_name`))
rna_K562_3.3=calculate_log2tpm(rna_K562_3.3,gene_length);
rna_K562_3.3 <- rna_K562_3.3[!duplicated(rna_K562_3.3$`gene_name`), ];rownames(rna_K562_3.3)=rna_K562_3.3$`gene_name`;rna_K562_3.3=subset(rna_K562_3.3,select=-c(`gene_name`))




rna_K562_4=read.delim("GSE280672_K562_CRISPRi_salmon.merged.gene_counts.txt", 
                      header = TRUE, 
                      stringsAsFactors = FALSE)
rna_K562_4=rna_K562_4[,c("gene_name","K562.CRISPRi.sgControl")];colnames(rna_K562_4)=c("gene_name","counts")
rna_K562_4$counts=gsub(",","",rna_K562_4$counts);rna_K562_4$counts=as.numeric(rna_K562_4$counts)
rna_K562_4=calculate_log2tpm(rna_K562_4,gene_length);
rna_K562_4 <- rna_K562_4[!duplicated(rna_K562_4$`gene_name`), ];rownames(rna_K562_4)=rna_K562_4$`gene_name`;rna_K562_4=subset(rna_K562_4,select=-c(`gene_name`))




rna_K562_5=read.delim("GSE280885_FPKM_anno.txt", 
                      header = TRUE, 
                      stringsAsFactors = FALSE)
rna_K562_5.1=rna_K562_5[,c("id","CON1")];rna_K562_5.1$CON1=log2(fpkm_to_tpm(rna_K562_5.1$CON1)+1)
rna_K562_5.2=rna_K562_5[,c("id","CON2")];rna_K562_5.2$CON2=log2(fpkm_to_tpm(rna_K562_5.2$CON2)+1)
rna_K562_5.3=rna_K562_5[,c("id","CON3")];rna_K562_5.3$CON3=log2(fpkm_to_tpm(rna_K562_5.3$CON3)+1)
colnames(rna_K562_5.1)=c("gene_name","counts");colnames(rna_K562_5.2)=c("gene_name","counts");colnames(rna_K562_5.3)=c("gene_name","counts")

rna_K562_5.1 <- rna_K562_5.1[!duplicated(rna_K562_5.1$`gene_name`), ];rownames(rna_K562_5.1)=rna_K562_5.1$`gene_name`;rna_K562_5.1=subset(rna_K562_5.1,select=-c(`gene_name`))
rna_K562_5.2 <- rna_K562_5.2[!duplicated(rna_K562_5.2$`gene_name`), ];rownames(rna_K562_5.2)=rna_K562_5.2$`gene_name`;rna_K562_5.2=subset(rna_K562_5.2,select=-c(`gene_name`))
rna_K562_5.3 <- rna_K562_5.3[!duplicated(rna_K562_5.3$`gene_name`), ];rownames(rna_K562_5.3)=rna_K562_5.3$`gene_name`;rna_K562_5.3=subset(rna_K562_5.3,select=-c(`gene_name`))


rna_K562_6=read.csv("GSE285084_CYTOR_KD_Gene_Counts.csv", 
                      header = TRUE, 
                      stringsAsFactors = FALSE)
rna_K562_6$Geneid=gsub("\\.[0-9]*","",rna_K562_6$Geneid)
library(clusterProfiler)
library(org.Hs.eg.db)
gene_symbols <- bitr(
  rna_K562_6$Geneid,                 # 输入的 ENSEMBL ID
  fromType = "ENSEMBL",        # 输入 ID 的类型
  toType = "SYMBOL",           # 输出的字段
  OrgDb = org.Hs.eg.db         # 使用的数据库
)
rna_K562_6=merge(rna_K562_6,gene_symbols,by.x="Geneid",by.y="ENSEMBL")
rna_K562_6.1=rna_K562_6[,c("SYMBOL","NC100102")];colnames(rna_K562_6.1)=c("gene_name","counts")
rna_K562_6.2=rna_K562_6[,c("SYMBOL","NC100103")];colnames(rna_K562_6.2)=c("gene_name","counts")
rna_K562_6.3=rna_K562_6[,c("SYMBOL","NC100104")];colnames(rna_K562_6.3)=c("gene_name","counts")

rna_K562_6.3=calculate_log2tpm(rna_K562_6.3,gene_length);
rna_K562_6.3 <- rna_K562_6.3[!duplicated(rna_K562_6.3$`gene_name`), ];rownames(rna_K562_6.3)=rna_K562_6.3$`gene_name`;rna_K562_6.3=subset(rna_K562_6.3,select=-c(`gene_name`))
rna_K562_6.2=calculate_log2tpm(rna_K562_6.2,gene_length);
rna_K562_6.2 <- rna_K562_6.2[!duplicated(rna_K562_6.2$`gene_name`), ];rownames(rna_K562_6.2)=rna_K562_6.2$`gene_name`;rna_K562_6.2=subset(rna_K562_6.2,select=-c(`gene_name`))
rna_K562_6.1=calculate_log2tpm(rna_K562_6.1,gene_length);
rna_K562_6.1 <- rna_K562_6.1[!duplicated(rna_K562_6.1$`gene_name`), ];rownames(rna_K562_6.1)=rna_K562_6.1$`gene_name`;rna_K562_6.1=subset(rna_K562_6.1,select=-c(`gene_name`))


rna_K562_more_list=list(rna_K562_1.1,rna_K562_1.2,rna_K562_1.3,rna_K562_2.1,rna_K562_2.2,rna_K562_2.3,rna_K562_3.1,rna_K562_3.2,rna_K562_3.3,rna_K562_4,rna_K562_5.1,rna_K562_5.2,rna_K562_5.3,rna_K562_6.1,rna_K562_6.2,rna_K562_6.3)
rna_K562_flux_list <- list()
library(METAFlux)
for (i in 1:length(rna_K562_more_list)){
  rna_K562_scores <- calculate_reaction_score(rna_K562_more_list[[i]])
  rna_K562_flux <- compute_flux(mras = rna_K562_scores, medium = human_blood)
  rna_K562_flux_list[[i]] <- rna_K562_flux
  gc()
}

K562_flux_more <- as.data.frame(do.call(cbind, rna_K562_flux_list));
pathway_show(human_gem,K562_flux_more,"K562")                        

K562_flux_total_total=cbind(K562_flux_total,K562_flux_more)
pathway_show(human_gem,K562_flux_total_total,"K562")






