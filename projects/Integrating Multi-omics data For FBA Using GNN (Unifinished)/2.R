setwd("C:/Users/86189/Desktop/neo_project")

library(METAFlux)
library(ggplot2)
library(tidyr)
# BiocManager::install("biomaRt")
# library(biomaRt)
library(org.Hs.eg.db)
# library(AnnotationDbi)
library(clusterProfiler)

# Load the data
epi_lung= read.table("output_E114.txt", header=FALSE, sep="\t");colnames(epi_lung)=c("ENSEMBL","Prob")
temp=sprintf("%011d", epi_lung$ENSEMBL);temp=paste("ENSG",temp,sep="");epi_lung$ENSEMBL=temp
epi_leukemia= read.table("output_E123.txt", header=FALSE, sep="\t");colnames(epi_leukemia)=c("ENSEMBL","Prob")
temp=sprintf("%011d", epi_leukemia$ENSEMBL);temp=paste("ENSG",temp,sep="");epi_leukemia$ENSEMBL=temp

# Get the gene symbols
result1 <- bitr(epi_lung$ENSEMBL, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
result2 <- bitr(epi_leukemia$ENSEMBL, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)

# Merge
epi_lung=merge(epi_lung,result1,by.x="ENSEMBL");rm(result1)
epi_leukemia=merge(epi_leukemia,result2,by.x="ENSEMBL");rm(result2)

# filter
epi_lung=epi_lung[!duplicated(epi_lung$SYMBOL), ]
epi_leukemia=epi_leukemia[!duplicated(epi_leukemia$SYMBOL), ]

## combine with RNA-seq
mini_meta=matrix(c("E114","A549","Lung","E123","K562","Leukemia"),ncol=3,byrow=TRUE)
colnames(mini_meta)=c("ID","cell_line","Tissue")

rna_lung=t(ALL_LINE_NCI60["A549",]);rna_lung=as.data.frame(rna_lung);
rna_lung$SYMBOL=rownames(rna_lung);colnames(rna_lung)=c("abundence","SYMBOL")
rna_leukemia=t(ALL_LINE_NCI60["K562",]);rna_leukemia=as.data.frame(rna_leukemia);
rna_leukemia$SYMBOL=rownames(rna_leukemia);colnames(rna_leukemia)=c("abundence","SYMBOL")

comb_lung=merge(epi_lung,rna_lung,by="SYMBOL");#rownames(comb_lung)=comb_lung$SYMBOL
comb_leukemia=merge(epi_leukemia,rna_leukemia,by="SYMBOL");#rownames(comb_leukemia)=comb_leukemia$SYMBOL


# 归一化
comb_lung$abundence_normalized <- (comb_lung$abundence - min(comb_lung$abundence)) / (max(comb_lung$abundence) - min(comb_lung$abundence))
comb_leukemia$abundence_normalized <- (comb_leukemia$abundence - min(comb_leukemia$abundence)) / (max(comb_leukemia$abundence) - min(comb_leukemia$abundence))

# Plot
## histogram
ggplot(comb_lung, aes(x=abundence)) + geom_histogram(binwidth=0.1, fill="blue", color="black") + ggtitle("Lung") + xlab("RNA-seq abundence") + ylab("Frequency")

## scientific plot


p=ggplot(comb_leukemia, aes(x=abundence)) + 
  geom_histogram(binwidth=0.1, fill="blue", color="white") + 
  ggtitle("Leukemia(K562) RNA-seq abundence distribution") +
  xlab("RNA-seq abundence") + ylab("Frequency")+
  ## remove background
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "grey"),
        axis.text.x = element_text(size = 6, colour = "black"),
        axis.text.y = element_text(size = 6, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),
        plot.title = element_text(size = 12, colour = "black"),
        legend.position = "none")
ggsave("figures//leukemia_abundence.png", p, width = 8, height = 6, dpi = 300)

ggplot(comb_lung, aes(x=abundence_normalized)) + geom_histogram(binwidth=0.01, fill="blue", color="black") + ggtitle("Lung_normalized") + xlab("RNA-seq normalized") + ylab("Frequency")
ggplot(comb_leukemia, aes(x=abundence_normalized)) + geom_histogram(binwidth=0.01, fill="blue", color="black") + ggtitle("Leukemia_normalized") + xlab("RNA-seq normalized") + ylab("Frequency")

p=ggplot(comb_leukemia, aes(x=Prob)) + 
  geom_histogram(binwidth=0.01, fill="blue", color="white") +
  ggtitle("Leukemia(K562) gene expression probability distribution") + xlab("probability(output of AttentiveChrome)") + ylab("Frequency")+
  ## remove background
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "grey"),
        axis.text.x = element_text(size = 6, colour = "black"),
        axis.text.y = element_text(size = 6, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black"),
        plot.title = element_text(size = 12, colour = "black"),
        legend.position = "none")
ggsave("figures//leukemia_prob.png", p, width = 8, height = 6, dpi = 300)






ggplot(comb_leukemia, aes(x=Prob)) + geom_histogram(binwidth=0.005, fill="blue", color="black") + ggtitle("K562/E123") + xlab("probability") + ylab("Frequency")

## correlation
correlation_lung <- cor(comb_lung$abundence_normalized, comb_lung$Prob, method = "spearman")
ggplot(comb_lung, aes(x = abundence_normalized, y = Prob)) +
  geom_point(color = "blue", size = 1,alpha=0.1) +  # 绘制散点
  geom_smooth(method = "lm", color = "navy", se = TRUE) +  # 添加线性回归线
  annotate("text", x = 0.7*max(comb_lung$abundence_normalized), y = max(comb_lung$Prob), 
           label = paste("Spearman r =", round(correlation_lung, 3)), 
           hjust = 0, vjust = 1, color = "black", size = 3) +  # 标注相关性
  theme_minimal() +  # 使用简洁主题
  labs(title = "Scatter Plot with Correlation",
       x = "RNA abundence",
       y = "Prob")+ylim(0,1)

correlation_leukemia <- cor(comb_leukemia$abundence_normalized, comb_leukemia$Prob, method = "spearman")
ggplot(comb_leukemia, aes(x = abundence_normalized, y = Prob)) +
  geom_point(color = "blue", size = 1,alpha=0.1) +  # 绘制散点
  geom_smooth(method = "lm", color = "navy", se = TRUE) +  # 添加线性回归线
  annotate("text", x = 0.7*max(comb_leukemia$abundence_normalized), y = max(comb_leukemia$Prob), 
           label = paste("Spearman r =", round(correlation_leukemia, 3)), 
           hjust = 0, vjust = 1, color = "black", size = 3) +  # 标注相关性
  theme_minimal() +  # 使用简洁主题
  labs(title = "Scatter Plot with Correlation",
       x = "RNA abundence_normalized",
       y = "Prob")+ylim(0,1)



################ only epi as input
input_lung=as.data.frame(epi_lung$Prob);
rownames(input_lung)=epi_lung$SYMBOL
scores<-calculate_reaction_score(input_lung)

flux_lung<-compute_flux(mras=scores,medium=human_blood) 
biomass_flux_lung=t(obj) %*% flux_lung

input_leukemia=as.data.frame(epi_leukemia$Prob);
rownames(input_leukemia)=epi_leukemia$SYMBOL
scores<-calculate_reaction_score(input_leukemia)
flux_leukemia<-compute_flux(mras=scores,medium=human_blood)
biomass_flux_leukemia=t(obj) %*% flux_leukemia


################ insert some irrelevance
lung_rna=na.omit(total[total$cancer_type=="Non-Small Cell Lung",])$cell_line
flux_rna_lung=data.frame(reaction=rownames(flux_lung))
for (i in 1:length(lung_rna)){
  rna=as.data.frame(ALL_LINE_NCI60[lung_rna[i],]);rna=t(rna)
  scores<-calculate_reaction_score(rna)
  temp=compute_flux(mras=scores,medium=human_blood)
  flux_rna_lung[[as.character(lung_rna[i])]] = as.vector(temp)
  #gc()
}
rownames(flux_rna_lung)=flux_rna_lung$reaction;flux_rna_lung=flux_rna_lung[,-1]
pathway_show(human_gem,flux_rna_lung,"rna lung")

leukemia_rna=na.omit(total[total$cancer_type=="Leukemia",])$cell_line
flux_rna_leukemia=data.frame(reaction=rownames(temp))
for (i in 1:length(leukemia_rna)){
  rna=as.data.frame(ALL_LINE_NCI60[leukemia_rna[i],]);rna=t(rna)
  scores<-calculate_reaction_score(rna)
  temp=compute_flux(mras=scores,medium=human_blood)
  flux_rna_leukemia[[as.character(leukemia_rna[i])]] = as.vector(temp)
  gc()
}
rownames(flux_rna_leukemia)=flux_rna_leukemia$reaction;flux_rna_leukemia=flux_rna_leukemia[,-1]
pathway_show(human_gem,flux_rna_leukemia,"rna leukemia")




