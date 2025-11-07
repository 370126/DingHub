setwd("C://Users//86189//Desktop//neo_project")
genes_GEM=read.table("GEM//genes.tsv", header=T, sep="\t")
genes_GEM=genes_GEM$genes

################### transcriptome data preparation ###################
library(tidyr)
library(dplyr)
# rna_K562_list
# rna_A549_list
names(rna_A549_list) <- paste0("Sample", 1:18)
names(rna_K562_list) <- paste0("Sample", 1:18)



filtered_list <- lapply(rna_A549_list, function(df) {
  df[df$SYMBOL %in% genes_GEM, c(1,2)]
})
genes_check <- sapply(filtered_list, function(df) {
  all(genes_GEM %in% df$SYMBOL)
})
print(genes_check)  # TRUE 表示该样本中包含所有基因，FALSE 否则
named_list <- lapply(filtered_list, function(df) {
  df_out <- df[, c("SYMBOL", "abundence")]
  colnames(df_out)[2] <- "abundence"
  return(df_out)
})
for (i in seq_along(named_list)) {
  colnames(named_list[[i]])[2] <- paste0("Sample_", i)
}
merged_df <- named_list[[1]]
for (i in 2:length(named_list)) {
  merged_df <- merge(merged_df, named_list[[i]], by = "SYMBOL", all = TRUE)
}
rna_A549_matrix=merged_df
rm(filtered_list);rm(named_list);rm(merged_matrix);rm(merged_df)


result1 <- bitr(rna_A549_matrix$SYMBOL, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Hs.eg.db)
rna_A549_matrix=merge(rna_A549_matrix,result1,by.x="SYMBOL");rm(result1);
View(rna_A549_matrix);rna_A549_matrix=subset(rna_A549_matrix, select=-c(SYMBOL))
# let the first column be the ESGN
rna_A549_matrix=rna_A549_matrix[,c(ncol(rna_A549_matrix),1:(ncol(rna_A549_matrix)-1))]

# check the rows that are all NAs
temp=rna_A549_matrix[,2:ncol(rna_A549_matrix)]
na_rows <- apply(temp, 1, function(x) all(is.na(x)))
# get the row names of the rows that are all NAs
na_row_names <- rna_A549_matrix[na_rows, 1];na_row_names
# remove the rows that are all NAs
rna_A549_matrix <- rna_A549_matrix[!na_rows, ]
# save into .csv files
write.csv(rna_A549_matrix, file="python/rna_A549_matrix.csv", row.names=F)




################### epigenome data preparation ###################

E114=read.table("all_E114_uniq.tsv", header=F, sep="\t")
E123=read.table("all_E123_uniq.tsv", header=F, sep="\t")
colnames(E114)=c("ENSEMBL","position","H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3", "Label")
colnames(E123)=c("ENSEMBL","position","H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3", "Label")
temp=sprintf("%011d", E114$ENSEMBL);temp=paste("ENSG",temp,sep="");E114$ENSEMBL=temp
temp=sprintf("%011d", E123$ENSEMBL);temp=paste("ENSG",temp,sep="");E123$ENSEMBL=temp
gc()

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyr)
library(dplyr)
E114_list=list()
for (i in 1:5) {
  temp=E114[,c(1,2,i+2)]
  colnames(temp)=c("ENSEMBL","position",colnames(E114)[i+2])
  temp <- temp %>%
    pivot_wider(
      names_from = position,
      values_from = colnames(E114)[i+2],
      names_sort = TRUE,  # 列按 1~100 顺序排列
    )
  # result1 <- bitr(temp$ENSEMBL, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  # Merge
  # temp=merge(temp,result1,by.x="ENSEMBL");rm(result1)
  # temp=subset(temp, select=-c(ENSEMBL))
  # let the first column be the gene symbol and order gene symbol
  # temp=temp[,c(ncol(temp),1:(ncol(temp)-1))];
  # temp=temp[order(temp[,1]),]
  # filter genes only in genes_GEM
  temp=temp[temp$ENSEMBL %in% genes_GEM,]
  E114_list[[i]]=temp
  names(E114_list)[i]=colnames(E114)[i+2]
}





E123_list=list()
for (i in 1:5) {
  temp=E123[,c(1,2,i+2)]
  colnames(temp)=c("ENSEMBL","position",colnames(E123)[i+2])
  temp <- temp %>%
    pivot_wider(
      names_from = position,
      values_from = colnames(E123)[i+2],
      names_sort = TRUE,  # 列按 1~100 顺序排列
    )
  # result1 <- bitr(temp$ENSEMBL, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  # temp=merge(temp,result1,by.x="ENSEMBL");rm(result1)
  # temp=subset(temp, select=-c(ENSEMBL))
  # let the first column be the gene symbol and order gene symbol
  # temp=temp[,c(ncol(temp),1:(ncol(temp)-1))];temp=temp[order(temp[,1]),]
  # filter genes only in genes_GEM
  temp=temp[temp$ENSEMBL %in% genes_GEM,]
  E123_list[[i]]=temp
  names(E123_list)[i]=colnames(E123)[i+2]
}


# save into .csv files
for (i in 1:5) {
  write.csv(E114_list[[i]], file=paste("python/E114_", names(E114_list)[i], ".csv", sep=""), row.names=F)
  write.csv(E123_list[[i]], file=paste("python/E123_", names(E123_list)[i], ".csv", sep=""), row.names=F)
}



epi_lung= read.table("output_E114.txt", header=FALSE, sep="\t");colnames(epi_lung)=c("ENSEMBL","Prob")
temp=sprintf("%011d", epi_lung$ENSEMBL);temp=paste("ENSG",temp,sep="");epi_lung$ENSEMBL=temp
epi_leukemia= read.table("output_E123.txt", header=FALSE, sep="\t");colnames(epi_leukemia)=c("ENSEMBL","Prob")
temp=sprintf("%011d", epi_leukemia$ENSEMBL);temp=paste("ENSG",temp,sep="");epi_leukemia$ENSEMBL=temp
# min-max normalization
epi_lung$Prob=(epi_lung$Prob-min(epi_lung$Prob))/(max(epi_lung$Prob)-min(epi_lung$Prob))
epi_leukemia$Prob=(epi_leukemia$Prob-min(epi_leukemia$Prob))/(max(epi_leukemia$Prob)-min(epi_leukemia$Prob))

# save into .csv files
write.csv(epi_lung, file="python/epi_lung.csv", row.names=F)
write.csv(epi_leukemia, file="python/epi_leukemia.csv", row.names=F)








## average other columns based on first column
E114_ave=aggregate(E114[,2:ncol(E114)], by=list(E114[,1]), FUN=mean)
E123_ave=aggregate(E123[,2:ncol(E123)], by=list(E123[,1]), FUN=mean)

E114_ave=subset(E114_ave, select=-c(position, Label))
E123_ave=subset(E123_ave, select=-c(position, Label))

temp=sprintf("%011d", E114_ave$Group.1);temp=paste("ENSG",temp,sep="");E114_ave$ENSEMBL=temp
temp=sprintf("%011d", E123_ave$Group.1);temp=paste("ENSG",temp,sep="");E123_ave$ENSEMBL=temp

# Get the gene symbols
library(clusterProfiler)
library(org.Hs.eg.db)
result1 <- bitr(E114_ave$ENSEMBL, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db);result2 <- bitr(E123_ave$ENSEMBL, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)

# Merge
E114_ave=merge(E114_ave,result1,by.x="ENSEMBL");rm(result1)
E123_ave=merge(E123_ave,result2,by.x="ENSEMBL");rm(result2)

E114_ave=subset(E114_ave, select=-c(Group.1, ENSEMBL))
E123_ave=subset(E123_ave, select=-c(Group.1, ENSEMBL))

# let the first column be the gene symbol
E114_ave=E114_ave[,c(ncol(E114_ave),1:(ncol(E114_ave)-1))]
E123_ave=E123_ave[,c(ncol(E123_ave),1:(ncol(E123_ave)-1))]
E114_ave=E114_ave[order(E114_ave[,1]),]
E123_ave=E123_ave[order(E123_ave[,1]),]
# save dataframe into .csv files
write.csv(E114_ave, file="python/E114_ave.csv", row.names=F)
write.csv(E123_ave, file="python/E123_ave.csv", row.names=F)

#####

## make gene_name be the first column
epi_lung=epi_lung[,c(ncol(epi_lung),1:(ncol(epi_lung)-1))]
epi_lung=epi_lung[order(epi_lung[,1]),]
epi_lung$Prob=(epi_lung$Prob-min(epi_lung$Prob))/(max(epi_lung$Prob)-min(epi_lung$Prob))

epi_leukemia=epi_leukemia[,c(ncol(epi_leukemia),1:(ncol(epi_leukemia)-1))]
epi_leukemia=epi_leukemia[order(epi_leukemia[,1]),]
# save into .csv files
write.csv(epi_lung, file="python/epi_lung.csv", row.names=F)
write.csv(epi_leukemia, file="python/epi_leukemia.csv", row.names=F)
