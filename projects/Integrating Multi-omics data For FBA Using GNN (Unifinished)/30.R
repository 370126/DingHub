setwd("C:/Users/86189/Desktop/neo_project")
A549_flux_total_total_total
K562_flux_total_total

A549_prob=as.data.frame(comb_lung$Prob);A549_prob$gene_name=comb_lung$SYMBOL
K562_prob=as.data.frame(comb_leukemia$Prob);K562_prob$gene_name=comb_leukemia$SYMBOL

# 定义目标细胞系名称
cell_lines <- c("K562", "A549")
rna_data_list <- lapply(cell_lines, function(cl) {
  # 提取并转换为数据框
  as.data.frame(do.call(cbind, lapply(LOI_list, function(mat) mat[, cl])))
})
names(rna_data_list) <- cell_lines  # 命名列表元素
rna_K562_raw=rna_data_list$K562;rna_K562_raw$gene_name=rownames(rna_K562_raw)
rna_A549_raw=rna_data_list$A549;rna_A549_raw$gene_name=rownames(rna_A549_raw)

rna_K562_list=list()
rna_A549_list=list()
for (i in 1:(length(rna_A549_raw)-1)){
  rna_A549_list[[i]]=data.frame(gene_name=rna_A549_raw$gene_name,abundence=rna_A549_raw[,i])
}
tmp=length(rna_A549_list)
for (i in 1:(length(rna_A549))){
  rna_A549_list[[tmp+i]]=data.frame(gene_name=rownames(rna_A549),abundence=rna_A549[,i])
}

for (i in 1:(length(rna_K562_raw)-1)){
  rna_K562_list[[i]]=data.frame(gene_name=rna_K562_raw$gene_name,abundence=rna_K562_raw[,i])
}

library(data.table)
#转换行名为列并重命名
rna_K562_merged_list <- lapply(seq_along(rna_K562_more_list), function(i) {
  rna_K562_merged <- as.data.table(rna_K562_more_list[[i]], keep.rownames = "gene_name")
  setnames(rna_K562_merged, "counts", paste0("sample", i, "_abundence"))
  return(rna_K562_merged)
})
rna_A549_merged_list <- lapply(seq_along(rna_A549_more_list), function(i) {
  dt <- as.data.table(rna_A549_more_list[[i]], keep.rownames = "gene_name")
  original_col <- setdiff(names(dt), "gene_name")
  setnames(dt, original_col, paste0("sample", i, "_abundence"))
  return(dt)
})

rna_A549_list=c(rna_A549_list,rna_A549_merged_list)
rna_K562_list=c(rna_K562_list,rna_K562_merged_list);rna_K562_list=rna_K562_list[1:18]


epi_lung=subset(epi_lung,select = -c(ENSEMBL));colnames(epi_lung)=c("Prob","gene_name")
epi_leukemia=subset(epi_leukemia,select = -c(ENSEMBL));colnames(epi_leukemia)=c("Prob","gene_name")

## merge epi data
for (i in 1:length(rna_A549_list)){
  rna_A549_list[[i]]=merge(rna_A549_list[[i]],epi_lung,by="gene_name")
}

for (i in 1:length(rna_K562_list)){
  rna_K562_list[[i]]=merge(rna_K562_list[[i]],epi_leukemia,by="gene_name")
}


## calculate z-post
for (i in 1:length(rna_A549_list)){
  #names(rna_A549_list[[i]])[names(rna_A549_list[[i]]) == "gene_name"] <- "SYMBOL"
  z_score=Bayesian(rna_A549_list[[i]])
  rna_A549_list[[i]]$z_post=z_score
}

for (i in 1:length(rna_K562_list)){
  #names(rna_K562_list[[i]])[names(rna_K562_list[[i]]) == "gene_name"] <- "SYMBOL"
  z_score=Bayesian(rna_K562_list[[i]])
  rna_K562_list[[i]]$z_post=z_score
}

## filter
for (i in 1:length(rna_A549_list)){
  #rna_A549_list[[i]]=subset(rna_A549_list[[i]],select = -c(abundence_adjusted))
  names(rna_A549_list[[i]])[2] <- "abundence"
  #rna_A549_list[[i]]=filter_rna(rna_A549_list[[i]])
  rna_A549_list[[i]]=handle_epigenetic_transcriptome_conflict(rna_A549_list[[i]], prob_col = "Prob", abundence_col = "abundence", prob_high_threshold = 0.65, prob_low_threshold = 0.2)
}

for (i in 1:length(rna_K562_list)){
  names(rna_K562_list[[i]])[2] <- "abundence"
  #rna_K562_list[[i]]=filter_rna(rna_K562_list[[i]])
  rna_K562_list[[i]]=handle_epigenetic_transcriptome_conflict(rna_K562_list[[i]], prob_col = "Prob", abundence_col = "abundence", prob_high_threshold = 0.85, prob_low_threshold = 0.2)
}









############################ check the comflict genes

## scatter plot of conflict genes
p=ggplot(rna_K562_list[[1]], aes(x = abundence, y = Prob,colour = Conflict_Type)) +
  geom_point(alpha=0.3,size=1) +
  # set palette for scientific journal
  scale_color_manual(values = c("Type1" = "blue", "Type2" = "red", "None" = "black")) +
  labs(title = "K562 Conflict Genes") +
  theme_minimal()
ggsave("figures//K562_conflict_genes.png", p, width = 6, height = 4, dpi = 400)


## type2
A549_type2_list <- lapply(rna_A549_list, function(df) {
  df[df$Conflict_Type == "Type2", "SYMBOL", drop = TRUE]
});names(A549_type2_list) <- paste("Sample", 1:length(A549_type2_list))
K562_type2_list <- lapply(rna_K562_list, function(df) {
  df[df$Conflict_Type == "Type2", "SYMBOL", drop = TRUE]
});names(K562_type2_list) <- paste("Sample", 1:length(K562_type2_list))

# 提取高频SYMBOL

high_freqency=0.5*(length(rna_A549_list))


gene_counts=table(unlist(A549_type2_list))
high_freq_type2_genes_A549 <- names(gene_counts[gene_counts >= high_freqency]);high_freq_type2_genes_A549 <- as.character(high_freq_type2_genes_A549)
cat("共发现", length(high_freq_type2_genes_A549), "个高频基因：\n")

gene_counts=table(unlist(K562_type2_list))
high_freq_type2_genes_K562 <- names(gene_counts[gene_counts >= high_freqency]);high_freq_type2_genes_K562 <- as.character(high_freq_type2_genes_K562)
cat("共发现", length(high_freq_type2_genes_K562), "个高频基因：\n")

## 可视化
ggplot(data.frame(Count=as.numeric(table(unlist(A549_type2_list)))), aes(x=Count)) +
       geom_histogram(binwidth=1, fill="steelblue") +
       scale_x_continuous(breaks=1:18) +
       labs(x="Number of datasets containing gene",y="Number of genes",title="A549 type2 symbol co-occurence") +
       theme_minimal()
p=ggplot(data.frame(Count=as.numeric(table(unlist(K562_type2_list)))), aes(x=Count)) +
  geom_histogram(binwidth=1, fill="steelblue") +
  scale_x_continuous(breaks=1:18) +
  labs(x="Number of datasets containing gene",y="Number of genes",title="K562 type2 symbol co-occurence") +
  theme_minimal()
ggsave("figures//K562_type2_symbol_co-occurence.png", p, width = 6, height = 4, dpi = 400)

## pheatmap
show_co_occurence(A549_type2_list,"A549 type2 Gene Presence/Absence")
p=show_co_occurence(K562_type2_list,"K562 type2 Gene Presence/Absence")
ggsave("figures//K562_type2_gene_presence_absence.png", p, width = 6, height = 4, dpi = 400)

## type1
A549_type1_list <- lapply(rna_A549_list, function(df) {
  df[df$Conflict_Type == "Type1", "SYMBOL", drop = TRUE]
});names(A549_type1_list) <- paste("Sample", 1:length(A549_type1_list))
K562_type1_list <- lapply(rna_K562_list, function(df) {
  df[df$Conflict_Type == "Type1", "SYMBOL", drop = TRUE]
});names(K562_type1_list) <- paste("Sample", 1:length(K562_type1_list))

# 提取高频SYMBOL（出现次数 >=10 次）

gene_counts=table(unlist(A549_type1_list))

high_freq_type1_genes_A549 <- names(gene_counts[gene_counts >= 10]);high_freq_type1_genes_A549 <- as.character(high_freq_type1_genes_A549)
cat("共发现", length(high_freq_type1_genes_A549), "个高频基因：\n")

gene_counts=table(unlist(K562_type1_list))
high_freq_type1_genes_K562 <- names(gene_counts[gene_counts >= 10]);high_freq_type1_genes_K562 <- as.character(high_freq_type1_genes_K562)
cat("共发现", length(high_freq_type1_genes_K562), "个高频基因：\n")

## 可视化
ggplot(data.frame(Count=as.numeric(table(unlist(A549_type1_list)))), aes(x=Count)) +
  geom_histogram(binwidth=1, fill="steelblue") +
  scale_x_continuous(breaks=1:18) +
  labs(x="Number of datasets containing gene",y="Number of genes",title="A549 type1 symbol co-occurence") +
  theme_minimal()
ggplot(data.frame(Count=as.numeric(table(unlist(K562_type1_list)))), aes(x=Count)) +
  geom_histogram(binwidth=1, fill="steelblue") +
  scale_x_continuous(breaks=1:18) +
  labs(x="Number of datasets containing gene",y="Number of genes",title="K562 type1 symbol co-occurence") +
  theme_minimal()

## pheatmap
show_co_occurence(A549_type1_list,"A549 type1 Gene Presence/Absence")
show_co_occurence(K562_type1_list,"K562 type1 Gene Presence/Absence")


## re-filter
rna_A549_list <- lapply(rna_A549_list, function(df) {
  df <- df %>%
    dplyr::mutate(
      abundence_adjusted = dplyr::case_when(
        SYMBOL %in% high_freq_type2_genes_K562 ~ abundence,  # 条件1：基因在高频列表中
        Conflict_Type == "Type2" ~ Prob * abundence,                 # 条件2：Conflict_Type为Type2且基因不在高频列表中
        TRUE ~ abundence                                     # 默认情况：保留原始abundence
      )
    )
  # 统计NA值的数量
  na_count <- sum(is.na(df$abundence_adjusted))
  cat("当前数据框中 abundence_adjusted 列的 NA 值数量为:", na_count, "\n")
  return(df)
})

rna_K562_list <- lapply(rna_K562_list, function(df) {
  df <- df %>%
    dplyr::mutate(
      abundence_adjusted = dplyr::case_when(
        SYMBOL %in% high_freq_type2_genes_K562 ~ abundence,  # 条件1：基因在高频列表中
        Conflict_Type == "Type2" ~ Prob * abundence,                 # 条件2：Conflict_Type为Type2且基因不在高频列表中
        TRUE ~ abundence                                     # 默认情况：保留原始abundence
      )
    )
  # 统计NA值的数量
  na_count <- sum(is.na(df$abundence_adjusted))
  cat("当前数据框中 abundence_adjusted 列的 NA 值数量为:", na_count, "\n")
  return(df)
})



## for example
library(ggVennDiagram)
# 生成交互式韦恩图
P=ggVennDiagram(K562_type2_list[9:12]) +
  scale_fill_gradient(low = "white", high = "red") +  # 颜色映射
  labs(title = "K562 Type2 Gene Overlap") +                # 添加标题
  theme(plot.title = element_text(hjust = 0.5))       # 居中标题
ggsave("figures//K562_type2_gene_overlap.png", P, width = 6, height = 4, dpi = 400)

############################### flux test

## raw
library(METAFlux)
flux_A549_raw=data.frame(reaction=flux_leukemia$Reaction_ID)
for (i in 1:length(rna_A549_list)){
  tmp=as.data.frame(rna_A549_list[[i]]$abundence);
  rownames(tmp)=rna_A549_list[[i]]$SYMBOL;
  scores<-calculate_reaction_score(tmp)
  flux_A549_raw[,i+1]=compute_flux(mras=scores,medium=human_blood)
  gc()
}
rownames(flux_A549_raw)<-flux_A549_raw$reaction
flux_A549_raw=subset(flux_A549_raw,select = -c(reaction))

flux_K562_raw=data.frame(reaction=flux_leukemia$Reaction_ID)
for (i in 1:length(rna_K562_list)){
  tmp=as.data.frame(rna_K562_list[[i]]$abundence);
  rownames(tmp)=rna_K562_list[[i]]$SYMBOL;
  scores<-calculate_reaction_score(tmp)
  flux_K562_raw[,i+1]=compute_flux(mras=scores,medium=human_blood)
  gc()
}
rownames(flux_K562_raw)=flux_K562_raw$reaction
flux_K562_raw=subset(flux_K562_raw,select = -c(reaction))

score_A549_raw=pathway_comp(human_gem,flux_A549_raw,"A549_raw")
score_K562_raw=pathway_comp(human_gem,flux_K562_raw,"K562_raw")
score_A549_raw$mean=apply(score_A549_raw,1,mean);score_A549_raw$sd=apply(score_A549_raw,1,sd)
score_A549_raw$cv=score_A549_raw$sd/score_A549_raw$mean
score_K562_raw$mean=apply(score_K562_raw,1,mean);score_K562_raw$sd=apply(score_K562_raw,1,sd)
score_K562_raw$cv=score_K562_raw$sd/score_K562_raw$mean
score_A549_raw_cv=data.frame(pathway=rownames(score_A549_raw),cv=score_A549_raw$cv)
score_K562_raw_cv=data.frame(pathway=rownames(score_K562_raw),cv=score_K562_raw$cv)



## z-post
library(METAFlux)
flux_A549_comb=data.frame(reaction=flux_leukemia$Reaction_ID)
for (i in 1:length(rna_A549_list)){
  tmp=as.data.frame(rna_A549_list[[i]]$z_post);
  rownames(tmp)=rna_A549_list[[i]]$SYMBOL;
  scores<-calculate_reaction_score(tmp)
  flux_A549_comb[,i+1]=compute_flux(mras=scores,medium=human_blood)
  gc()
}
rownames(flux_A549_comb)<-flux_A549_comb$reaction
flux_A549_comb=subset(flux_A549_comb,select = -c(reaction))


flux_K562_comb=data.frame(reaction=flux_leukemia$Reaction_ID)
for (i in 1:length(rna_K562_list)){
  tmp=as.data.frame(rna_K562_list[[i]]$z_post);
  rownames(tmp)=rna_K562_list[[i]]$SYMBOL;
  scores<-calculate_reaction_score(tmp)
  flux_K562_comb[,i+1]=compute_flux(mras=scores,medium=human_blood)
  gc()
}
rownames(flux_K562_comb)=flux_K562_comb$reaction
flux_K562_comb=subset(flux_K562_comb,select = -c(reaction))

pathway_show(human_gem,flux_A549_comb,"A549_Bayesian")
pathway_show(human_gem,flux_K562_comb,"K562_Bayesian")




## filtered
library(METAFlux)
flux_A549_adjusted=data.frame(reaction=flux_leukemia$Reaction_ID)
for (i in 1:length(rna_A549_list)){
  temp=na.omit(rna_A549_list[[i]])
  tmp=as.data.frame(temp$abundence_adjusted);
  rownames(tmp)=temp$SYMBOL;
  scores<-calculate_reaction_score(tmp)
  flux_A549_adjusted[,i+1]=compute_flux(mras=scores,medium=human_blood)
  gc()
}
rownames(flux_A549_adjusted)<-flux_A549_adjusted$reaction
flux_A549_adjusted=subset(flux_A549_adjusted,select = -c(reaction))


flux_K562_adjusted=data.frame(reaction=flux_leukemia$Reaction_ID)
for (i in 1:length(rna_K562_list)){
  temp=na.omit(rna_K562_list[[i]])
  tmp=as.data.frame(temp$abundence_adjusted);
  rownames(tmp)=temp$SYMBOL;
  scores<-calculate_reaction_score(tmp)
  flux_K562_adjusted[,i+1]=compute_flux(mras=scores,medium=human_blood)
  gc()
}
rownames(flux_K562_adjusted)=flux_K562_adjusted$reaction
flux_K562_adjusted=subset(flux_K562_adjusted,select = -c(reaction))

pathway_show(human_gem,flux_A549_adjusted,"A549_Bayesian")
pathway_show(human_gem,flux_K562_adjusted,"K562_Bayesian")


############################### cv test

# 25个最具代表性的代谢通路
selected_pathways <- c(
  # 能量代谢
  "Glycolysis / Gluconeogenesis",
  "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",
  "Oxidative phosphorylation",
  "Pentose phosphate pathway",
  "Fatty acid oxidation",
  
  # 生物合成
  "Purine metabolism",
  "Pyrimidine metabolism",
  "Fatty acid biosynthesis (even-chain)",
  "Cholesterol biosynthesis 1 (Bloch pathway)",
  "Glutathione metabolism",
  
  # 氧化应激
  "ROS detoxification",
  "Ascorbate and aldarate metabolism",
  "Heme degradation",
  
  # 信号传导
  "Phosphatidylinositol phosphate metabolism",
  "Arachidonic acid metabolism",
  "Eicosanoid metabolism",
  "Prostaglandin biosynthesis",
  
  # 治疗靶点
  "Nicotinate and nicotinamide metabolism",  # NAD+代谢
  "Xenobiotics metabolism",                 # 化疗耐药
  "Drug metabolism",                        # 药物响应
  "Folate metabolism",                      # 抗叶酸治疗
  "Sphingolipid metabolism",                # 凋亡调控
  "Ubiquinone synthesis",                   # 线粒体靶向
  "Vitamin E metabolism",                   # 抗氧化治疗
  "Heme synthesis"                          # 铁死亡相关
)


## filter
score_A549_adjusted=pathway_comp(human_gem,flux_A549_adjusted,"A549_adjusted")
score_K562_adjusted=pathway_comp(human_gem,flux_K562_adjusted,"K562_adjusted")
## compute the pathway activity variation for each row
score_A549_adjusted$mean=apply(score_A549_adjusted,1,mean);score_A549_adjusted$sd=apply(score_A549_adjusted,1,sd)
score_A549_adjusted$cv=score_A549_adjusted$sd/score_A549_adjusted$mean
score_K562_adjusted$mean=apply(score_K562_adjusted,1,mean);score_K562_adjusted$sd=apply(score_K562_adjusted,1,sd)
score_K562_adjusted$cv=score_K562_adjusted$sd/score_K562_adjusted$mean
score_A549_adjusted_cv=data.frame(pathway=rownames(score_A549_adjusted),cv=score_A549_adjusted$cv)
score_K562_adjusted_cv=data.frame(pathway=rownames(score_K562_adjusted),cv=score_K562_adjusted$cv)

draw_cv_sorted(score_A549_adjusted_cv,"A549_adjusted",score_K562_adjusted_cv,"K562_adjusted",3)
draw_cv_sorted(score_A549_adjusted_cv,"A549_adjusted",score_A549_raw_cv,"A549 raw",3)
draw_cv_sorted(score_K562_adjusted_cv,"K562_adjusted",score_K562_raw_cv,"K562 raw",3)

draw_cv_sorted_select(score_A549_raw_cv,"A549",score_K562_raw_cv,"K562",3,selected_pathways)
draw_cv_sorted_select(score_A549_adjusted_cv,"A549 adjusted",score_K562_adjusted_cv,"K562 adjusted",3,selected_pathways)
draw_cv_sorted_select(score_A549_adjusted_cv,"A549_adjusted",score_A549_raw_cv,"A549 raw",3,selected_pathways)
draw_cv_sorted_select(score_K562_adjusted_cv,"K562_adjusted",score_K562_raw_cv,"K562 raw",3,selected_pathways)




## z-post
score_A549_comb=pathway_comp(human_gem,flux_A549_comb,"A549_comb")
score_K562_comb=pathway_comp(human_gem,flux_K562_comb,"K562_comb")
## compute the pathway activity variation for each row
score_A549_comb$mean=apply(score_A549_comb,1,mean);score_A549_comb$sd=apply(score_A549_comb,1,sd)
score_A549_comb$cv=score_A549_comb$sd/score_A549_comb$mean

score_K562_comb$mean=apply(score_K562_comb,1,mean);score_K562_comb$sd=apply(score_K562_comb,1,sd)
score_K562_comb$cv=score_K562_comb$sd/score_K562_comb$mean

score_A549_comb_cv=data.frame(pathway=rownames(score_A549_comb),cv=score_A549_comb$cv)
score_K562_comb_cv=data.frame(pathway=rownames(score_K562_comb),cv=score_K562_comb$cv)

draw_cv_sorted(score_A549_raw_cv,"A549 raw",score_K562_raw_cv,"K562 raw",4)
draw_cv_sorted(score_A549_comb_cv,"A549_comb",score_K562_comb_cv,"K562_comb",4)
draw_cv_sorted(score_A549_comb_cv,"A549_comb",score_A549_raw_cv,"A549 raw",4)
draw_cv_sorted(score_K562_comb_cv,"K562_comb",score_K562_raw_cv,"K562 raw",4)

draw_cv_sorted_select(score_A549_raw_cv,"A549",score_K562_raw_cv,"K562",top_number = 3,selected_pathways)
draw_cv_sorted_select(score_A549_comb_cv,"A549 comb",score_K562_comb_cv,"K562 comb",3,selected_pathways)
draw_cv_sorted_select(score_A549_comb_cv,"A549_comb",score_A549_raw_cv,"A549 raw",3,selected_pathways)
draw_cv_sorted_select(score_K562_comb_cv,"K562_comb",score_K562_raw_cv,"K562 raw",3,selected_pathways)




## draw
cv_significance_list=list(A549=data.frame(raw=score_A549_raw[selected_pathways,]$cv,comb=score_A549_comb[selected_pathways,]$cv,adjusted=score_A549_adjusted[selected_pathways,]$cv),
                          K562=data.frame(raw=score_K562_raw[selected_pathways,]$cv,comb=score_K562_comb[selected_pathways,]$cv,adjusted=score_K562_adjusted[selected_pathways,]$cv))

combined_data <- map2_dfr(cv_significance_list, names(cv_significance_list), ~ process_cell_line(.x, .y))
stat_results <- combined_data %>%
  group_by(cell_line) %>%
  group_modify(~ {
    pairwise_wilcox_test(
      data = .,
      formula = value ~ treatment,
      paired = TRUE,
      p.adjust.method = "BH"
    ) %>%
      add_xy_position(x = "treatment") %>%
      filter(!(group1 == "comb" & group2 == "adjusted")) # 示例筛选条件
  })

stat_results_filtered <- stat_results %>%
  filter(p.adj.signif != "ns") %>%  # 过滤非显著结果
  arrange(cell_line, xmin) %>%      # 按位置排序
  group_by(cell_line) %>%
  mutate(
    y.position = max(combined_data$value) * (1 + 0.1 * row_number())  # 动态调整高度
  )

# 2. 创建基础图形
p <- ggplot(combined_data, aes(x = treatment, y = value, fill = treatment)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.6) +
  geom_line(aes(group = subject), color = "gray50", alpha = 0.4) +
  scale_fill_manual(values = treatment_colors) +
  facet_wrap(~ cell_line, ncol = 2, scales = "free_y") +
  labs(title = "CV Comparisons",
       x = "Treatment",
       y = "CV Value") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

# 3. 手动添加显著性标记（按细胞系分步添加）
for(cl in unique(stat_results_filtered$cell_line)){
  sub_stats <- stat_results_filtered %>% 
    filter(cell_line == cl) %>%
    mutate(y.position = max(combined_data[combined_data$cell_line == cl, "value"]) * (1.05 + 0.1*row_number()))
  if(nrow(sub_stats) > 0){
    p <- p + 
      ggpubr::stat_pvalue_manual(
        data = sub_stats,
        label = "p.adj.signif",
        tip.length = 0.01,
        bracket.size = 0.6,
        xmin = "xmin",
        xmax = "xmax",
        y.position = "y.position",
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }
}

p <- p + 
  coord_cartesian(ylim = c(
    min(combined_data$value), 
    max(combined_data$value) * 1.2  # 留出20%空间
  ))

print(p)

# 保存图像
ggsave("figures//CV_Comparisons.png", p, width = 12, height = 8, dpi = 300)





p=draw_cv_sorted(score_A549_raw_cv,"A549",score_K562_raw_cv,"K562",3)
print(p)
# Save the plot
ggsave(
  filename = "figures//cv_sorted.jpg", # Recommended: PDF for vector graphics
  plot = p,
  width = 7,  # Adjust width based on journal requirements (e.g., single/double column)
  height = 5, # Adjust height as needed
  units = "in", # Units for width and height: "in", "cm", "mm", or "px"
  dpi = 600   # Dots per inch: 300-600 is typical for print. Higher for finer details.
)


# ggsave("figures//CV_sorted.png", p, width = 10, height = 8, dpi = 300)
p=draw_cv_sorted_select(score_A549_raw_cv,"A549",score_K562_raw_cv,"K562",5,selected_pathways)
print(p)
#ggsave("figures//CV_sorted_selected.png", p, width = 10, height = 8, dpi = 300)
ggsave(
  filename = "figures//cv_sorted_selected.jpg", # Recommended: PDF for vector graphics
  plot = p,
  width = 7,  # Adjust width based on journal requirements (e.g., single/double column)
  height = 5, # Adjust height as needed
  units = "in", # Units for width and height: "in", "cm", "mm", or "px"
  dpi = 600   # Dots per inch: 300-600 is typical for print. Higher for finer details.
)
