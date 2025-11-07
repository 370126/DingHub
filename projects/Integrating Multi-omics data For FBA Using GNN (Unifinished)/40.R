## raw
test_A549_raw=flux_A549_raw[,];test_A549_raw$Reaction_ID=rownames(test_A549_raw)
test_K562_raw=flux_K562_raw[,];test_K562_raw$Reaction_ID=rownames(test_K562_raw)

test_A549_raw=merge(test_A549_raw,experi_A549,by="Reaction_ID")
rownames(test_A549_raw)=test_A549_raw$Reaction_ID;test_A549_raw=subset(test_A549_raw,select=-c(Reaction_ID))
test_K562_raw=merge(test_K562_raw,experi_K562,by="Reaction_ID")
rownames(test_K562_raw)=test_K562_raw$Reaction_ID;test_K562_raw=subset(test_K562_raw,select=-c(Reaction_ID))

cor_A549_raw=rep(0,length(colnames(test_A549_raw))-1)
for (i in 1:length(cor_A549_raw)) {
  cor_A549_raw[i]=cor(as.double(test_A549_raw[,i]),as.double(test_A549_raw$experimental_flux),method="spearman")
}
cor_K562_raw=rep(0,length(colnames(test_K562_raw))-1)
for (i in 1:length(cor_K562_raw)){
  cor_K562_raw[i]=cor(as.double(test_K562_raw[,i]),as.double(test_K562_raw$experimental_flux),method="spearman")
}

## comb
test_A549_comb=flux_A549_comb[,];test_A549_comb$Reaction_ID=rownames(test_A549_comb)
test_K562_comb=flux_K562_comb[,];test_K562_comb$Reaction_ID=rownames(test_K562_comb)

test_A549_comb=merge(test_A549_comb,experi_A549,by="Reaction_ID")
rownames(test_A549_comb)=test_A549_comb$Reaction_ID;test_A549_comb=subset(test_A549_comb,select=-c(Reaction_ID))
test_K562_comb=merge(test_K562_comb,experi_K562,by="Reaction_ID")
rownames(test_K562_comb)=test_K562_comb$Reaction_ID;test_K562_comb=subset(test_K562_comb,select=-c(Reaction_ID))

cor_A549_comb=rep(0,length(colnames(test_A549_comb))-1)
for (i in 1:length(cor_A549_comb)) {
  cor_A549_comb[i]=cor(as.double(test_A549_comb[,i]),as.double(test_A549_comb$experimental_flux),method="spearman")
}
cor_K562_comb=rep(0,length(colnames(test_K562_comb))-1)
for (i in 1:length(cor_K562_comb)){
  cor_K562_comb[i]=cor(as.double(test_K562_comb[,i]),as.double(test_K562_comb$experimental_flux),method="spearman")
}

## adjusted
test_A549_adjusted=flux_A549_adjusted[,];test_A549_adjusted$Reaction_ID=rownames(test_A549_adjusted)
test_K562_adjusted=flux_K562_adjusted[,];test_K562_adjusted$Reaction_ID=rownames(test_K562_adjusted)

test_A549_adjusted=merge(test_A549_adjusted,experi_A549,by="Reaction_ID")
rownames(test_A549_adjusted)=test_A549_adjusted$Reaction_ID;test_A549_adjusted=subset(test_A549_adjusted,select=-c(Reaction_ID))
test_K562_adjusted=merge(test_K562_adjusted,experi_K562,by="Reaction_ID")
rownames(test_K562_adjusted)=test_K562_adjusted$Reaction_ID;test_K562_adjusted=subset(test_K562_adjusted,select=-c(Reaction_ID))

cor_A549_adjusted=rep(0,length(colnames(test_A549_adjusted))-1)
for (i in 1:length(cor_A549_adjusted)) {
  cor_A549_adjusted[i]=cor(as.double(test_A549_adjusted[,i]),as.double(test_A549_adjusted$experimental_flux),method="spearman")
}
cor_K562_adjusted=rep(0,length(colnames(test_K562_adjusted))-1)
for (i in 1:length(cor_K562_adjusted)){
  cor_K562_adjusted[i]=cor(as.double(test_K562_adjusted[,i]),as.double(test_K562_adjusted$experimental_flux),method="spearman")
}













## barplot correlation
library(ggplot2)
library(tidyr)
cor_raw_df=data.frame(A549_raw=cor_A549_raw,K562_raw=cor_K562_raw,sample=paste("sample",1:length(cor_A549_raw),sep=""));cor_raw_df=pivot_longer(cor_raw_df,cols = 1:2,names_to = "cell line",values_to = "correlation")
ggplot(cor_raw_df,aes(x=sample,y=correlation,fill=`cell line`))+
  geom_bar(stat="identity",position="dodge")+
  coord_cartesian(ylim = c(0.3, 0.8)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title="correlation between predicted flux and experimental flux",x="sample",y="correlation")+
  ## remove background
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "grey"))
ggsave("figures//correlation_raw.png",p,width = 12, height = 8, dpi = 300)



cor_comb_df=data.frame(A549_comb=cor_A549_comb,K562_comb=cor_K562_comb,sample=paste("sample",1:length(cor_A549_comb),sep=""));cor_comb_df=pivot_longer(cor_comb_df,cols = 1:2,names_to = "cell line",values_to = "correlation")
ggplot(cor_comb_df,aes(x=sample,y=correlation,fill=`cell line`))+
  geom_bar(stat="identity",position="dodge")+
  coord_cartesian(ylim = c(0, 0.8)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title="correlation between predicted flux and experimental flux (comb)",x="reaction",y="correlation")

cor_adjusted_df=data.frame(A549_adjusted=cor_A549_adjusted,K562_adjusted=cor_K562_adjusted,sample=paste("sample",1:length(cor_A549_adjusted),sep=""));cor_adjusted_df=pivot_longer(cor_adjusted_df,cols = 1:2,names_to = "cell line",values_to = "correlation")
ggplot(cor_adjusted_df,aes(x=sample,y=correlation,fill=`cell line`))+
  geom_bar(stat="identity",position="dodge")+
  coord_cartesian(ylim = c(0, 0.8)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title="correlation between predicted flux and experimental flux (adjusted)",x="reaction",y="correlation")

cor_significance_list=list(A549=data.frame(raw=cor_A549_raw,comb=cor_A549_comb,adjusted=cor_A549_adjusted),
                           K562=data.frame(raw=cor_K562_raw,comb=cor_K562_comb,adjusted=cor_K562_adjusted))
cor_significance_list








###################### visualization ######################



library(tidyverse)
library(rstatix)
library(ggpubr)
# 定义处理函数
process_cell_line <- function(data, cell_name) {
  # 转换为长格式并添加样本ID
  data_long <- data %>%
    mutate(subject = 1:n()) %>%
    pivot_longer(cols = -subject,
                 names_to = "treatment",
                 values_to = "value") %>%
    mutate(treatment = factor(treatment),
           cell_line = cell_name)
  
  return(data_long)
}

# 处理两个细胞系数据
combined_data <- map2_dfr(cor_significance_list, names(cor_significance_list), ~ process_cell_line(.x, .y))

# 1. 计算统计结果（过滤非显著值）
stat_results_filtered <- combined_data %>%
  group_by(cell_line) %>%
  group_modify(~ {
    # Friedman检验（保留整体结果）
    friedman_res <- friedman_test(data = ., formula = value ~ treatment | subject)
    
    # 事后检验并过滤非显著结果
    pairwise_wilcox_test(
      data = .,
      formula = value ~ treatment,
      paired = TRUE,
      p.adjust.method = "BH"
    ) %>%
      add_xy_position(x = "treatment") %>%
      filter(p.adj.signif != "ns") %>%  # 关键过滤步骤
      mutate(
        friedman_statistic = friedman_res$statistic,
        friedman_p = friedman_res$p
      )
  })

# 2. 可视化（自动适配显著性标记）
p <- ggplot(combined_data, aes(x = treatment, y = value, fill = treatment)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.6) +
  geom_line(aes(group = subject), color = "gray50", alpha = 0.4) +
  scale_fill_manual(values = treatment_colors) +
  stat_pvalue_manual(
    stat_results_filtered,
    label = "p.adj.signif",
    tip.length = 0.01,
    bracket.size = 0.6,
    step.increase = 0.1,
    facet.by = "cell_line",  # 确保分面匹配
    inherit.aes = FALSE
  ) +
  facet_wrap(~ cell_line, ncol = 2, scales = "free_y") +
  labs(title = "Correlation Comparisons",
       x = "data type",
       y = "correlation value (Spearman)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

print(p)
# 保存图像
ggsave("figures//correlation_comparison.png", p, width = 12, height = 8, dpi = 300)
