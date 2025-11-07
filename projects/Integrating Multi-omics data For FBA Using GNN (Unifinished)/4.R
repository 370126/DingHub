setwd("C:/Users/86189/Desktop/neo_project")
#load("workspace.RData")
gc()

############ primitive validation of the models
library(readxl)
#experi <- read_excel("aaz1482_data_file_s5.xlsx",skip = 1)
#colnames(experi)=gsub("-| |\\(.*\\)|\\/.*","",colnames(experi))

experi_more <- read_excel("1218595databases1_corrected(1).xls",sheet = "CORE data")
experi_more=experi_more[experi_more$`Method (a)`=="HILIC" | experi_more$`Method (a)`=="BGA",]

selected_metabolites <- c(
  "adenosine",             # 能量代谢(ATP)、信号传递
  "alanine",               # 糖异生、氮代谢
  "arginine",              # 尿素循环、NO合成、免疫调节
  "asparagine",            # 氨基酸代谢、肿瘤代谢
  "aspartate",             # TCA循环、嘧啶合成
  "glutamine",             # 氮载体、肿瘤代谢关键营养
  "glutamate",             # 神经递质、氨基酸代谢枢纽
  "glycine",               # 单碳代谢、血红素合成
  "choline",               # 磷脂代谢、甲基供体
  "phosphocholine",        # 膜磷脂代谢关键中间体
  "creatine",              # 能量缓冲系统(磷酸肌酸)
  "glucose",               # 核心能量代谢底物
  "isoleucine",            # 支链氨基酸(BCAA)、蛋白质合成
  "lactate",               # 乳酸循环、代谢酸碱平衡
  "leucine",               # BCAA、mTOR信号调控
  "lysine",                # 组蛋白修饰(乙酰化)、蛋白质合成
  "methionine",            # 甲基供体(SAM合成)
  "serine",                # 单碳代谢、核苷酸合成
  "threonine",             # 免疫球蛋白糖基化、代谢调控
  "tryptophan",            # NAD+合成、神经递质前体
  "tyrosine",              # 神经递质/激素合成
  "valine",                # BCAA、能量代谢
  "taurine",               # 抗氧化、渗透调节
  "spermidine",            # 自噬诱导、DNA稳定
  "citrulline",            # 尿素循环、精氨酸再生
  "betaine",               # 甲基化反应、渗透保护
  "homocysteine"           # 甲基循环、心血管疾病标志
)

experi_more=experi_more[experi_more$`Metabolite (b)` %in% selected_metabolites,]

original_colnames <- gsub("\\.\\.\\.\\d+$", "", colnames(experi_more[,2:length(colnames(experi_more))]))

exchange_reaction=human_gem[human_gem$SUBSYSTEM=="Exchange/demand reactions",c("ID","EQUATION")]
exchange_reaction$Metabolite=gsub("\\[s\\].*","",exchange_reaction$EQUATION)
exchange_reaction$Metabolite=gsub("D-lactate","lactate",exchange_reaction$Metabolite)

experi_reactions=data.frame(Metabolite=experi_more$`Metabolite (b)`);
experi_reactions=merge(experi_reactions,exchange_reaction,by="Metabolite")

experi=experi_more[experi_more$`Metabolite (b)` %in% exchange_reaction$Metabolite,]

library(dplyr)
library(purrr)
library(tibble)
# 取平均
metab_col <- experi %>% select(`Metabolite (b)`)
processed_data <- experi %>%
  select(-`Metabolite (b)`) %>%
  split.default(gsub("\\.{3}\\d+$", "", names(.))) %>%
  imap(~ {
    if (ncol(.x) > 1) {
      rowMeans(.x, na.rm = TRUE)
    } else {
      .x[[1]]
    }
  }) %>%
  as_tibble(.name_repair = "unique") %>% 
  rename_with(~ gsub("\\.{3}\\d+$", "", .))  # 安全替换后缀
experi_mean <- bind_cols(experi["Metabolite (b)"], processed_data)
colnames(experi_mean)=gsub("-| |\\(.*\\)|\\/.*","",colnames(experi_mean))

#experi_mean <- subset(experi_mean, select = -c(Calibrated, Method))
temp=experi_reactions[,c("Metabolite","ID")]
experi_integrated <- merge(experi_mean,temp,by="Metabolite")


################## weighted mean
experi_K562=cbind(experi_integrated[,c("ID")],experi_integrated[,"K562"])
colnames(experi_K562)=c("Reaction_ID","experimental_flux")
flux_leukemia$Reaction_ID=rownames(flux_leukemia)
test_leukemia=merge(flux_leukemia,experi_K562,by="Reaction_ID")
rownames(test_leukemia)=test_leukemia$Reaction_ID;test_leukemia=test_leukemia[,-1]

correlation=rep(0,length(colnames(test_leukemia)))
for (i in 1:length(colnames(test_leukemia))){
  correlation[i]=cor(as.double(test_leukemia[,i]),as.double(test_leukemia$experimental_flux),method="spearman")
}
test_leukemia_correlation=rbind(test_leukemia,correlation)

experi_A549=cbind(experi_integrated[,c("ID")],experi_integrated[,"A549"])
colnames(experi_A549)=c("Reaction_ID","experimental_flux")
flux_lung$Reaction_ID=rownames(flux_leukemia)
test_lung=merge(flux_lung,experi_A549,by="Reaction_ID")
rownames(test_lung)=test_lung$Reaction_ID;test_lung=test_lung[,-1]

correlation=rep(0,length(colnames(test_lung)))
for (i in 1:length(colnames(test_lung))){
  correlation[i]=cor(as.double(test_lung[,i]),as.double(test_lung$experimental_flux),method="spearman")
}

test_lung_correlation=rbind(test_lung,correlation)





## visualization
library(ggplot2)
library(tidyr)
#test_long_leukemia=pivot_longer(test_leukemia,cols=1:(length(colnames(test_leukemia))-1),names_to="Weight",values_to="Flux")
#gplot(test_long_leukemia, aes(x=Flux,y=experimental_flux,colour = Weight)) + geom_point() + theme_minimal()
#test_long_lung=pivot_longer(test_lung,cols=1:(length(colnames(test_lung))-1),names_to="Weight",values_to="Flux")
#ggplot(test_long_lung, aes(x=Flux,y=experimental_flux,colour = Weight)) + geom_point() + theme_minimal()

## Bar plot
test_weighted_total <- cbind(
  t(test_lung_correlation[nrow(test_lung_correlation), 1:(ncol(test_lung_correlation) - 1)]),
                          t(test_leukemia_correlation[nrow(test_leukemia_correlation), 1:(ncol(test_leukemia_correlation) - 1)])
  )
                            
                             
colnames(test_weighted_total)=c("Lung","Leukemia");
test_weighted_total=as.data.frame(test_weighted_total)
test_weighted_total$weight=as.numeric(rownames(test_weighted_total))
test_weighted_total_long=pivot_longer(as.data.frame(test_weighted_total),cols=1:2,names_to="Cell_line",values_to="Correlation")

ggplot(data=test_weighted_total_long,aes(x=weight,y=Correlation,fill=Cell_line)) + 
  geom_bar(stat="identity",position="dodge")+
  coord_cartesian(ylim = c(0.3, 0.8)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title="correlation between predicted flux and experimental flux",x="weight of transcriptom",y="correlation")+
  ## remove background
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "grey"))
  #geom_bar(stat="identity",position="dodge") + theme_minimal()+labs(title="Weighted mean")

################## Bayesian model
library(METAFlux)
bay_lung=as.data.frame(comb_lung$z_post_prob);rownames(bay_lung)=comb_lung$SYMBOL
scores<-calculate_reaction_score(bay_lung)
flux_bay_lung<-compute_flux(mras=scores,medium=human_blood)
flux_bay_lung=as.data.frame(flux_bay_lung)
gc()
bay_leukemia=as.data.frame(comb_leukemia$z_post_prob);rownames(bay_leukemia)=comb_leukemia$SYMBOL
scores<-calculate_reaction_score(bay_leukemia)
flux_bay_leukemia<-compute_flux(mras=scores,medium=human_blood)
flux_bay_leukemia=as.data.frame(flux_bay_leukemia)
gc()


flux_bay_lung$Reaction_ID=rownames(flux_leukemia)
test_bay_lung=merge(flux_bay_lung,experi_A549,by="Reaction_ID")
rownames(test_bay_lung)=test_bay_lung$Reaction_ID;test_bay_lung=test_bay_lung[,-1]
cor(as.double(test_bay_lung$`comb_lung$z_post_prob`),as.double(test_bay_lung$experimental_flux),method="spearman")

flux_bay_leukemia$Reaction_ID=rownames(flux_leukemia)
test_bay_leukemia=merge(flux_bay_leukemia,experi_K562,by="Reaction_ID")
rownames(test_bay_leukemia)=test_bay_leukemia$Reaction_ID;test_bay_leukemia=test_bay_leukemia[,-1]
cor(as.double(test_bay_leukemia$`comb_leukemia$z_post_prob`),as.double(test_bay_leukemia$experimental_flux),method="spearman")


################## Multiply
flux_multi_leukemia$Reaction_ID=rownames(flux_leukemia)
flux_multi_leukemia=subset(flux_multi_leukemia,select=-c(reaction))
test_multi_leukemia=merge(flux_multi_leukemia,experi_K562,by="Reaction_ID")
rownames(test_multi_leukemia)=test_multi_leukemia$Reaction_ID;test_multi_leukemia=test_multi_leukemia[,-1]
correlation=rep(0,length(colnames(test_multi_leukemia)))
for (i in 1:length(colnames(test_multi_leukemia))){
  correlation[i]=cor(as.double(test_multi_leukemia[,i]),as.double(test_multi_leukemia$experimental_flux),method="spearman")
}
test_multi_leukemia_correlation=rbind(test_multi_leukemia,correlation)

flux_multi_lung$Reaction_ID=rownames(flux_leukemia);
flux_multi_lung=subset(flux_multi_lung,select=-c(reaction))
test_multi_lung=merge(flux_multi_lung,experi_A549,by="Reaction_ID")
rownames(test_multi_lung)=test_multi_lung$Reaction_ID;test_multi_lung=test_multi_lung[,-1]
correlation=rep(0,length(colnames(test_multi_lung)))
for (i in 1:length(colnames(test_multi_lung))){
  correlation[i]=cor(as.double(test_multi_lung[,i]),as.double(test_multi_lung$experimental_flux),method="spearman")
}
test_multi_lung_correlation=rbind(test_multi_lung,correlation)

## Bar plot
test_multi_total <- cbind(
  t(test_multi_lung_correlation[nrow(test_multi_lung_correlation), 1:(ncol(test_multi_lung_correlation) - 1)]),
  t(test_multi_leukemia_correlation[nrow(test_multi_leukemia_correlation), 1:(ncol(test_multi_leukemia_correlation) - 1)])
)
colnames(test_multi_total)=c("Lung","Leukemia");
test_multi_total=as.data.frame(test_multi_total)
test_multi_total$weight=as.numeric(rownames(test_multi_total))
test_multi_total_long=pivot_longer(as.data.frame(test_multi_total),cols=1:2,names_to="Cell_line",values_to="Correlation")
ggplot(data=test_multi_total_long,aes(x=weight,y=Correlation,fill=Cell_line)) + geom_bar(stat="identity",position="dodge") + theme_minimal()+labs(title="Multiply")


############################################
