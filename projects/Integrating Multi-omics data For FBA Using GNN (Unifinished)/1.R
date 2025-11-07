#install.packages("devtools")
#devtools::install_github('KChen-lab/METAFlux')

setwd("C:/Users/86189/Desktop/neo_project")
library(METAFlux)
library(ggplot2)
library(tidyr)
library(ggpubr)

Hgem <- METAFlux:::Hgem
obj=Hgem$Obj
# 找到 Hgem$Obj 中非零项的索引
#non_zero_indices <- which(Hgem$Obj != 0)
# 提取对应的反应名称
#non_zero_reactions <- Hgem$Reaction[non_zero_indices]
# 输出非零项对应的反应名称
#print(Hgem$Reaction[non_zero_indices])
#print(Hgem$pathway[non_zero_indices])


# 读取数据
ALL_LINE= read.table("OmicsExpressionProteinCodingGenesTPMLogp1.csv", header = T, sep = ",")
ALL=ALL_LINE
colnames(ALL_LINE) <- gsub("\\.\\.[0-9]*\\.", "", colnames(ALL_LINE))

meta=read.csv("Model.csv", header = T, sep = ",")
ID2NAME=meta[,c(1,4)];colnames(ID2NAME) <- c("X", "Cellline")

NCI2CT=read.table("NCI60.csv", header = F,sep=',');
NCI2CT[,1]=gsub("-| |\\(.*\\)|\\/.*","",NCI2CT[,1])
colnames(NCI2CT) <- c("cell_line", "cancer_type");NCI2CT$cell_line=toupper(NCI2CT$cell_line)

# ALL_LINE1=ALL_LINE #备份
rownames(ALL_LINE)=ALL_LINE$X;ALL_LINE <- ALL_LINE[, -1]
rownames(ALL_LINE) <- ID2NAME$Cellline[match(rownames(ALL_LINE), ID2NAME$X)]
lines=as.data.frame(rownames(ALL_LINE))

## NCI_60
NCI_60=c(
  "A549/ATCC", "BT-549", "COLO 205", "HCC-2998", "HCT-116", "HCT-15", "HT29", "KM12", 
  "SW-620", "SF-268", "SF-295", "SF-539", "SNB-19", "SNB-75", "U251", "LOX IMVI", 
  "MALME-3M", "MCF7", "MDA-MB-231", "HS 578T", "MDA-MB-435", "SK-MEL-2", 
  "SK-MEL-28", "SK-MEL-5", "UACC-257", "UACC-62", "IGROV1", "OVCAR-3", "OVCAR-4", 
  "OVCAR-5", "OVCAR-8", "NCI/ADR-RES", "SK-OV-3", "786-0", "A498", "ACHN", "CAKI-1", 
  "RXF 393", "SN12C", "TK-10", "UO-31", "PC-3", "DU-145", "MCF10A", "MCF10F", 
  "MDA-MB-468", "T-47D", "HL-60(TB)", "K-562", "MOLT-4", "RPMI-8226", "SR", 
  "CCRF-CEM", "KARPAS-299", "SUP-T1", "NCI-H23", "NCI-H226", "NCI-H322M", "NCI-H460", 
  "NCI-H522", "HOP-62", "HOP-92"
)
NCI_60=gsub("-| |\\(.*\\)|\\/.*","",NCI_60);NCI_60
ALL_LINE_NCI60 <- ALL_LINE[rownames(ALL_LINE) %in% NCI_60,]

## 正常的细胞系做对照
normal=meta[meta$OncotreePrimaryDisease=="Non-Cancerous",]$StrippedCellLineName
ALL_LINE_Normal <- ALL_LINE[rownames(ALL_LINE) %in% normal,]
set.seed(666)  # 设置随机种子，保证结果可重复
ALL_LINE_Normal <- ALL_LINE_Normal[sample(nrow(ALL_LINE_Normal), 30), ]


data("human_gem")
data("human_blood")
## 表观遗传相关反应
SAM_Reactions=c("HMR_4241","HMR_8025","HMR_8026","HMR_8027","HMR_7765",'HMR_3875')
CoA_Reactions=c("HMR_1090","HMR_4916")

length1=length(ALL_LINE_NCI60[,1])
flux_NCI60=data.frame(demethylation=rep(NA, length1), methylation1=rep(NA, length1),
                      methylation2=rep(NA, length1),methylation3=rep(NA, length1),
                      SAM_transport=rep(NA, length1),aceCOA_transport=rep(NA, length1),
                      SAM_synthesis=rep(NA, length1),COA_transport=rep(NA, length1),
                      biomass_flux=rep(NA, length1),type=rep("NCI60", length1))
for (i in 1:length(ALL_LINE_NCI60[,1])){
  Cell_line <- rownames(ALL_LINE_NCI60)[i]
  cell_line_data <- ALL_LINE_NCI60[i,]
  cell_line_data <- t(cell_line_data)
  scores<-calculate_reaction_score(cell_line_data)
  flux<-compute_flux(mras=scores,medium=human_blood) 
  biomass_flux=t(obj) %*% flux
  flux_NCI60$cell_line[i] <- Cell_line
  flux_NCI60$biomass_flux[i] <- biomass_flux
  flux_NCI60$demethylation[i] <- flux["HMR_4241",]
  flux_NCI60$methylation1[i] <- flux["HMR_8025",]
  flux_NCI60$methylation2[i] <- flux["HMR_8026",]
  flux_NCI60$methylation3[i] <- flux["HMR_8027",]
  flux_NCI60$SAM_transport[i] <- flux["HMR_7765",]
  flux_NCI60$SAM_synthesis[i] <- flux["HMR_3875",]
  flux_NCI60$COA_transport[i] <- flux["HMR_4916",]
  flux_NCI60$aceCOA_transport[i] <- flux["HMR_1090",]
}

length2=length(ALL_LINE_Normal[,1])
flux_normal=data.frame(demethylation=rep(NA, length2), methylation1=rep(NA, length2),
                      methylation2=rep(NA, length2),methylation3=rep(NA, length2),
                      SAM_transport=rep(NA, length2),SAM_synthesis=rep(NA, length2),
                      aceCOA_transport=rep(NA, length2),COA_transport=rep(NA, length2),
                      biomass_flux=rep(NA, length2),type=rep("Normal", length2))
for (i in 1:length(ALL_LINE_Normal[,1])){
  Cell_line <- rownames(ALL_LINE_Normal)[i]
  cell_line_data <- ALL_LINE_Normal[i,]
  cell_line_data <- t(cell_line_data)
  scores<-calculate_reaction_score(cell_line_data)
  flux<-compute_flux(mras=scores,medium=human_blood) 
  biomass_flux=t(obj) %*% flux
  flux_normal$cell_line[i] <- Cell_line
  flux_normal$biomass_flux[i] <- biomass_flux
  flux_normal$demethylation[i] <- flux["HMR_4241",]
  flux_normal$methylation1[i] <- flux["HMR_8025",]
  flux_normal$methylation2[i] <- flux["HMR_8026",]
  flux_normal$methylation3[i] <- flux["HMR_8027",]
  flux_normal$SAM_transport[i] <- flux["HMR_7765",]
  flux_normal$SAM_synthesis[i] <- flux["HMR_3875",]
  flux_normal$COA_transport[i] <- flux["HMR_4916",]
  flux_normal$aceCOA_transport[i] <- flux["HMR_1090",]
}






#####
total= rbind(flux_NCI60, flux_normal);total=na.omit(total);total_ori=total
total_norm <- as.data.frame(lapply(total[,1:9], function(x) (x - mean(x)) / sd(x)));total_norm=cbind(total_norm,total[,c(10,11)])
total=total_norm;rm(total_norm)
total$methylation_sum=total$methylation1+total$methylation2+total$methylation3-total$demethylation+total$SAM_transport+total$SAM_synthesis

#save(total, file = "total.RData")


total_merge=merge(total,NCI2CT,by="cell_line",all.x=TRUE);
total=total_merge;rm(total_merge)
## 各种癌症
#breast<-na.omit(breast);breast=rbind(breast,total[total$type=="Normal",])
lung=total[total$cancer_type=="Non-Small Cell Lung",];lung=na.omit(lung);lung=rbind(lung,(total[total$type=="Normal",])[1:length(lung$type)+2,])
CNS=total[total$cancer_type=="CNS",];CNS=na.omit(CNS);CNS=rbind(CNS,(total[total$type=="Normal",])[1:length(CNS$type)+2,])
melanoma=total[total$cancer_type=="Melanoma",];melanoma=na.omit(melanoma);melanoma=rbind(melanoma,total[total$type=="Normal",])
colon=total[total$cancer_type=="Colon",];colon=na.omit(colon);colon=rbind(colon,total[total$type=="Normal",][1:length(colon$type)+2,])
ovarian=total[total$cancer_type=="Ovarian",];ovarian=na.omit(ovarian); ovarian=rbind(ovarian,total[total$type=="Normal",][1:length(ovarian$type)+2,])
renal=total[total$cancer_type=="Renal",];renal=na.omit(renal);renal=rbind(renal,total[total$type=="Normal",][1:length(renal$type)+2,])
prostate=total[total$cancer_type=="Prostate",];prostate=na.omit(prostate);prostate=rbind(prostate,total[total$type=="Normal",][1:length(prostate$type)+2,])
leukemia=total[total$cancer_type=="Leukemia",];leukemia=na.omit(leukemia);leukemia=rbind(leukemia,total[total$type=="Normal",][1:length(leukemia$type)+2,])

# 总的
total_long=pivot_longer(total, cols = c(demethylation, methylation1, methylation2, methylation3, SAM_transport, SAM_synthesis, aceCOA_transport, COA_transport,methylation_sum), names_to = "reaction", values_to = "flux")
ggplot(data=total_long, aes(x=reaction, y=flux, fill=type)) + 
  geom_boxplot(outlier.shape = NA) + ylim(-3,+3)+ 
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = "Flux of NCI60 and Normal Cell Lines", x = "Cell Line Type", y = "Flux") + 
  scale_fill_manual(values = c("NCI60" = "pink", "Normal" = "skyblue"))+
  stat_compare_means( method = "t.test" ,label = "p.signif",na.rm = TRUE,hide.ns = TRUE)
ggplot(data=total_long,aes(x=type,y=flux,fill=type)) + geom_boxplot(outlier.shape = NA) +
  facet_wrap(~reaction) + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = "Flux of NCI60 and Normal Cell Lines", x = "Cell Line Type", y = "Flux") + 
  scale_fill_manual(values = c("NCI60" = "pink", "Normal" = "skyblue"))+ylim(-3,3)+stat_compare_means(method = "t.test", label = "p.signif",hide.ns = TRUE)
ggplot(data=total,aes(x=type,y=biomass_flux,fill=type)) + geom_boxplot() + theme_minimal() + labs(title = "Biomass Flux of NCI60 and Normal Cell Lines", x = "Cell Line Type", y = "Flux") + scale_fill_manual(values = c("NCI60" = "pink", "Normal" = "skyblue"))+stat_compare_means(method = "t.test", label = "p.signif")



ggplot(data=total,aes(x=type,y=methylation1,fill=type)) + geom_boxplot() + theme_minimal() + 
  labs(title = "Methylation1 Flux of NCI60 and Normal Cell Lines", x = "Cell Line Type", y = "Flux") + 
  scale_fill_manual(values = c("NCI60" = "pink", "Normal" = "skyblue"))+
  stat_compare_means(method = "t.test", label = "p.signif")
ggplot(data=total,aes(x=type,y=methylation2,fill=type)) + geom_boxplot() + theme_minimal() + 
  labs(title = "Methylation2 Flux of NCI60 and Normal Cell Lines", x = "Cell Line Type", y = "Flux") + 
  scale_fill_manual(values = c("NCI60" = "pink", "Normal" = "skyblue"))+
  stat_compare_means(method = "t.test", label = "p.signif")
ggplot(data=total,aes(x=type,y=methylation3,fill=type)) + geom_boxplot() + theme_minimal() + 
  labs(title = "Methylation3 Flux of NCI60 and Normal Cell Lines", x = "Cell Line Type", y = "Flux") + 
  scale_fill_manual(values = c("NCI60" = "pink", "Normal" = "skyblue"))+
  stat_compare_means(method = "t.test", label = "p.signif")
ggplot(data=total,aes(x=type,y=SAM_transport,fill=type)) + geom_boxplot() + theme_minimal() +
  labs(title = "SAM Transport Flux of NCI60 and Normal Cell Lines", x = "Cell Line Type", y = "Flux") + 
  scale_fill_manual(values = c("NCI60" = "pink", "Normal" = "skyblue"))+
  stat_compare_means(method = "t.test", label = "p.signif")
ggplot(data=total,aes(x=type,y=SAM_synthesis,fill=type)) + geom_boxplot() + theme_minimal() +
  labs(title = "SAM Synthesis Flux of NCI60 and Normal Cell Lines", x = "Cell Line Type", y = "Flux") + 
  scale_fill_manual(values = c("NCI60" = "pink", "Normal" = "skyblue"))+
  stat_compare_means(method = "t.test", label = "p.signif")
ggplot(data=total,aes(x=type,y=methylation_sum,fill=type)) + geom_boxplot() + theme_minimal() +
  labs(title = "Methylation Flux of NCI60 and Normal Cell Lines", x = "Cell Line Type", y = "Flux") + 
  scale_fill_manual(values = c("NCI60" = "pink", "Normal" = "skyblue"))+
  stat_compare_means(method = "t.test", label = "p.signif")
ggplot(data=total,aes(x=type,y=biomass_flux,fill=type)) + geom_boxplot() + theme_minimal() +
  labs(title = "Biomass Flux of NCI60 and Normal Cell Lines", x = "Cell Line Type", y = "Flux") + 
  scale_fill_manual(values = c("NCI60" = "pink", "Normal" = "skyblue"))+
  stat_compare_means(method = "t.test", label = "p.signif")



total_long_melanoma=pivot_longer(melanoma, cols = c(demethylation, methylation1, methylation2, methylation3, SAM_transport, SAM_synthesis, aceCOA_transport, COA_transport,biomass_flux,methylation_sum), names_to = "reaction", values_to = "flux")
ggplot(data=total_long_melanoma, aes(x=reaction, y=flux, fill=type)) + geom_boxplot(outlier.shape = NA) + ylim(-3,3)+theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = "Flux of Melanoma and Normal Cell Lines", x = "Cell Line Type", y = "Flux") + scale_fill_manual(values = c("NCI60" = "pink", "Normal" = "skyblue"))+
  stat_compare_means(method = "t.test", label = "p.signif",group.by = "type",hide.ns = TRUE)
ggplot(data=melanoma,aes(x=type,y=biomass_flux,fill=type)) + geom_boxplot() + theme_minimal() + labs(title = "Biomass Flux of Melanoma and Normal Cell Lines", x = "Cell Line Type", y = "Flux") + scale_fill_manual(values = c("NCI60" = "pink", "Normal" = "skyblue"))+stat_compare_means(method = "t.test", label = "p.signif")

total_long_colon=pivot_longer(colon, cols = c(demethylation, methylation1, methylation2, methylation3, SAM_transport, SAM_synthesis, aceCOA_transport, COA_transport,biomass_flux), names_to = "reaction", values_to = "flux")
ggplot(data=total_long_colon, aes(x=reaction, y=flux, fill=type)) + geom_boxplot(outlier.shape = NA) + ylim(-3,3)+theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = "Flux of Colon and Normal Cell Lines", x = "Cell Line Type", y = "Flux") + scale_fill_manual(values = c("NCI60" = "pink", "Normal" = "skyblue"))+
  stat_compare_means(method = "t.test", label = "p.signif",group.by = "type",hide.ns = TRUE)
ggplot(data=colon,aes(x=type,y=biomass_flux,fill=type)) + geom_boxplot() + theme_minimal() + labs(title = "Biomass Flux of Colon and Normal Cell Lines", x = "Cell Line Type", y = "Flux") + scale_fill_manual(values = c("NCI60" = "pink", "Normal" = "skyblue"))+stat_compare_means(method = "t.test", label = "p.signif")

total_long_ovarian=pivot_longer(ovarian, cols = c(demethylation, methylation1, methylation2, methylation3, SAM_transport, SAM_synthesis, aceCOA_transport, COA_transport,biomass_flux), names_to = "reaction", values_to = "flux")
ggplot(data=total_long_ovarian, aes(x=reaction, y=flux, fill=type)) + geom_boxplot(outlier.shape = NA) + ylim(-3,3)+theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = "Flux of Ovarian and Normal Cell Lines", x = "Cell Line Type", y = "Flux") + scale_fill_manual(values = c("NCI60" = "pink", "Normal" = "skyblue"))+
  stat_compare_means(method = "t.test", label = "p.signif",group.by = "type",hide.ns = TRUE)
ggplot(data=ovarian,aes(x=type,y=biomass_flux,fill=type)) + geom_boxplot() + theme_minimal() + labs(title = "Biomass Flux of Ovarian and Normal Cell Lines", x = "Cell Line Type", y = "Flux") + scale_fill_manual(values = c("NCI60" = "pink", "Normal" = "skyblue"))+stat_compare_means(method = "t.test", label = "p.signif")

total_long_lung=pivot_longer(lung, cols = c(demethylation, methylation1, methylation2, methylation3, SAM_transport, SAM_synthesis, aceCOA_transport, COA_transport,biomass_flux), names_to = "reaction", values_to = "flux")
ggplot(data=total_long_lung, aes(x=reaction, y=flux, fill=type)) + geom_boxplot(outlier.shape = NA) + ylim(-3,3)+theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = "Flux of Lung and Normal Cell Lines", x = "Cell Line Type", y = "Flux") + scale_fill_manual(values = c("NCI60" = "pink", "Normal" = "skyblue"))+
  stat_compare_means(method = "t.test", label = "p.signif",group.by = "type",hide.ns = TRUE)
ggplot(data=lung,aes(x=type,y=biomass_flux,fill=type)) + geom_boxplot() + theme_minimal() + labs(title = "Biomass Flux of Non-small Cell Lung and Normal Cell Lines", x = "Cell Line Type", y = "Flux") + scale_fill_manual(values = c("NCI60" = "pink", "Normal" = "skyblue"))+stat_compare_means(method = "t.test", label = "p.signif")








