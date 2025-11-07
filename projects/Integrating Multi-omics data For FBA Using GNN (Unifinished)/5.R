setwd("C:/Users/86189/Desktop/neo_project")

leukemia_lines <- meta[
  (meta$OncotreeLineage == "Lymphoid" | meta$OncotreeLineage == "Myeloid") & 
    meta$OncotreePrimaryDisease == "Non-Cancerous",
  "StrippedCellLineName"
]
leukemia_rna = ALL_LINE[leukemia_lines,];
leukemia_rna=na.omit(leukemia_rna);leukemia_rna=t(leukemia_rna);leukemia_rna=as.data.frame(leukemia_rna)
leukemia_mutated=leukemia_rna;
mutated_genes <- c("IDH1", "IDH2","FLT3")
leukemia_mutated[rownames(leukemia_mutated) %in% mutated_genes, ] <- 0
colnames(leukemia_mutated) <- gsub("$","_IDH1/2_silenced",colnames(leukemia_mutated))
leukemia_together=cbind(leukemia_rna,leukemia_mutated)
leukemia_together=leukemia_together[,order(colnames(leukemia_together))]

flux_leukemia_mutated=data.frame(reaction=rownames(flux_leukemia))
library(METAFlux)
for (i in 1:ncol(leukemia_together)){
  scores<-calculate_reaction_score(leukemia_together[i]);
  temp=compute_flux(mras=scores,medium=human_blood);
  flux_leukemia_mutated[[as.character(colnames(leukemia_together)[i])]] = as.vector(temp)
  gc()
}
rownames(flux_leukemia_mutated)=flux_leukemia_mutated$reaction;flux_leukemia_mutated=flux_leukemia_mutated[,-1]
pathway_show(human_gem,flux_leukemia_mutated,"leukemia mutated")
pathway_show_mutated(human_gem,flux_leukemia_mutated,"leukemia mutated")

histon_reactions=c("HMR_4241","HMR_8025","HMR_8026","HMR_8027","HMR_7765",'HMR_3875',"HMR_1090","HMR_4916","HMR_4097","biomass_human")
histon_reactions_name=c("Histone demethylation","Histone methylation1","Histone methylation2","Histone methylation3","SAM_transport","SAM_synthesis","COA_transport","aceCOA_transport","aceCOA_synthesis","biomass_human")
histon_reactions=data.frame(reaction=histon_reactions,name=histon_reactions_name)
flux_histon_reactions=flux_leukemia_mutated[histon_reactions$reaction,];rownames(flux_histon_reactions)=histon_reactions$name
pheatmap::pheatmap(flux_histon_reactions,cluster_cols = F,color = rev(mapal),scale = "row",fontsize = 6,main = "Histone modification related Reactions")
i=1
diff_histon <- data.frame(reaction = rownames(flux_histon_reactions))
while (i <= ncol(flux_histon_reactions)) {
  temp=(flux_histon_reactions[,i+1]-flux_histon_reactions[,i])/abs(flux_histon_reactions[,i]);
  diff_histon[[colnames(flux_histon_reactions)[i]]] = temp;
  i=i+2;
}
rownames(diff_histon) <- diff_histon$reaction;diff_histon <- diff_histon[,-1]
pheatmap::pheatmap(diff_histon,cluster_cols = F,cluster_rows = T,color = rev(mapal),scale = "column",fontsize = 6,main = "Histone modification related Reactions change rate")





#####################
mutation_seq=list(raw="EMPTY_GENE",IDH1=c("IDH1"),IDH2=c("IDH2"),DNMT3A=c("DNMT3A"),TET2="TET2",
                  IDH12=c("IDH1","IDH2"),IDH12DNMT3A=c("IDH1","IDH2","DNMT3A"),IDH12DNMT3ATET2=c("IDH1","IDH2","DNMT3A","TET2"),IDH12TET2=c("IDH1","IDH2","TET2"),DNMT3ATET2=c("DNMT3A","TET2"),
                  add_NPM1=c("IDH1","IDH2","DNMT3A","TET2","NPM1"),add_ASXL1=c("IDH1","IDH2","DNMT3A","TET2","ASXL1"),add_TP53=c("IDH1","IDH2","DNMT3A","TET2","TP53"),all_mutation=c("IDH1","IDH2","DNMT3A","TET2","NPM1","ASXL1","TP53"))
leukemia_trajectory_list=list(HCC1739BL=data.frame(Reaction_ID=rownames(flux_leukemia)),
                                                   IM95=data.frame(Reaction_ID=rownames(flux_leukemia)))
library(METAFlux)
for (i in 1:length(leukemia_trajectory_list)) {
  sample_name <- names(leukemia_trajectory_list)[i]
  rna <- as.data.frame(leukemia_rna[, sample_name, drop = FALSE])
  rownames(rna) <- rownames(leukemia_rna)
  for (j in 1:length(mutation_seq)) {
    mutation_name <- names(mutation_seq)[j]
    target_genes <- mutation_seq[[j]]
    rna_mutated <- rna
    # 检查基因是否存在并赋值
    existing_genes <- intersect(rownames(rna_mutated), target_genes)
    if (length(existing_genes) > 0) {
      rna_mutated[existing_genes, ] <- 0
    } else {
      warning(paste("No genes found for mutation:", mutation_name))
    }
    print(rna_mutated[target_genes,])
    scores <- calculate_reaction_score(as.matrix(rna_mutated))
    temp <- compute_flux(mras = scores, medium = human_blood)
    if (length(temp) == nrow(leukemia_trajectory_list[[i]])) {
      leukemia_trajectory_list[[i]][[mutation_name]] <- temp
    } else {
      stop(paste("Dimension mismatch for mutation:", mutation_name))
    }
  }
}


for (i in 1:length(leukemia_trajectory_list)){
  rownames(leukemia_trajectory_list[[i]])=rownames(flux_leukemia)
  leukemia_trajectory_list[[i]]=subset(leukemia_trajectory_list[[i]],select=-c(Reaction_ID))
  colnames(leukemia_trajectory_list[[i]])=gsub("\\[.*","",colnames(leukemia_trajectory_list[[i]]))
}





######################### PCA #########################
tmp=leukemia_trajectory_list[[1]]
tmp=cbind(tmp,flux_K562_comb[1:18])

filter_traj <- function(tmp,num_features=10000){
tmp=t(tmp)
variances <- apply(tmp, 2, var)
keep_columns <- order(variances, decreasing = TRUE)[1:num_features]
tmp_filtered <- tmp[, keep_columns]
tmp<- scale(tmp_filtered, center = TRUE, scale = TRUE) 
return(tmp)
}

pca_result <- prcomp(tmp)
library(factoextra)
# 绘制前两个主成分的样本分布
fviz_pca_ind(pca_result, 
             col.ind = "cos2",   # 颜色表示样本对主成分的余弦相似度
             gradient.cols = c("blue",  "red"),
             alpha.var = 0.5,
             repel = TRUE)


tmp=leukemia_trajectory_list[[1]]
tmp=tmp[histon_reactions$reaction,]
tmp=cbind(tmp,flux_A549_comb[histon_reactions$reaction,])
tmp=t(tmp)
tmp<- scale(tmp, center = TRUE, scale = TRUE) 
pca_result <- prcomp(tmp)
library(factoextra)
# 绘制前两个主成分的样本分布
fviz_pca_ind(pca_result, 
             col.ind = "cos2",   # 颜色表示样本对主成分的余弦相似度
             gradient.cols = c("blue",  "red"),
             alpha.var = 0.5,
             repel = TRUE)



######################### UMAP #########################

## data preparation
# all fluxes
traj_HCC1739BL_all=leukemia_trajectory_list[[1]];traj_IM95_all=leukemia_trajectory_list[[2]]
traj_HCC1739BL_all=cbind(traj_HCC1739BL_all,flux_K562_comb[1:18]);traj_IM95_all=cbind(traj_IM95_all,flux_K562_comb[1:18])
# histone modification related fluxes
traj_HCC1739BL_histon=traj_HCC1739BL_all[histon_reactions$reaction,];traj_IM95_histon=traj_IM95_all[histon_reactions$reaction,]
traj_HCC1739BL_histon=cbind(traj_HCC1739BL_histon,flux_A549_comb[histon_reactions$reaction,]);traj_IM95_histon=cbind(traj_IM95_histon,flux_A549_comb[histon_reactions$reaction,])
# p=pheatmap::pheatmap(traj_HCC1739BL_histon,cluster_cols = T,border_color = NA,treeheight_row = 0,treeheight_col = 20,
#                      cellwidth = 7,cellheight = 5,
#                    color = rev(mapal),scale = "row",fontsize = 6,main = "Histone modification related Reactions: HCC1739BL")
# ggsave("figures//HCC1739BL_Histon.png",p,width=6,height=3.5,dpi=400)
# p=pheatmap::pheatmap(traj_IM95_histon,cluster_cols = T,border_color = NA,treeheight_row = 0,treeheight_col = 20,
#                      cellwidth = 7,cellheight = 5,
#                    color = rev(mapal),scale = "row",fontsize = 6,main = "Histone modification related Reactions: IM95")
# ggsave("figures//IM95_Histon.png",p,width=6,height=3.5,dpi=400)

# random sample
set.seed(888)
random_reactions=sample(human_gem$ID,10)
traj_HCC1739BL_random=traj_HCC1739BL_all[random_reactions,]
traj_HCC1739BL_random=cbind(traj_HCC1739BL_random,flux_K562_comb[random_reactions,])

# filter traj
traj_HCC1739BL_all=filter_traj(traj_HCC1739BL_all)
traj_IM95_all=filter_traj(traj_IM95_all)
traj_HCC1739BL_histon=filter_traj(traj_HCC1739BL_histon,length(histon_reactions$reaction))
traj_IM95_histon=filter_traj(traj_IM95_histon,length(histon_reactions$reaction))
traj_HCC1739BL_random=filter_traj(traj_HCC1739BL_random,length(random_reactions))
traj_HCC1739BL_random=traj_HCC1739BL_random[1:32,]
## UMAP
graph_HCC1739BL_all=umap_cluter(traj_HCC1739BL_all,tissue_name="HCC1739BL_all")
graph_IM95_all=umap_cluter(traj_IM95_all,tissue_name="IM95_all")
graph_HCC1739BL_histon=umap_cluter(traj_HCC1739BL_histon,tissue_name="HCC1739BL_histon")
graph_IM95_histon=umap_cluter(traj_IM95_histon,tissue_name="IM95_histon")
graph_HCC1739BL_random=umap_cluter(traj_HCC1739BL_random,tissue_name="HCC1739BL_random")

graph_HCC1739BL_random[[2]]

ggsave(
  filename = "figures//HCC1739BL_all_umap.jpg", # Recommended: PDF for vector graphics
  plot = graph_HCC1739BL_all[[2]],
  width = 9,  # Adjust width based on journal requirements (e.g., single column is often 3.5 inches, double 7 inches)
  height = 6, # Adjust height as needed to prevent label overlap or squishing
  units = "in", # Units for width and height: "in", "cm", "mm", or "px"
  dpi = 600   # Dots per inch: 300-600 is typical for print. Higher for finer details.
)
# # save high quality figure
# ggsave("figures//HCC1739BL_all_umap.png",graph_HCC1739BL_all[[2]],width=8,height=6,dpi=300)
# ggsave("figures//IM95_all_umap.png",graph_IM95_all[[2]],width=8,height=6,dpi=300)
# ggsave("figures//HCC1739BL_histon_umap.png",graph_HCC1739BL_histon[[2]],width=8,height=6,dpi=300)
# ggsave("figures//IM95_histon_umap.png",graph_IM95_histon[[2]],width=8,height=6,dpi=300)
# ggsave("figures//HCC1739BL_random_umap.png",graph_HCC1739BL_random[[2]],width=8,height=6,dpi=300)


## compare with raw data
# all fluxes
traj_HCC1739BL_all_raw=leukemia_trajectory_list[[1]];traj_IM95_all_raw=leukemia_trajectory_list[[2]]
traj_HCC1739BL_all_raw=cbind(traj_HCC1739BL_all_raw,flux_K562_raw[1:18]);traj_IM95_all_raw=cbind(traj_IM95_all_raw,flux_K562_raw[1:18])

# filter traj
traj_HCC1739BL_all_raw=filter_traj(traj_HCC1739BL_all_raw)
traj_IM95_all_raw=filter_traj(traj_IM95_all_raw)

## UMAP
graph_HCC1739BL_all_raw=umap_cluter(traj_HCC1739BL_all_raw,tissue_name="HCC1739BL_all_raw")

graph_IM95_all_raw=umap_cluter(traj_IM95_all_raw,tissue_name="IM95_all_raw")
graph_HCC1739BL_all_raw[[2]]


ggsave(
  filename = "figures//HCC1739BL_all_raw_umap.jpg", # Recommended: PDF for vector graphics
  plot = graph_HCC1739BL_all_raw[[2]],
  width = 9,  # Adjust width based on journal requirements (e.g., single column is often 3.5 inches, double 7 inches)
  height = 6, # Adjust height as needed to prevent label overlap or squishing
  units = "in", # Units for width and height: "in", "cm", "mm", or "px"
  dpi = 600   # Dots per inch: 300-600 is typical for print. Higher for finer details.
)
ggsave("figures//HCC1739BL_all_raw_umap.png",graph_HCC1739BL_all_raw[[2]],width=8,height=6,dpi=300)
ggsave("figures//IM95_all_raw_umap.png",graph_IM95_all_raw[[2]],width=8,height=6,dpi=300)



colnames(flux_K562_raw)=paste("K562_sample",c(1:length(flux_K562_raw)),sep="_")
colnames(flux_A549_raw)=paste("A549_sample",c(1:length(flux_A549_raw)),sep="_")
p <- pathway_show(human_gem,flux_K562_raw,"K562 Pathway Activity")
# save large figure
ggsave("figures//K562_Pathway_Activity.png",p,width=7,height=12,dpi=400)
p <- pathway_show(human_gem,flux_A549_raw,"A549 Pathway Activity")
ggsave("figures//A549_Pathway_Activity.png",p,width=7,height=12,dpi=400)









