setwd("C:/Users/20145/Desktop/westlake_summer/");getwd()
library(openxlsx)
library(tidyverse)
complex_df=read.xlsx("./data/complex_df_with_p_values.xlsx")
distance_df <- read.csv("./data/distance_df.csv", header = TRUE, stringsAsFactors = FALSE)
distance_weight_df <- distance_df[,paste0(c("eta","epsilon","xc","scaling","xc_epsilon"),"_rela")]
rownames(distance_weight_df)=distance_df$gene


## p-value calculating functions
generate_null_dist_para <- function(n,N=10000,weight_df,para){
  para_col=paste0(para,"_rela")
  reduced_df=weight_df[,para_col,drop = FALSE]
  gene_list=as.character(rownames(reduced_df))
  result=matrix(NA,nrow = N,ncol=ncol(reduced_df))
  for (i in 1:N){
    gene_null=sample(gene_list, n, replace = FALSE)
    temp=reduced_df[gene_null,]
    scores=sum(temp)
    result[i,]=scores
  }
  colnames(result)=colnames(reduced_df)
  result=as.data.frame(result)
  return(result)
}


cal_p_SOFT_para <- function(genes,para,distance_weight_df){
  df=distance_weight_df
  para_name=paste0(para,"_rela")
  n=length(genes)
  temp=na.omit(distance_weight_df[genes,para_name])
  scores=sum(temp)
  N_sample=20000
  scores_null=generate_null_dist_para(n,N_sample,distance_weight_df,para)
  counts <- sum(scores_null < scores)
  p_value <- (counts+1) / (N_sample+1)
  return(p_value)
}


find_supercomplex <- function(gene_seed,para,network_matrix,distance_weight_df,max_iter=25,edge_cutoff=0.1){
  para_name=paste0(para,"_rela")
  # edge_cutoff =0.1
  genes_all=rownames(distance_weight_df)
  genes_mat=rownames(network_matrix)
  supercomplex <- tibble(
    iter = 0:max_iter,
    genes = vector("list", max_iter+1),
    p_value = rep(NA_real_, max_iter+1),
    size = rep(NA_real_, max_iter+1)
  )
  supercomplex$genes[1]=list(gene_seed)
  supercomplex$p_value[1] <- cal_p_SOFT_para(supercomplex$genes[[1]], para, distance_weight_df)
  supercomplex$size[1]=length(gene_seed)
  
  candidates <- c()
  for (i in 1:max_iter){
    print(paste0("Iteration ",i-1,": current p-value = ",supercomplex$p_value[i],"; size = ",supercomplex$size[i]))
    # find new neighbors
    outter_nodes=c()
    if (i==1){
      outter_nodes=supercomplex$genes[[i]]
      # print(paste0("Start with ", length(outter_nodes)," seeds"))
    }
    else{
      outter_nodes=setdiff(supercomplex$genes[[i]],supercomplex$genes[[i-1]])
    }
    # find new neighbors
    new_neighbors=c()
    for (j in 1:length(outter_nodes)){
      seed=outter_nodes[j]
      if (seed %in% genes_mat){
        temp=network_matrix[seed,]
        temp=na.omit(names(temp[temp>edge_cutoff]))
        new_neighbors=na.omit(c(new_neighbors, temp))
        new_neighbors=intersect(new_neighbors, genes_all)
      }
    }
    new_neighbors=unique(new_neighbors)
    # print(paste0("Found ",length(new_neighbors)," new neighbors"))
    # update candidates
    candidates <- setdiff(unique(c(candidates, new_neighbors)), supercomplex$genes[[i]])
    ## choose whether to add new nodes to the complex
    temp=distance_weight_df[candidates, ]
    temp$gene=rownames(temp)
    candidates=temp[order(temp[[para_name]]), "gene"]
    best_gene=na.omit(candidates[1])
    temp=na.omit(union(supercomplex$genes[[i]], best_gene))
    p_candidate <- cal_p_SOFT_para(temp, para, distance_weight_df)
    if (supercomplex$p_value[i]-p_candidate > 0 & length(best_gene)!=0) {
      # print(paste0("Adding new gene: ", best_gene))
      # add into supercomplex
      supercomplex$genes[[i+1]] <- union(supercomplex$genes[[i]], best_gene)
      supercomplex$p_value[i+1] <- p_candidate
      supercomplex$size[i+1] <- length(supercomplex$genes[[i+1]])
      candidates <- setdiff(candidates, best_gene)
    }
    else{
      # print("Nothing could be improved! Quiting...")
      return(na.omit(supercomplex))
    }
  }
  return(supercomplex)
}
temp=find_supercomplex("HAP4","xc",genetic_all_mat,distance_weight_df,10,0.2)


## plot
library(ggplot2)
library(ggforce)
library(ggnewscale)
library(ggrepel)
library(patchwork)

plot_supercomplex <- function(old_genes,new_genes,para){
p <- ggplot(data.frame(x = seq(0, 2, length.out = 1000)),aes(x=x))+
  geom_hline(yintercept=1, color = "grey",linewidth = 0.3)+
  geom_vline(xintercept=1, color = "grey",linewidth = 0.3)+
  geom_line(aes(y = eta_x(x)),size = 0.3) +
  geom_line(aes(y = epsilon_x(x)),size = 0.3) +
  geom_line(aes(y = xc_x(x)),size = 0.3) +
  geom_point(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,shape=para),
             size=0.5,alpha=0.5)+
  scale_shape_manual(values=c(1,2,3), name="Parameters")+
  theme(
    panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(size = 10),     # 调整主标题大小
    plot.subtitle = element_text(size = 8),  # 调整副标题大小
    axis.title = element_text(size = 6),     # 调整所有坐标轴标题大小
    axis.text = element_text(size = 5)       # 调整所有坐标轴标签大小
  )
p1 <- p +
  geom_mark_hull(data=subset(lifes_set_alpha_bp_short_df, gene %in% old_genes),
                 aes(x=mean_rela, y=steepness_rela), size=0.1,
                 alpha=0.2, expand=0.05,show.legend = FALSE,fill="pink",color="red") +
  geom_point(data=subset(lifes_set_alpha_bp_short_df, gene %in% old_genes),
             aes(x=mean_rela, y=steepness_rela), alpha=1,color="red",size =0.75 )+
  geom_text_repel(data=subset(lifes_set_alpha_bp_short_df, gene %in% old_genes),
                  aes(x=mean_rela, y=steepness_rela, label=gene),color="red",
                  size=1.5, nudge_x = 0.05, nudge_y = 0.05,
                  segment.size = 0.2,max.overlaps = 20, box.padding = 0.1,
                  show.legend = FALSE) +
  xlim(0,2)+ylim(0,2)+
  coord_fixed(ratio=1)+
  theme(legend.position = "none")+
  labs(x="Lifespan Relative to WT", y="Steepness Relative to WT",
       title=paste0("OLD complex p-value: ",round(cal_p_SOFT_para(old_genes,para,distance_weight_df),4)," (",para,")"))

p2 <- p + 
  geom_mark_hull(data=subset(lifes_set_alpha_bp_short_df, gene %in% new_genes),
                 aes(x=mean_rela, y=steepness_rela), size=0.1,
                 alpha=0.2, expand=0.05,show.legend = FALSE,fill="pink",color="red") +
  geom_point(data=subset(lifes_set_alpha_bp_short_df, gene %in% new_genes),
             aes(x=mean_rela, y=steepness_rela), alpha=1,color="red",size=0.75)+
  geom_text_repel(data=subset(lifes_set_alpha_bp_short_df, gene %in% new_genes),
                  aes(x=mean_rela, y=steepness_rela, label=gene),color="red",
                  size=1.5, nudge_x = 0.05, nudge_y = 0.05,
                  segment.size = 0.2,max.overlaps = 20, box.padding = 0.1,
                  show.legend = FALSE) +
  xlim(0,2)+ylim(0,2)+
  coord_fixed(ratio=1)+
  # theme(legend.position = "none")+
  labs(x="Lifespan Relative to WT", y="Steepness Relative to WT",
       title=paste0("SUPER complex p-value: ",round(cal_p_SOFT_para(new_genes,para,distance_weight_df),4)," (",para,")"),
       )
p3=p1 / p2
return(p3)
}


# find_supercomplex("VMA16", "epsilon", genetic_all_mat, distance_weight_df)

### for each genes & complexes as seeds
eta_length=length(eta_genes)
epsilon_length=length(epsilon_genes)
xc_length=length(xc_genes)
xc_epsilon_length=length(xc_epsilon_genes)
scaling_length=length(scaling_genes)
all_length = eta_length + epsilon_length + xc_length + xc_epsilon_length + scaling_length

## single genes
supercomplex_single_genes_df = tibble(
  gene=c(eta_genes, epsilon_genes, xc_genes, xc_epsilon_genes, scaling_genes),
  hard_cluster=c(rep("eta", eta_length),
         rep("epsilon", epsilon_length),
         rep("xc", xc_length),
         rep("xc_epsilon", xc_epsilon_length),
         rep("scaling", scaling_length)),
  supercomplex=NA_real_,
  p_value=NA_real_,
  size=NA_real_,
  plot=NA_real_
)

# start growing every seeds
for (i in 1:all_length){
  gene_seed=supercomplex_single_genes_df$gene[i]
  para=supercomplex_single_genes_df$hard_cluster[i]
  print(paste0("Growing supercomplex for ",gene_seed," with parameter ",para))
  supercomplex <- find_supercomplex(gene_seed, para, genetic_all_mat, distance_weight_df)
  supercomplex_single_genes_df$supercomplex[i] <- list(supercomplex$genes[[length(supercomplex$genes)]])
  supercomplex_single_genes_df$p_value[i] <- supercomplex$p_value[length(supercomplex$p_value)]
  supercomplex_single_genes_df$size[i] <- supercomplex$size[length(supercomplex$size)]
  supercomplex_single_genes_df$plot[i] <- list(plot_supercomplex(
    old_genes = supercomplex_single_genes_df$gene[i],
    new_genes = supercomplex_single_genes_df$supercomplex[[i]],
    para = para
  ))
}

## complexes
supercomplex_complex_df=tibble(
  Complex_ac=c(complex_list$eta$Complex_ac,complex_list$epsilon$Complex_ac,
               complex_list$xc$Complex_ac,complex_list$xc_epsilon$Complex_ac,
               complex_list$scaling$Complex_ac),
  cluster=c(rep("eta", nrow(complex_list$eta)),
         rep("epsilon", nrow(complex_list$epsilon)),
         rep("xc", nrow(complex_list$xc)),
         rep("xc_epsilon", nrow(complex_list$xc_epsilon)),
         rep("scaling", nrow(complex_list$scaling))),
  gene_seed=c(complex_list$eta$genes,complex_list$epsilon$genes,
              complex_list$xc$genes,complex_list$xc_epsilon$genes,
              complex_list$scaling$genes),
  old_p_value=c(complex_list$eta$SOFT_eta_p,complex_list$epsilon$SOFT_epsilon_p,
                complex_list$xc$SOFT_xc_p,complex_list$xc_epsilon$SOFT_xc_epsilon_p,
                complex_list$scaling$SOFT_scaling_p),
  old_size=NA_real_,
  supercomplex=NA_real_,
  p_value=NA_real_,
  size=NA_real_,
  plot=NA_real_
)

# start growing seeds
for (i in 1:nrow(supercomplex_complex_df)) {
  complex=supercomplex_complex_df$Complex_ac[i]
  para=supercomplex_complex_df$cluster[i]
  seeds=unlist(supercomplex_complex_df$gene_seed[i])
  supercomplex_complex_df$old_size[i] <- length(seeds)
  print(paste0("Growing supercomplex for ",complex," with parameter ",para))
  supercomplex <- find_supercomplex(seeds, para, genetic_all_mat, distance_weight_df)
  supercomplex_complex_df$supercomplex[i] <- list(supercomplex$genes[[length(supercomplex$genes)]])
  supercomplex_complex_df$p_value[i] <- supercomplex$p_value[length(supercomplex$p_value)]
  supercomplex_complex_df$size[i] <- supercomplex$size[length(supercomplex$size)]
  supercomplex_complex_df$plot[i] <- list(plot_supercomplex(
    old_genes = seeds,
    new_genes = supercomplex_complex_df$supercomplex[[i]],
    para = para
  ))
}

# save results
# 
# # remove plot col
# supercomplex_single_genes_save_df=subset(supercomplex_single_genes_df, select = -c(plot))
# supercomplex_complex_save_df=subset(supercomplex_complex_df, select = -c(plot))
# # save into .xlsx
# write.xlsx(supercomplex_single_genes_save_df, "./data/supercomplex_single_genes_save_df.xlsx")
# write.xlsx(supercomplex_complex_save_df, "./data/supercomplex_complex_save_df.xlsx")
