## load functions
setwd("C:/Users/20145/Desktop/westlake_summer/");getwd()
source("./SR_fitting.R")
## load data
setwd("C:/Users/20145/Desktop/westlake_summer/database/genetic_interaction");getwd()
library(data.table)
library(readxl)
library(tidyverse)
# distance_df$hard_cluster <- distance_short_df$hard_cluster
# write.csv(distance_df, "../../data/distance_df.csv", row.names = FALSE)
# write.csv(lifes_set_alpha_bp_df,"../../data/lifes_set_alpha_bp_df.csv",row.names = FALSE)
# write.csv(counts_alpha,"../../data/counts_alpha.csv",row.names = FALSE)
# write into .xlsx
library(openxlsx)
# write.xlsx(lifes_set_alpha_complex_df, "../../data/lifes_set_alpha_complex.xlsx", row.names = FALSE)

# NEW distance df
distance_df <- read.csv("./data/distance_df_new.csv", header = TRUE, stringsAsFactors = FALSE)

counts_alpha=read.csv("./data/counts_alpha.csv",header=TRUE, stringsAsFactors = FALSE)
lifes_set_alpha_bp_df <- fread("./data/lifes_set_alpha_bp_df.csv", sep="," ,header=TRUE, stringsAsFactors = FALSE)
lifes_set_alpha_bp_df=as.data.frame(lifes_set_alpha_bp_df);head(lifes_set_alpha_bp_df)

simu_long=read.csv("./data/simu_long.csv")

lifes_set_alpha_bp_df$hard_cluster <- distance_df$hard_cluster[match(lifes_set_alpha_bp_df$gene, distance_df$gene)]
distance_df$num_of_experi <- counts_alpha$num_of_experi[match(distance_df$gene, counts_alpha$gene)]

temp=lifes_set_alpha_bp_df %>%
  group_by(gene) %>%
  # calculate the standard deviation
  summarize(
    lifespan_sd = sd(mean_rela_bootstrap, na.rm = TRUE),
    steepness_sd = sd(steepness_rela_bootstrap, na.rm = TRUE),
  )
# match with distance_df
distance_df$lifespan_sd <- temp$lifespan_sd[match(distance_df$gene, temp$gene)]
distance_df$steepness_sd <- temp$steepness_sd[match(distance_df$gene, temp$gene)]
rm(temp)
gc()

### find closest genes & sub-matrices
threshold <- 0.1
min_length <- 0.32
min_sd <- 0.05

# genes itself
eta_genes = subset(distance_df, eta_rela <= threshold & length >= min_length & lifespan_sd <= min_sd & steepness_sd <= min_sd)$gene
epsilon_genes = subset(distance_df, epsilon_rela <= threshold & length >= min_length & lifespan_sd <= min_sd & steepness_sd <= min_sd)$gene
xc_genes = subset(distance_df, xc_rela <= threshold & length >= min_length & lifespan_sd <= min_sd & steepness_sd <= min_sd)$gene
xc_epsilon_genes=subset(distance_df, xc_epsilon_rela <= threshold & length >= min_length & lifespan_sd <= min_sd & steepness_sd <= min_sd)$gene
scaling_genes = subset(distance_df, scaling_rela <= threshold & length >= min_length & lifespan_sd <= min_sd & steepness_sd <= min_sd)$gene
close_genes= unique(c(eta_genes, xc_epsilon_genes,scaling_genes))

# close_genes=top20_genes
# close_genes=top10_genes
# close_genes="SOD2"
# plot on the graph
library(ggplot2)
library(ggforce)
library(ggnewscale)
library(ggrepel)
ggplot(data.frame(x = seq(0, 2, length.out = 1000)),aes(x=x))+
  geom_hline(yintercept=1, color = "grey",linewidth = 0.3)+
  geom_vline(xintercept=1, color = "grey",linewidth = 0.3)+
  geom_line(aes(y = eta_x(x)),size = 0.3) +
  geom_line(aes(y = epsilon_x(x)),size = 0.3) +
  geom_line(aes(y = xc_x(x)),size = 0.3) +
  geom_point(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,shape=para),
             size=0.75,alpha=0.5)+
  scale_shape_manual(values = c(1, 2, 3), name = "simu para") +
  geom_point(data=subset(distance_df, gene %in% close_genes),
             aes(x=mean_rela, y=steepness_rela, 
                 color=hard_cluster), alpha=1,size=1.2)+
  geom_text_repel(data=subset(distance_df, gene %in% close_genes),
                 aes(x=mean_rela, y=steepness_rela, label=gene,color=hard_cluster),
                 size=1.5, nudge_x = 0.05, nudge_y = 0.05,
                 segment.size = 0.2,max.overlaps = 20, box.padding = 0.1,
                 show.legend = FALSE) +
  scale_color_brewer(palette = "Set1", name = "genes closest to") +
  xlim(0,1.75)+ylim(0,1.75)+
  #theme(legend.position = "none")+
  labs(x="Lifespan Relative to WT", y="Steepness Relative to WT",
       title=paste("Relative distance threshold=",threshold,", Max SD=",min_sd)) +
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


# confidence interval
label <- function(df,x,y){
  temp=df %>%
    group_by(gene) %>%
    summarize(
      # choose 95% confidence interval as the label position
      x = quantile(!!sym(x), probs = 0.90, na.rm = TRUE),
      y = quantile(!!sym(y), probs = 0.90, na.rm = TRUE),
    )
  temp$hard_cluster=df$hard_cluster[match(temp$gene, df$gene)]
  return(temp)
}

# select genes with xc_rela_new top10
close_genes=distance_df$gene[order(distance_df$eta_rela_new)][1:10]

label_df= label(subset(lifes_set_alpha_bp_df, gene %in% close_genes), "mean_rela_bootstrap", "steepness_rela_bootstrap")
ggplot(subset(lifes_set_alpha_bp_df,gene %in% close_genes), 
       aes(x = mean_rela_bootstrap, y = steepness_rela_bootstrap)) +
  geom_hline(yintercept=1, color = "grey",linewidth = 0.3)+
  geom_vline(xintercept=1, color = "grey",linewidth = 0.3)+
  geom_line(data=data.frame(x=seq(0, 3, length.out = 1000), y=eta_x(seq(0, 3, length.out = 1000))),
            aes(x=x, y=y), size=0.3 ) +
  geom_line(data=data.frame(x=seq(0, 3, length.out = 1000), y=epsilon_x(seq(0, 3, length.out = 1000))),
            aes(x=x, y=y), size=0.3) +
  geom_line(data=data.frame(x=seq(0, 3, length.out = 1000), y=xc_x(seq(0, 3, length.out = 1000))),
            aes(x=x, y=y), size=0.3) +
  geom_point(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,shape=para),
            size=0.75,alpha=0.5,color="black")+
  scale_shape_manual(values = c(1, 2, 3)) +
  new_scale_color() +
  geom_point(size = 0.01, alpha = 0.02, color="gray" ,show.legend = FALSE) +
  stat_ellipse(aes(color=hard_cluster,group=gene), level= 0.95, size = 0.5) +
  scale_color_brewer(palette = "Set1", name = "genes closest to") +
  geom_text_repel(data=label_df, aes(x = x, y = y, label = gene, color=hard_cluster),
                  size = 2, nudge_x = 0.05, nudge_y = 0.05,
                  segment.size = 0.1, max.overlaps = 20, box.padding = 0.1,
                  show.legend = FALSE) +
  xlim(0, 2) + ylim(0, 2) +
  coord_fixed(ratio=1) +
  labs(x="Lifespan Relative to WT", y="Steepness Relative to WT",
       title = "With 95% CI" ) +
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







### gene interaction data
# gene interaction score
# dt <- fread("./Data File S1. Raw genetic interaction datasets_ Pair-wise interaction format/SGA_NxN.txt",header=TRUE);dt=as.data.frame(dt)
# dt=subset(dt,`P-value`<0.05)
# dt=dt[,c(2,4,6,7)]
# to upper all string element in firsts 2 columns
# dt <- dt %>%
#   mutate(across(1:2, toupper))
# genetic_interaction_df=dt
# gc()
# length(unique(genetic_interaction_df$`Query allele name`))
# length(unique(genetic_interaction_df$`Array allele name`))
# intersect(close_genes,c(genetic_interaction_df$`Query allele name`,genetic_interaction_df$`Array allele name`))


# gene profile similarity
dt <- fread("./Data File S3. Genetic interaction profile similarity matrices/Global_Similarity_Network.txt", header = FALSE);dt=as.data.frame(dt)
# first 2 rows and first 2 columns are the names of the genes
mat <- as.matrix(dt[-c(1, 2), -c(1, 2)]);mat=apply(mat, 2, as.numeric)
rownames(mat) <- toupper(dt[-c(1, 2), 1]);colnames(mat) <- toupper(dt[1, -c(1, 2)])
diag(mat) <- NA
genetic_all_mat <- as.matrix(mat)
print(nrow(genetic_all_mat))
# dt <- fread("./Data File S3. Genetic interaction profile similarity matrices/Essential_Similarity_Network.txt", header = FALSE);dt=as.data.frame(dt)
# mat <- as.matrix(dt[-c(1, 2), -c(1, 2)]);mat=apply(mat, 2, as.numeric)
# rownames(mat) <- toupper(dt[-c(1, 2), 1]);colnames(mat) <- toupper(dt[1, -c(1, 2)])
# genetic_essential_mat <- as.matrix(mat)
# print(nrow(genetic_essential_mat))
# dt <- fread("./Data File S3. Genetic interaction profile similarity matrices/Nonessential_Similarity_Network.txt", header = FALSE);dt=as.data.frame(dt)
# mat <- as.matrix(dt[-c(1, 2), -c(1, 2)]);mat=apply(mat, 2, as.numeric)
# rownames(mat) <- toupper(dt[-c(1, 2), 1]);colnames(mat) <- toupper(dt[1, -c(1, 2)])
# genetic_nonessential_mat <- as.matrix(mat)
# print(nrow(genetic_nonessential_mat))
rm(dt, mat)
gc()
## pleiotropic genes
genetic_pleiotropic_df=read_xlsx("./Data File S7_Pleiotropic gene analysis.xlsx",col_names = TRUE)


# check their existence  ratio in matrix
sum(close_genes %in% rownames(genetic_all_mat))/length(close_genes)
sum(eta_genes %in% rownames(genetic_all_mat))/length(eta_genes)
sum(xc_epsilon_genes %in% rownames(genetic_all_mat))/length(xc_epsilon_genes)
sum(scaling_genes %in% rownames(genetic_all_mat))/length(scaling_genes)
# sub-matrices
eta_mat=genetic_all_mat[intersect(eta_genes, rownames(genetic_all_mat)), ,drop = FALSE ]
xc_epsilon_mat=genetic_all_mat[intersect(xc_epsilon_genes, rownames(genetic_all_mat)), ,drop = FALSE ]
scaling_mat=genetic_all_mat[intersect(scaling_genes, rownames(genetic_all_mat)), ,drop = FALSE ]
xc_mat=genetic_all_mat[intersect(xc_genes, rownames(genetic_all_mat)), ,drop = FALSE ]
epsilon_mat=genetic_all_mat[intersect(epsilon_genes, rownames(genetic_all_mat)), ,drop = FALSE ]

# filter out genes with edges
find_first_neighbors <- function(eta_mat){
  edge_cutoff =0.2 # as described in literature PCC=0.2
  eta_edges_list=list()
  for (i in 1:nrow(eta_mat)){
    temp=na.omit(eta_mat[i,]);
    temp=names(temp[temp>edge_cutoff])
    eta_edges_list[[i]]=intersect(temp,distance_df$gene)
    names(eta_edges_list)[i]=rownames(eta_mat)[i]
  }
  distance_eta_df=subset(distance_df, gene %in% c(unlist(eta_edges_list),rownames(eta_mat)))
  rownames(distance_eta_df) <- distance_eta_df$gene
  distance_eta_df$first_neighbor <- NA
  distance_eta_df$is_root <- NA
  for (i in 1:nrow(eta_mat)){
    gene <- rownames(eta_mat)[i]
    gene_list <- eta_edges_list[[i]]
    for (gene_aff in gene_list){
      distance_eta_df[gene_aff, "first_neighbor"] <- gene
      distance_eta_df[gene_aff, "is_root"] <- "first_neighbor"
    }
    distance_eta_df[gene, "first_neighbor"] <- gene
    distance_eta_df[gene, "is_root"] <- "root"
  }
  return(distance_eta_df)
}

plot_neighbors <- function(df,name){
  library(ggsci)
  edges <- df %>%
    filter(is_root != "root") %>%
    left_join(df %>% select(gene, x = mean_rela, y = steepness_rela),
              by = "gene") %>%
    rename(x_gene = x, y_gene = y) %>%
    left_join(df %>% select(gene, x = mean_rela, y = steepness_rela),
              by = c("first_neighbor" = "gene")) %>%
    rename(x_neighbor = x, y_neighbor = y)
  
  ggplot(df, aes(x = mean_rela, y = steepness_rela)) +
    geom_hline(yintercept=1, color = "grey",linewidth = 0.3)+
    geom_vline(xintercept=1, color = "grey",linewidth = 0.3)+
    geom_line(data=data.frame(x = seq(0, 2, length.out = 1000)),aes(x=x, y = eta_x(x)),size = 0.3) +
    geom_line(data=data.frame(x = seq(0, 2, length.out = 1000)),aes(x=x, y = epsilon_x(x)),size = 0.3) +
    geom_line(data=data.frame(x = seq(0, 2, length.out = 1000)),aes(x=x, y = xc_x(x)),size = 0.3) +
    geom_point(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,shape=para),
               size=0.75,alpha=0.5, show.legend = FALSE)+
    scale_shape_manual(values = c(1, 2, 3), name = "simu para") +
    geom_segment(data = edges,
                 aes(x = x_gene, y = y_gene,
                     xend = x_neighbor, yend = y_neighbor),
                 alpha = 0.25,size=0.25, show.legend = FALSE) +
    geom_point(aes(color = first_neighbor, size = is_root)) +
    # scale_color_frontiers()+
    scale_size_manual(values = c(1, 2), name = "is root") +
    # geom_text(aes(label = gene), vjust = -1, size = 1.5,alpha=0.5,
    #           nudge_x = 0.01, nudge_y = 0.01,
    #           segment.size = 0.1, max.overlaps = 20, box.padding = 0.1,
    #           show.legend = FALSE) +
    xlim(0,1.8)+ylim(0.2,1.6)+
    coord_fixed(ratio=1) +
    labs(x = "Relative lifespan", y = "Relative steepness",
         title = paste("first neighbor of ", name)) +
    # remove grid lines
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    )
}

distance_scaling_df = find_first_neighbors(scaling_mat)
distance_xc_epsilon_df = find_first_neighbors(xc_epsilon_mat)
distance_xc_df=find_first_neighbors(xc_mat)
distance_epsilon_df=find_first_neighbors(epsilon_mat)
plot_neighbors(distance_scaling_df,"scaling genes")
plot_neighbors(distance_xc_epsilon_df,"xc_epsilon genes")
plot_neighbors(distance_epsilon_df,"epsilon genes")
plot_neighbors(distance_xc_df,"xc genes")







