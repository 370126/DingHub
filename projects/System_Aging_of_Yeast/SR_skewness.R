# View(lifes_set_alpha_bp_short_df)

skew_all <- function(x){
  1.1749*x^(-1.6962)-0.1901
}

skew_eta <- function(x){
  1.5209*x^(-1.3955)-0.5324
}

skew_epsilon <- function(x){
  1.1698*x^(-1.7136)-0.1714
}

skew_xc <- function(x){
  1.1172*x^(-1.7441)-0.1510
}

# plot the function and the simulation data
library(ggplot2)
ggplot(data.frame(x = seq(0, 3, length.out = 1000)), aes(x = x)) +
  geom_line(aes(y = skew_eta(x)),size = 0.3) +
  geom_line(aes(y = skew_epsilon(x)),size = 0.3) +
  geom_line(aes(y = skew_xc(x)),size = 0.3) +
  # geom_line(aes(y = skew_all(x)),size = 0.5,color="purple",alpha=0.5) +
  coord_cartesian(xlim=c(0, 3), ylim=c(-5, 9))+
  geom_point(data=simu_long,aes(x=steepness_rela,y=skewness_rela,shape=para),size=1)+
  labs(x = "Relative steepness", y = "Relative skewness",
       title = "Functions fitting saperately of simulation data") +
  # use hollow shapes for the points
  scale_shape_manual(values = c(1, 2, 3)) +
  theme_minimal()

ggplot(data.frame(x = seq(0, 3, length.out = 1000)), aes(x = x)) +
  geom_line(aes(y = skew_all(x)),size = 0.5,color="black",alpha=0.5) +
  coord_cartesian(xlim=c(0, 3), ylim=c(-5, 9))+
  geom_point(data=simu_long,aes(x=steepness_rela,y=skewness_rela,shape=para),size=1)+
  labs(x = "Relative steepness", y = "Relative skewness",
       title = "Functions fitting together of simulation data") +
  # use hollow shapes for the points
  scale_shape_manual(values = c(1, 2, 3)) +
  theme_minimal()

### TEST EXPERIMENTAL DATA
steep <- function(vec){
  vec=as.numeric(vec)
  vec=sort(vec,decreasing=TRUE)
  l=length(vec)
  vec_s=vec
  steepness_vec=1/(sd(vec_s)/mean(vec_s))
  return(steepness_vec)
}

## wild type
# without bootstrap
library(tidyverse)
lifes_WT_alpha_df <- tibble(
  gene = rep(NA, length(yeast_roi_alpha$set_genotype)),
  lifes = vector("list", length(yeast_roi_alpha$set_genotype)),
  steepness_rela = NA_real_,
  mean_lifes_rela = NA_real_,
  skewness_rela = NA_real_
)
for (i in seq_along(yeast_roi_alpha$set_genotype)) {
  lifes_WT_alpha_df$gene[i] <- yeast_roi_alpha$set_genotype[i]
  vec <- as.character(yeast_roi_alpha$ref_lifespans[i])
  vec <- as.numeric(unlist(strsplit(vec, ",")))
  lifes_WT_alpha_df$lifes[[i]] <- vec
  lifes_WT_alpha_df$steepness_rela[i] <- steep(vec)/steepness_WT_alpha
  lifes_WT_alpha_df$mean_lifes_rela[i] <- mean(vec)/mean_lifespan_WT_alpha
  lifes_WT_alpha_df$skewness_rela[i] <- skew(vec)/skewness_WT_alpha
}

# skewness v.s steepness
library(ggplot2)
ggplot(lifes_WT_alpha_df, aes(x = steepness_rela, y = skewness_rela)) +
  geom_vline(xintercept = 1, color = "black",alpha=0.2,size=0.2) +
  geom_hline(yintercept = 1, color = "black",alpha=0.2,size=0.2) +
  geom_line(data=data.frame(x = seq(0, 3, length.out = 1000)),
            aes(x = x, y = skew_all(x)), size = 0.2, color="black", alpha=1) +
  geom_point(data=simu_long,aes(x=steepness_rela,y=skewness_rela,shape=para),
             size=1,alpha=1) +
  scale_shape_manual(values = c(1, 2, 3)) +
  geom_point(size = 1, alpha = 0.2,color="red",stroke=0.5,show.legend = FALSE) +
  # stat_smooth(method = 'lm', se = TRUE, size = 0.75, color="red") +
  coord_cartesian(ylim=c(-5, 9))+
  labs(x = "Relative steepness", y = "Relative skewness",
       title = "Wild type experiment") +
  # remove grid lines
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )




# 
# # skewness v.s lifespan
# ggplot(lifes_WT_alpha_df, aes(x = mean_lifes_rela, y = skewness_rela)) +
#   stat_smooth(data = simu_long,aes(x = lifespan_rela, y = skewness_rela, group = para), 
#               method="gam",se=FALSE,color="black",size=0.2,fullrange = TRUE) +
#   geom_point(data=simu_long,aes(x=lifespan_rela,y=skewness_rela,shape=para),
#              size=1,alpha=1) +
#   scale_shape_manual(values = c(1, 2, 3)) +
#   geom_point(aes(color=gene),size = 0.5, alpha = 0.2,show.legend = FALSE) +
#   xlim(0, 3) + ylim(-3, 10) +
#   labs(x = "Relative lifespan", y = "Relative skewness",
#        title = "Wild type bootstrap results") +
#   theme_minimal()
# # steepness v.s lifespan
# ggplot(lifes_WT_alpha_df, aes(x = mean_lifes_rela, y = steepness_rela)) +
#   geom_line(data=data.frame(x = seq(0, 3, length.out = 1000)),
#             aes(x = x, y = eta_x(x)), size = 0.2, alpha=1,color="black") +
#   geom_line(data=data.frame(x = seq(0, 3, length.out = 1000)),
#             aes(x = x, y = epsilon_x(x)), size = 0.2, alpha=1,color="black") +
#   geom_line(data=data.frame(x = seq(0, 3, length.out= 1000)),
#             aes(x = x, y = xc_x(x)), size = 0.2, alpha=1,color="black") +
#   geom_point(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,shape=para),
#              size=1,alpha=1) +
#   scale_shape_manual(values = c(1, 2, 3)) +
#   geom_point(aes(color=gene),size = 0.5, alpha = 0.2,show.legend = FALSE) +
#   xlim(0, 3) + ylim(0.5, 5) +
#   labs(x = "Relative lifespan", y = "Relative steepness",
#        title = "Wild type bootstrap results") +
#   theme_minimal()

## bootstrap
N_bootstrap=1000;
n=floor(0.1*length(lifes_WT_alpha))
lifes_WT_alpha_bp_df=data.frame(lifespan_rela=rep(NA_real_, N_bootstrap),
                             steepness_rela=rep(NA_real_, N_bootstrap),
                             skewness_rela=rep(NA_real_, N_bootstrap))
for (i in 1:N_bootstrap) {
  # resample with replacement
  sample_indices <- sample(seq_along(lifes_WT_alpha), n, replace = TRUE)
  lifes_WT_alpha_bp_df$lifespan_rela[i] <- mean(lifes_WT_alpha[sample_indices], na.rm = TRUE)/mean_lifespan_WT_alpha
  lifes_WT_alpha_bp_df$steepness_rela[i] <- steep(lifes_WT_alpha[sample_indices])/steepness_WT_alpha
  lifes_WT_alpha_bp_df$skewness_rela[i] <- skew(lifes_WT_alpha[sample_indices])/skewness_WT_alpha
}
# plot
ggplot(lifes_WT_alpha_bp_df, aes(x = steepness_rela, y = skewness_rela)) +
  geom_vline(xintercept = 1, color = "black",alpha=0.2,size=0.2) +
  geom_hline(yintercept = 1, color = "black",alpha=0.2,size=0.2) +
  geom_line(data=data.frame(x = seq(0, 3, length.out = 1000)),
            aes(x = x, y = skew_all(x)), size = 0.2, color="black", alpha=1) +
  geom_point(data=simu_long,aes(x=steepness_rela,y=skewness_rela,shape=para),
             size=1,alpha=1) +
  scale_shape_manual(values = c(1, 2, 3)) +
  geom_point(size = 1, alpha = 0.2,color="grey50",stroke=0.5,show.legend = FALSE) +
  stat_ellipse(level= 0.95, size = 0.2,color="red" ) +
  # stat_smooth(method = 'lm', se = TRUE, size = 0.75, color="red") +
  coord_cartesian(ylim=c(-5, 9))+
  labs(x = "Relative steepness", y = "Relative skewness",
       title = "Wild type bootstrap") +
  # remove grid lines
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )




## experimental
# skewness v.s steepness
label <- function(df,x,y){
  df %>%
    group_by(gene) %>%
    summarize(
      x = 1*min(!!sym(x), na.rm = TRUE),
      y = 1*min(!!sym(y), na.rm = TRUE)
    )
}
# plot
# steepness v.s skewness
label_df= label(subset(lifes_set_alpha_bp_df, gene %in% top20_genes), "steepness_rela_bootstrap", "skewness_rela_bootstrap")
library(ggrepel)
ggplot(subset(lifes_set_alpha_bp_df,gene %in% top20_genes), 
       aes(x = steepness_rela_bootstrap, y = skewness_rela_bootstrap,group=gene,color=gene)) +
  geom_vline(xintercept = 1, color = "black",alpha=0.2,size=0.2) +
  geom_hline(yintercept = 1, color = "black",alpha=0.2,size=0.2) +
  geom_line(data=data.frame(x = seq(0, 3, length.out = 1000)),
            aes(x = x, y = skew_all(x)), size = 0.2, color="black", alpha=1) +
  geom_point(data=simu_long,aes(x=steepness_rela,y=skewness_rela,shape=para),
             size=1,alpha=1,color="black")+
  scale_shape_manual(values = c(1, 2, 3)) +
  geom_point(size = 0.01, alpha = 0.01,show.legend = FALSE) +
  stat_ellipse(level= 0.95, size = 0.4,color="red" ) +
  geom_text_repel(data=label_df, aes(x = x, y = y, label = gene,color=gene),
                  size = 2, show.legend = FALSE)+
  coord_cartesian(xlim=c(0, 3), ylim=c(-5, 9))+
  labs(x = "Relative steepness", y = "Relative skewness",
       title = "Experiments bootstrap results") +
  # remove grid lines
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )

# # skewness v.s lifespan
# label_df= label(subset(lifes_set_alpha_bp_df, gene %in% top10_genes), "mean_rela_bootstrap", "skewness_rela_bootstrap")
# ggplot(subset(lifes_set_alpha_bp_df,gene %in% top10_genes), 
#        aes(x = mean_rela_bootstrap, y = skewness_rela_bootstrap,group=gene,color=gene)) +
#   stat_smooth(data = simu_long,aes(x = lifespan_rela, y = skewness_rela, group = para), 
#             method="gam",se=FALSE,color="black",size=0.2,fullrange = TRUE) +
#   geom_point(data=simu_long,aes(x=lifespan_rela,y=skewness_rela,shape=para),
#              size=1,alpha=1,color="black") +
#   scale_shape_manual(values = c(1, 2, 3)) +
#   geom_point(size = 0.01, alpha = 0.1,show.legend = FALSE) +
#   stat_ellipse(level= 0.95, size = 0.2, color="red" ) +
#   geom_text_repel(data=label_df, aes(x = x, y = y, label = gene,color=gene),
#                   size = 2, show.legend = FALSE)+
#   xlim(0, 3) + ylim(-3, 8) +
#   labs(x = "Relative lifespan", y = "Relative skewness",
#        title = "Experiments bootstrap results") +
#   theme_minimal()
# # steepness v.s lifespan
# label_df= label(subset(lifes_set_alpha_bp_df, gene %in% top10_genes), "mean_rela_bootstrap", "steepness_rela_bootstrap")
# ggplot(subset(lifes_set_alpha_bp_df,gene %in% top10_genes), 
#        aes(x = mean_rela_bootstrap, y = steepness_rela_bootstrap,group=gene,color=gene)) +
#   geom_line(data=data.frame(x = seq(0, 2, length.out = 1000)),
#             aes(x = x, y = eta_x(x)), size = 0.2, alpha=1,color="black") +
#   geom_line(data=data.frame(x = seq(0, 2, length.out = 1000)),
#             aes(x = x, y = epsilon_x(x)), size = 0.2, alpha=1,color="black") +
#   geom_line(data=data.frame(x = seq(0, 2, length.out= 1000)),
#             aes(x = x, y = xc_x(x)), size = 0.2, alpha=1,color="black") +
#   geom_point(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,shape=para),
#              size=1,alpha=1,color="black") +
#   scale_shape_manual(values = c(1, 2, 3)) +
#   geom_point(size = 0.01, alpha = 0.1,show.legend = FALSE) +
#   stat_ellipse(level= 0.95, size = 0.2,color="red") +
#   geom_text_repel(data=label_df, aes(x = x, y = y, label = gene,color=gene),
#                   size = 2, show.legend = FALSE)+
#   xlim(0, 2) + ylim(0.5, 2) +
#   labs(x = "Relative lifespan", y = "Relative steepness",
#        title = "Experiments bootstrap results") +
#   theme_minimal()

    
    
    