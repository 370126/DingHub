# source functions
source("./SR_fitting.R")
# equations
# (x-x0)+(y(x)-y0)*y'=0
library(Deriv)
eta_eq <- function(x,x0,y0){
  x-x0 + Deriv(eta_x)(x)*(eta_x(x)-y0)
}
epsilon_eq <- function(x,x0,y0){
  x-x0 + Deriv(epsilon_x)(x)*(epsilon_x(x)-y0)
}
xc_eq <- function(x,x0,y0){
  x-x0 + Deriv(xc_x)(x)*(xc_x(x)-y0)
}
### find root for every data point
lifes_set_alpha_bp_short_df=map_df(unique(lifes_set_alpha_bp_df$gene),function(gene){
  mean_mean_rela_bp=mean(lifes_set_alpha_bp_df$mean_rela_bootstrap[lifes_set_alpha_bp_df$gene==gene])
  mean_steepness_rela_bp=mean(lifes_set_alpha_bp_df$steepness_rela_bootstrap[lifes_set_alpha_bp_df$gene==gene])
  mean_skewness_rela_bp=mean(lifes_set_alpha_bp_df$skewness_rela_bootstrap[lifes_set_alpha_bp_df$gene==gene])
  tibble(gene=gene,
         mean_rela=mean_mean_rela_bp,
         steepness_rela=mean_steepness_rela_bp,
         skewness_rela=mean_skewness_rela_bp)
})

# plot with FOB1 x0,y0
eta_eq_plot <- function(x0, y0) {
  x_vals <- seq(-2, 8, length.out = 1000)
  y_vals <- eta_eq(x_vals, x0, y0)
  data.frame(x = x_vals, y = y_vals)
}
epsilon_eq_plot <- function(x0, y0) {
  x_vals <- seq(-2, 8, length.out = 1000)
  y_vals <- epsilon_eq(x_vals, x0, y0)
  data.frame(x = x_vals, y = y_vals)
}
xc_eq_plot <- function(x0, y0) {
  x_vals <- seq(-2, 8, length.out = 1000)
  y_vals <- xc_eq(x_vals, x0, y0)
  data.frame(x = x_vals, y = y_vals)
}
ggplot(eta_eq_plot(0.1308, 4.9324), aes(x = x, y = y)) +
  geom_line(size = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red",size=0.1) +
  ylim(-10,20) +
  labs(x = "x", y = "eta_eq(x)", title = "eta_eq Function with FOB1 Parameters") +
  theme_minimal()
ggplot(epsilon_eq_plot(0.1006, 2.4711), aes(x = x, y = y)) +
  geom_line(size = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red",size=0.1) +
  ylim(-10,20) +
  labs(x = "x", y = "epsilon_eq(x)", title = "epsilon_eq Function with FOB1 Parameters") +
  theme_minimal()
ggplot(xc_eq_plot(0.1006, 2.4711), aes(x = x, y = y)) +
  geom_line(size = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red",size=0.1) +
  ylim(-10,10) +
  labs(x = "x", y = "xc_eq(x)", title = "xc_eq Function with FOB1 Parameters") +
  theme_minimal()

distance_df <- data.frame(gene = lifes_set_alpha_bp_short_df$gene,
                          eta = NA, epsilon = NA, xc = NA)
library(rootSolve)
for (i in 1:length(lifes_set_alpha_bp_short_df$gene)) {
  x0 <- lifes_set_alpha_bp_short_df$mean_rela[i]
  y0 <- lifes_set_alpha_bp_short_df$steepness_rela[i]
  search_interval <- c(-2, 6)
  n = 50
  all_roots_eta <- uniroot.all(
    f = eta_eq,
    interval = search_interval,
    n = n,
    x0 = x0,
    y0 = y0
  )
  if(length(all_roots_eta) == 0) {all_roots_eta <- NA}
  all_roots_epsilon <- uniroot.all(
    f = epsilon_eq,
    interval = search_interval,
    n = n,
    x0 = x0,
    y0 = y0
  )
  if(length(all_roots_epsilon) == 0) {all_roots_epsilon <- NA}
  all_roots_xc <- uniroot.all(
    f = xc_eq,
    interval = search_interval,
    n = n,
    x0 = x0,
    y0 = y0
  )
  if(length(all_roots_xc) == 0) {all_roots_xc <- NA}
  eta_dis <- ifelse(!is.na(all_roots_eta), min(sqrt((eta_x(all_roots_eta) - y0)^2 + (all_roots_eta - x0)^2)), NA)
  epsilon_dis <- ifelse(!is.na(all_roots_epsilon), min(sqrt((epsilon_x(all_roots_epsilon) - y0)^2 + (all_roots_epsilon - x0)^2)), NA)
  xc_dis <- ifelse(!is.na(all_roots_xc), min(sqrt((xc_x(all_roots_xc) - y0)^2 + (all_roots_xc - x0)^2)), NA)
  distance_df$gene[i] <- lifes_set_alpha_bp_short_df$gene[i]
  distance_df$eta[i] <- eta_dis
  distance_df$epsilon[i] <- epsilon_dis
  distance_df$xc[i] <- xc_dis
  # display the progress
  if (i %% 10 == 0) {
    print(paste("Processed gene:", lifes_set_alpha_bp_short_df$gene[i], " ", i, "/", length(lifes_set_alpha_bp_short_df$gene)))
  }
}
# check out rows with NA values
distance_df_NA = distance_df[!complete.cases(distance_df),];length(distance_df_NA$gene)
temp=subset(lifes_set_alpha_bp_short_df, gene %in% distance_df_NA$gene)
library(ggplot2)
ggplot(data.frame(x = seq(-2, 2, length.out = 1000)), aes(x = x)) +
  geom_hline(yintercept=1, color = "black",linewidth = 0.3)+
  geom_line(aes(y = eta_x(x)),size = 0.3) +
  geom_line(aes(y = epsilon_x(x)),size = 0.3) +
  geom_line(aes(y = xc_x(x)),size = 0.3) +
  xlim(0, 2) + ylim(0,3)+
  geom_point(data=simu_long,aes(x=lifespan_rela,y=steepness_rela,shape=para),size=1)+
  geom_point(data=temp,aes(x=mean_rela,y=steepness_rela,shape=gene),size=1.3) +
  coord_fixed(ratio = 1) +
  labs(x = "Relative lifespan", y = "Relative steepness",
       title = "Functions fitting of simulation data") +
  scale_shape_manual(values = c(1,2,3))+
  theme_minimal()
# add the distance to scaling line
distance_df=merge(distance_df, lifes_set_alpha_bp_short_df, by="gene", all.x=TRUE)
distance_df$scaling=abs(distance_df$steepness_rela-1)
distance_df$length=sqrt((distance_df$steepness_rela-1)^2 + (distance_df$mean_rela-1)^2)
# relative distance
distance_df$eta_rela=distance_df$eta/distance_df$length
distance_df$epsilon_rela=distance_df$epsilon/distance_df$length
distance_df$xc_rela=distance_df$xc/distance_df$length
distance_df$scaling_rela=distance_df$scaling/distance_df$length

# plot the histogram of distances on one single figure but separately
feature_df <- distance_df[, c("gene","eta_rela", "epsilon_rela", "xc_rela", "scaling_rela")]
feature_long=pivot_longer(feature_df,
                          cols = !gene,
                          names_to = "to",
                          values_to = "distance_rela")

ggplot(feature_long, aes(x = distance_rela, color = to, fill = to)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Density Distribution of the relative distances to each parameters")
rm(feature_long)
gc()














###################################### PCA & DBSCN by distances
# 
# ## relative distance PCA
# rownames(feature_df) <- distance_df$gene
# feature_df=subset(feature_df,select=-gene)
# pca_result <- prcomp(feature_df, center = FALSE, scale. = FALSE)
# # Convert PCA result to a data frame for plotting
# pca_df <- as.data.frame(pca_result$x)
# pca_df$gene <- rownames(pca_df)
# eigenvalues <- pca_result$sdev^2
# total_variance <- sum(eigenvalues)
# variance_explained_ratio <- eigenvalues / total_variance
# pc1_variance <- paste0(round(variance_explained_ratio[1] * 100, 2), "%")
# pc2_variance <- paste0(round(variance_explained_ratio[2] * 100, 2), "%")
# loadings <- as.data.frame(pca_result$rotation[, 1:2])
# loadings$feature <- rownames(loadings)
# arrow_scale <- 1
# ggplot(pca_df, aes(x = PC1, y = PC2, label = gene)) +
#   geom_point(aes(color = gene), size = 2, alpha=0.5, stroke=0) +
#   geom_segment(data = loadings,
#                aes(x = 0, y = 0,
#                    xend = PC1 * arrow_scale,
#                    yend = PC2 * arrow_scale),
#                arrow = arrow(length = unit(0.2, "cm")),
#                color = "red", size = 0.5) +
#   geom_text(data = loadings,
#             aes(x = PC1 * arrow_scale * 1.1,
#                 y = PC2 * arrow_scale * 1.1,
#                 label = feature),
#             color = "red", size = 3) +
#   labs(title = paste0("PCA of Relative Distances to Scaling Line (PC1: ", pc1_variance, ", PC2: ", pc2_variance, ")"),
#        x = paste0("Principal Component 1 (", pc1_variance, ")"),
#        y = paste0("Principal Component 2 (", pc2_variance, ")")) +
#   theme_minimal() +
#   theme(legend.position = "none")
# 
# ## absolute distance PCA
# rownames(feature_ab_df) <- distance_df$gene
# feature_ab_df=subset(feature_ab_df,select=-gene)
# pca_result_ab <- prcomp(feature_ab_df, center = FALSE, scale. = FALSE)
# # Convert PCA result to a data frame for plotting
# pca_df_ab <- as.data.frame(pca_result_ab$x)
# pca_df_ab$gene <- rownames(pca_df_ab)
# eigenvalues_ab <- pca_result_ab$sdev^2
# total_variance_ab <- sum(eigenvalues_ab)
# variance_explained_ratio_ab <- eigenvalues_ab / total_variance_ab
# pc1_variance_ab <- paste0(round(variance_explained_ratio_ab[1] * 100, 2), "%")
# pc2_variance_ab <- paste0(round(variance_explained_ratio_ab[2] * 100, 2), "%")
# loadings_ab <- as.data.frame(pca_result_ab$rotation[, 1:2])
# loadings_ab$feature <- rownames(loadings_ab)
# arrow_scale_ab <- 1
# ggplot(pca_df_ab, aes(x = PC1, y = PC2, label = gene)) +
#   geom_point(aes(color = gene), size = 2, alpha=0.5, stroke=0) +
#   geom_segment(data = loadings_ab,
#                aes(x = 0, y = 0,
#                    xend = PC1 * arrow_scale_ab,
#                    yend = PC2 * arrow_scale_ab),
#                arrow = arrow(length = unit(0.2, "cm")),
#                color = "red", size = 0.5) +
#   geom_text(data = loadings_ab,
#             aes(x = PC1 * arrow_scale_ab * 1.1,
#                 y = PC2 * arrow_scale_ab * 1.1,
#                 label = feature),
#             color = "red", size = 3) +
#   labs(title = paste0("PCA of Absolute Distances to Scaling Line (PC1: ", pc1_variance_ab, ", PC2: ", pc2_variance_ab, ")"),
#        x = paste0("Principal Component 1 (", pc1_variance_ab, ")"),
#        y = paste0("Principal Component 2 (", pc2_variance_ab, ")")) +
#   theme_minimal() +
#   theme(legend.position = "none")
# 




# 
# ##### clustering by DBSCAN
# DBS <- function (feature_df, k, h){ 
# library(dbscan)
# feature_matrix <- as.matrix(feature_df)  # 转换为矩阵
# # -----------------------------
# # Step 1. 绘制 kNN 距离图（用于手动找拐点）
# # -----------------------------
# k <- k  # minPts
# kNNdistplot(feature_matrix, k = k)
# abline(h = h,col="red", lty = 2)  # 临时参考线
# h_good=h
# # # -----------------------------
# # # Step 2. 自动尝试 eps 值，选择最优聚类效果
# # # -----------------------------
# # # 定义 eps 值的范围（你也可以扩展）
# # eps_values <- seq(0.01, 0.2, by = 0.001)
# 
# # # 记录每个 eps 对应的非噪声簇数
# # results <- data.frame(eps = eps_values, clusters = NA)
# 
# # for (i in seq_along(eps_values)) {
# #   res <- dbscan(feature_matrix, eps = eps_values[i], minPts = k)
# #   num_clusters <- length(unique(res$cluster[res$cluster != 0]))  # 排除噪声
# #   results$clusters[i] <- num_clusters
# # }
# # # 找出非噪声簇最多时对应的 eps（你也可以改成轮廓系数最优）
# # best_eps <- results$eps[which.max(results$clusters)]
# # cat("🌟 最优 eps 值为：", best_eps, "，对应的聚类数为：", max(results$clusters), "\n")
# 
# # -----------------------------
# # Step 3. 用最佳 eps 进行 DBSCAN 聚类
# # -----------------------------
# final_result <- dbscan(feature_matrix, eps = h_good, minPts = k)
# cluster_labels <- factor(final_result$cluster)
# return(cluster_labels)
# }
# cluster_labels <- DBS(feature_df,10,0.06)
# # -----------------------------
# # Step 4. PCA 降维并可视化聚类结果
# # -----------------------------
# pca_df$cluster <- cluster_labels
# ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
#   geom_point(size = 2, alpha = 0.7) +
#   geom_segment(data = loadings,
#                aes(x = 0, y = 0,
#                    xend = PC1 * arrow_scale,
#                    yend = PC2 * arrow_scale),
#                arrow = arrow(length = unit(0.2, "cm")),
#                color = "red", size = 0.5) +
#   geom_text(data = loadings,
#             aes(x = PC1 * arrow_scale * 1.1, 
#                 y = PC2 * arrow_scale * 1.1, 
#                 label = feature),
#             color = "red", size = 3) +
#   labs(title = paste0("PCA of Relative Distances to Scaling Line"),
#        x = paste0("Principal Component 1 (", pc1_variance, ")"),
#        y = paste0("Principal Component 2 (", pc2_variance, ")")) +
#   theme_minimal()
# 
# 
# 
# ###################################### PCA & DBSCN by statistical characteristics
# # review 2-D
# p4=ggplot() +
#   geom_hline(yintercept=1, color = "black",linewidth = 0.7,alpha=0.5)+
#   geom_point(data = lifes_set_alpha_bp_short_df[lifes_set_alpha_bp_short_df$gene %in% top20_genes, ], 
#              aes(x = mean_rela, y = steepness_rela, color = gene), 
#              size = 2, alpha = 0.5) +
#   new_scale_color() +
#   geom_line(data = simu_long, aes(x = lifespan_rela, y = steepness_rela, color = para, group = para), size = 0.7) +
#   # geom_point(data = simu_long, aes(x = lifespan_rela, y = steepness_rela, color = para, group = para, shape = para)) +
#   scale_color_brewer(palette = "Set1", name = "Parameters") +
#   scale_shape_discrete(guide = "none") +
#   xlim(0.2, 2.0) +
#   theme_minimal() +
#   labs(x = "Relative Mean Lifespan", y = "Relative Steepness",
#        title='Simulation of parameter experiments & experimental data (zoom into 3 genes)');p4
# 
# 
# 
# 
# ## plot 3-D graph
# install.packages("plotly")
# 
# library(plotly)
# plot_ly(lifes_set_alpha_bp_short_df ,
#                x = ~mean_rela,
#                y = ~steepness_rela, 
#                z = ~skewness_rela, 
#                color = ~gene, 
#                marker = list(size = 1, opacity = 0.8)) %>%
#   layout(scene = list(xaxis = list(title = 'lifespan'),
#                       yaxis = list(title = 'steepness'),
#                       zaxis = list(title = 'skewness')),
#          title = "lifespan v.s. steepness v.s. skewness") 
# # 20 genes
# plot_ly(lifes_set_alpha_bp_short_df[lifes_set_alpha_bp_short_df$gene %in% top20_genes,] ,
#         x = ~mean_rela,
#         y = ~steepness_rela, 
#         z = ~skewness_rela, 
#         color = ~gene, 
#         marker = list(size = 2, opacity = 0.8)) %>%
#   layout(scene = list(xaxis = list(title = 'lifespan'),
#                       yaxis = list(title = 'steepness'),
#                       zaxis = list(title = 'skewness')),
#          title = "lifespan v.s. steepness v.s. skewness") 
# 
# write.csv(lifes_set_alpha_bp_short_df,"lifes_set_alpha_bp_short_df.csv",row.names = FALSE)
# write.csv(simu_long,"simu_long.csv",row.names = FALSE)

