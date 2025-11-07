umap_cluter <- function(tmp,tissue_name=NA){
  ######################### UMAP #########################
  # 加载包
  library(umap)
  library(ggplot2)
  library(ggrepel)
  
  
  # 生成分组标签（根据新的行号定义）
  group_labels <- rep(NA, nrow(tmp))  # 初始化分组向量
  group_labels[1] <- "normal_tissue"
  group_labels[2:5] <- "Single_Pre_leukemic_Initiating_Mutations"
  group_labels[6:10] <- "Several_Pre_leukemic_Initiating_Mutations"
  group_labels[11:14] <- "Additional_Mutations"
  group_labels[15:nrow(tmp)] <- "Cancer_tissue"
  group_levels <- c("normal_tissue",
                    "Single_Pre_leukemic_Initiating_Mutations",
                    "Several_Pre_leukemic_Initiating_Mutations",
                    "Additional_Mutations",
                    "Cancer_tissue")
  
  # 设置UMAP参数（示例：平衡模式）
  umap_config <- umap.defaults
  umap_config$n_neighbors <- 10
  umap_config$min_dist <- 0.05
  umap_config$random_state <- 123
  umap_result <- umap(tmp,umap_config)
  
  color_palette <- c(
    "normal_tissue" = "#4DAF4A",          # 绿色
    "Single_Pre_leukemic_Initiating_Mutations" = "#E41A1C",    # 红色
    "Several_Pre_leukemic_Initiating_Mutations" = "#FF7F00",   # 橙色
    "Additional_Mutations" = "#377EB8",   # 蓝色
    "Cancer_tissue" = "#984EA3"           # 紫色
  )
  
  # 构建绘图数据框（含样本名称和分组）
  umap_df <- data.frame(
    UMAP1 = umap_result$layout[,1],
    UMAP2 = umap_result$layout[,2],
    Sample = rownames(tmp)  # 假设样本名在PCA结果的行名中
  )
  umap_df$Group <- factor(group_labels, levels = group_levels)
  
  # 增强型UMAP可视化
  p1=ggplot(umap_df, aes(x = UMAP1, y = UMAP2, 
                      color = Group)) +  # 添加形状区分
    geom_point(size = 2, alpha = 0.75, stroke = 0.8) +
    geom_text_repel(
      aes(label = Sample),
      size = 2.5,
      box.padding = 0.35,
      max.overlaps = 25,
      segment.color = "grey80",
      min.segment.length = 0.2,
      force = 1.5
    ) +
    scale_color_manual(values = color_palette) +
    labs(
      x = "UMAP1",
      y = "UMAP2"
    ) +
    theme_classic() +
    theme(
      text = element_text(size = 4),
      axis.title = element_text(face = "bold", size = 5),
      legend.position = "right",
      legend.key.size = unit(0.8, "lines"),
      legend.spacing.y = unit(0.2, "cm"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
    ) +
    guides(
      color = guide_legend(
        title = "Biological Groups",
        override.aes = list(size = 3, alpha = 1)
      )
    )+
    # 标题
    ggtitle(paste0("UMAP: ",tissue_name))
  
  
  ######################## DBSCAN聚类 ########################
  library(ggplot2)
  library(ggrepel)
  library(KernSmooth)
  set.seed(123)
  dbscan_res <- dbscan(umap_df[, c("UMAP1", "UMAP2")], 
                       eps = 0.5,       # 根据数据调整
                       minPts = 4)      # 根据数据调整
  
  # 将聚类结果加入数据框（噪声点标记为0）
  umap_df$Cluster <- factor(
    ifelse(dbscan_res$cluster == 0, 
           "Outlier", 
           paste0("Cluster", dbscan_res$cluster))
  )
  
  shape_mapping <- c(
    "Outlier" = 4,      # 灰色叉号
    "Cluster1" = 16,    # 圆形
    "Cluster2" = 17,    # 三角形
    "Cluster3" = 15,     # 方形
    "Cluster4" = 18,    # 菱形
    "Cluster5" = 8,     # 六边形
    "Cluster6" = 9,     # 五角星
    "Cluster7" = 10,    # 十字形
    "Cluster8" = 11    # 加号
  )
  
  # 核密度估计
  dens_est <- bkde2D(umap_df[, c("UMAP1", "UMAP2")], 
                     bandwidth = c(0.5, 0.5), 
                     gridsize = c(100, 100))
  dens_df <- expand.grid(UMAP1 = dens_est$x1, UMAP2 = dens_est$x2)
  dens_df$Density <- as.vector(dens_est$fhat)

  p2=ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
    # 绘制密度等高线
    geom_contour(data = dens_df, 
                 aes(x = UMAP1, y = UMAP2, z = Density),
                 color = "grey90",
                 bins = 4,  # 等高线层数
                 linewidth = 0.5)+
    # 绘制所有点（颜色=生物学分组，形状=聚类）
    geom_point(
      aes(color = Group, 
          shape = Cluster),
      size = 2, 
      alpha = 0.7,
      stroke = 0.8     # 形状边框粗细
    ) +
    # 强制噪声点显示为灰色叉号
    geom_point(
      data = subset(umap_df, Cluster == "Outlier"),
      aes(x = UMAP1, y = UMAP2),
      shape = 4,        # 叉号
      color = "grey50", # 灰色
      size = 2,
      stroke = 0.8
    ) +
    
    # 标签标注（可选）
    geom_text_repel(
      aes(label = Sample,color=Group),
      size = 2,
      box.padding = 0.3,
      max.overlaps = 20,
      segment.color = "grey35",
      segment.size = 0.1,
    ) +
    
    # 颜色和形状映射
    scale_color_manual(values = color_palette) +
    scale_shape_manual(
      values = shape_mapping,
      guide = guide_legend(override.aes = list(color = "black"))) + # 形状图例用黑色
    labs(
      x = "UMAP1",
      y = "UMAP2"
    ) +
    theme_classic() +
    theme(
      text = element_text(family = "Times New Roman", size = 8),
      axis.title = element_text(face = "bold"),
      legend.position = "right",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
    )+
    # 标题
    ggtitle(paste0("UMAP Clustering with DBSCAN: ",tissue_name))
  a=list(p1,p2)
  p2
  return(a)
}







umap_cluter <- function(tmp, tissue_name = NA) {
  # Load necessary libraries
  library(umap)
  library(ggplot2)
  library(ggrepel)
  library(dbscan) # Make sure dbscan is explicitly loaded for dbscan function
  library(KernSmooth) # For bkde2D
  library(dplyr) # For filter and other data manipulation
  
  ######################### UMAP Visualization (Biological Groups) #########################
  
  # Generate group labels based on row indices
  group_labels <- rep(NA, nrow(tmp))
  group_labels[1] <- "Normal Tissue"
  group_labels[2:5] <- "Single Pre-leukemic Mutations"
  group_labels[6:10] <- "Several Pre-leukemic Mutations"
  group_labels[11:14] <- "Additional Mutations"
  group_labels[15:nrow(tmp)] <- "Cancer Tissue"
  
  # Define factor levels to ensure consistent ordering in legend and plots
  group_levels <- c("Normal Tissue",
                    "Single Pre-leukemic Mutations",
                    "Several Pre-leukemic Mutations",
                    "Additional Mutations",
                    "Cancer Tissue")
  
  # Set UMAP parameters (example: balanced mode, adjust as needed)
  umap_config <- umap.defaults
  umap_config$n_neighbors <- 10
  umap_config$min_dist <- 0.05
  umap_config$random_state <- 123
  umap_result <- umap(tmp, umap_config)
  
  # Define a publication-friendly color palette
  color_palette <- c(
    "Normal Tissue" = "#4DAF4A",                       # Green
    "Single Pre-leukemic Mutations" = "#E41A1C",       # Red
    "Several Pre-leukemic Mutations" = "#FF7F00",      # Orange
    "Additional Mutations" = "#377EB8",                # Blue
    "Cancer Tissue" = "#984EA3"                        # Purple
  )
  
  # Build plotting data frame (including sample names and groups)
  umap_df <- data.frame(
    UMAP1 = umap_result$layout[, 1],
    UMAP2 = umap_result$layout[, 2],
    Sample = rownames(tmp), # Assuming sample names are in tmp's row names
    Group = factor(group_labels, levels = group_levels) # Add Group directly here
  )
  
  # *** CORRECTED: Create a subset for labeling, excluding "Cancer Tissue" ***
  # Now 'labels_df' correctly includes the 'Group' column from the start.
  labels_df <- umap_df %>% filter(Group != "Cancer Tissue")
  
  # UMAP Plot 1: Colored by Biological Groups
  p1 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
    geom_point(size = 3, alpha = 0.8, stroke = 0.8) +
    # Use the filtered data for geom_text_repel
    geom_text_repel(
      data = labels_df, # Only label non-Cancer Tissue samples
      aes(label = Sample), # No need for 'color=Group' here as it's already mapped globally
      size = 3,
      box.padding = 0.4,
      point.padding = 0.4,
      max.overlaps = Inf,
      segment.color = "grey60",
      segment.size = 0.3,
      min.segment.length = 0.1,
      force = 2,
      show.legend = FALSE
    ) +
    scale_color_manual(values = color_palette) +
    labs(
      x = "UMAP1",
      y = "UMAP2",
      #title = paste0("UMAP Projection of Samples", ifelse(!is.na(tissue_name), paste0(" (", tissue_name, ")"), "")),
      color = "Biological Group",
      #caption = "Fig. XA: UMAP projection of samples, colored by predefined biological groups. Samples are labeled for identification, excluding 'Cancer Tissue' for clarity."
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.caption = element_text(hjust = 0, size = 10, face = "italic"),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      legend.key.size = unit(1, "lines"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      panel.grid.major = element_line(linewidth = 0.3, color = "grey95"),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = 0.6, color = "black")
    ) +
    guides(
      color = guide_legend(
        title = "Biological Group",
        override.aes = list(size = 4, alpha = 1, stroke = 1)
      )
    )
  
  ######################## DBSCAN Clustering (UMAP-based) ########################
  
  set.seed(123) # For reproducibility of DBSCAN
  dbscan_res <- dbscan(umap_df[, c("UMAP1", "UMAP2")],
                       eps = 0.5,
                       minPts = 4)
  
  # Add clustering result to data frame (noise points marked as 0)
  umap_df$Cluster <- factor(
    ifelse(dbscan_res$cluster == 0,
           "Outlier",
           paste0("Cluster ", dbscan_res$cluster))
  )
  
  # *** NEW LABELS_DF FOR P2: Re-filter after Cluster is added, keeping Group column for color aesthetic ***
  labels_df_p2 <- umap_df %>% filter(Group != "Cancer Tissue")
  
  # Dynamically generate shapes for clusters, ensuring Outlier is mapped to a specific shape
  unique_clusters <- sort(unique(umap_df$Cluster))
  
  if("Outlier" %in% unique_clusters) {
    unique_clusters <- c("Outlier", setdiff(unique_clusters, "Outlier"))
  }
  
  base_shapes <- c(4, 16, 17, 15, 18, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2, 3)
  shape_mapping <- setNames(head(base_shapes, length(unique_clusters)), unique_clusters)
  
  if ("Outlier" %in% names(shape_mapping)) {
    shape_mapping["Outlier"] <- 4
  }
  
  # Kernel Density Estimation for background contours
  dens_est <- bkde2D(umap_df[, c("UMAP1", "UMAP2")],
                     bandwidth = c(0.5, 0.5),
                     gridsize = c(100, 100))
  dens_df <- expand.grid(UMAP1 = dens_est$x1, UMAP2 = dens_est$x2)
  dens_df$Density <- as.vector(dens_est$fhat)
  
  # UMAP Plot 2: DBSCAN Clustering with Density Contours
  p2 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
    # Density contour lines
    geom_contour(data = dens_df,
                 aes(x = UMAP1, y = UMAP2, z = Density),
                 color = "grey80",
                 bins = 5,
                 linewidth = 0.6) +
    # Plot all points, colored by Group, shaped by Cluster
    geom_point(
      aes(color = Group,
          shape = Cluster),
      size = 3,
      alpha = 0.8,
      stroke = 0.8
    ) +
    # Use the filtered data for geom_text_repel
    geom_text_repel(
      data = labels_df_p2, # Pass the filtered data frame here
      aes(label = Sample, color = Group), # 'Group' is now correctly present in labels_df_p2
      size = 3,
      box.padding = 0.4,
      point.padding = 0.4,
      max.overlaps = Inf,
      segment.color = "grey60",
      segment.size = 0.3,
      min.segment.length = 0.1,
      force = 2,
      show.legend = FALSE
    ) +
    # Color and shape mappings
    scale_color_manual(values = color_palette) +
    scale_shape_manual(
      values = shape_mapping,
      guide = guide_legend(title = "DBSCAN Cluster", override.aes = list(color = "black", size = 4, alpha = 1, stroke = 1))
    ) +
    labs(
      x = "UMAP1",
      y = "UMAP2",
      #title = paste0("UMAP with DBSCAN Clustering", ifelse(!is.na(tissue_name), paste0(" (", tissue_name, ")"), "")),
      #subtitle = paste0("DBSCAN parameters: eps = ", umap_config$min_dist, ", minPts = ", umap_config$n_neighbors),
      color = "Biological Group",
      shape = "DBSCAN Cluster",
      #caption = "Fig. XB: UMAP projection with DBSCAN clustering. Density contours (grey lines) indicate regions of higher point density. Points are colored by biological group and shaped by their assigned DBSCAN cluster (X indicates outlier points). Sample labels for 'Cancer Tissue' are hidden for clarity."
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 12),
      plot.caption = element_text(hjust = 0, size = 10, face = "italic"),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      legend.key.size = unit(1, "lines"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      panel.grid.major = element_line(linewidth = 0.3, color = "grey95"),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = 0.6, color = "black")
    ) +
    guides(
      color = guide_legend(title = "Biological Group", override.aes = list(size = 4, alpha = 1, stroke = 1)),
      shape = guide_legend(title = "DBSCAN Cluster", override.aes = list(color = "black", size = 4, alpha = 1, stroke = 1))
    )
  
  return(list(p1_umap_groups = p1, p2_umap_dbscan = p2))
}
