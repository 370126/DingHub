draw_cv_sorted <- function(df1, name1, df2, name2, top_number) {
  library(tidyverse)
  library(ggrepel)
  
  # 对每个数据集按 CV 降序排序，并添加排名和来源标记
  df1_sorted <- df1 %>%
    arrange(desc(cv)) %>%
    mutate(
      rank = row_number(),
      source = name1
    )
  
  df2_sorted <- df2 %>%
    arrange(desc(cv)) %>%
    mutate(
      rank = row_number(),
      source = name2
    )
  
  # 合并数据
  combined_data <- bind_rows(df1_sorted, df2_sorted)
  
  # 计算两组数据的平均CV值
  mean_cv1 <- mean(df1_sorted$cv, na.rm = TRUE)
  mean_cv2 <- mean(df2_sorted$cv, na.rm = TRUE)
  
  # 提取每个数据集 CV 最大的前 top_number 名通路
  top5_labels <- combined_data %>%
    group_by(source) %>%
    slice_head(n = top_number) %>%
    ungroup()
  
  # 定义需要特殊标注的通路
  special_pathways <- c(
    "Glycolysis / Gluconeogenesis",
    "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",
    "Glutathione metabolism"
  )
  
  # 提取特殊通路数据
  special_labels <- combined_data %>%
    filter(pathway %in% special_pathways)
  
  # 绘制图形
  ggplot(combined_data, aes(x = rank, y = cv, color = source)) +
    geom_hline(
      yintercept = 0.3,
      color = "red",
      linetype = "dashed",
      size = 0.5
    ) +
    annotate(
      "text",
      x = min(combined_data$rank) * 1.1,
      y = 0.2,
      label = "CV = 0.3",
      color = "red",
      size = 3
    )+
    # # 添加平均CV参考线
    # geom_hline(
    #   yintercept = c(mean_cv1, mean_cv2),
    #   color = c("#E69F00", "#56B4E9"),  # 使用不同颜色区分
    #   linetype = "dashed",
    #   linewidth = 0.8
    # ) +
    # 添加参考线标注
    # annotate(
    #   "text",
    #   x = -Inf,  # 左侧对齐
    #   y = mean_cv1 + 0.05,  # 微调Y轴位置
    #   label = paste0(name1, " mean\nCV = ", round(mean_cv1, 2)),
    #   color = "#E69F00",
    #   size = 2.5,
    #   hjust = -0.1,
    #   vjust = 0
    # ) +
    # annotate(
    #   "text",
    #   x = -Inf,
    #   y = mean_cv2 - 0.05,
    #   label = paste0(name2, " mean\nCV = ", round(mean_cv2, 2)),
    #   color = "#56B4E9",
    #   size = 2.5,
    #   hjust = -0.1,
    #   vjust = 1
    # ) +
    # 数据点
    geom_point(size = 1.25, alpha = 0.5) +
    # 前top_number通路标签
    geom_text_repel(
      data = top5_labels,
      aes(label = pathway),
      size = 2,
      box.padding = 0.3,
      segment.color = "grey50",
      segment.size = 0.3,
      min.segment.length = 0.1,
      show.legend = FALSE
    ) +
    # 特殊通路标签
    geom_text_repel(
      data = special_labels,
      aes(label = paste0(pathway, "\nCV = ", round(cv, 2))),
      size = 2,
      color = "black",
      box.padding = 0.8,
      segment.color = "grey40",
      segment.size = 0.3,
      nudge_x = 0.5,
      nudge_y = 0.5,
      show.legend = FALSE
    ) +
    labs(
      x = "Rank (Descending by CV)",
      y = "Coefficient of Variation",
      color = "Dataset",
      title = "Metabolic Pathway Variability Comparison"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = c(0.8, 0.9),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey50"),
      panel.grid.major = element_line(linewidth = 0.2),
      panel.grid.minor = element_blank()
    )
}

draw_cv_sorted <- function(df1, name1, df2, name2, top_number) {
  library(tidyverse)
  library(ggrepel)
  # You might also consider RColorBrewer for more specific palettes
  # library(RColorBrewer) 
  
  # 对每个数据集按 CV 降序排序，并添加排名和来源标记
  df1_sorted <- df1 %>%
    arrange(desc(cv)) %>%
    mutate(
      rank = row_number(),
      source = name1
    )
  
  df2_sorted <- df2 %>%
    arrange(desc(cv)) %>%
    mutate(
      rank = row_number(),
      source = name2
    )
  
  # 合并数据
  combined_data <- bind_rows(df1_sorted, df2_sorted)
  
  # 提取每个数据集 CV 最大的前 top_number 名通路
  top_labels <- combined_data %>%
    group_by(source) %>%
    slice_head(n = top_number) %>%
    ungroup()
  
  # 定义需要特殊标注的通路
  special_pathways <- c(
    "Glycolysis / Gluconeogenesis",
    "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",
    "Glutathione metabolism"
  )
  
  # 提取特殊通路数据
  special_labels <- combined_data %>%
    filter(pathway %in% special_pathways)
  
  # 绘制图形
  ggplot(combined_data, aes(x = rank, y = cv, color = source)) +
    # Reference line for CV = 0.3
    geom_hline(
      yintercept = 0.3,
      color = "red",
      linetype = "dashed",
      size = 0.6 # Increased line thickness
    ) +
    # Annotation for the reference line
    annotate(
      "text",
      x = max(combined_data$rank) * 0.1, # Position annotation towards the right
      y = 0.3 + 0.05, # Slightly above the line
      label = "CV = 0.3",
      color = "red",
      size = 3.5, # Slightly larger text
      hjust = 1 # Right-align text
    ) +
    # Data points
    geom_point(size = 1.8, alpha = 0.6) + # Increased point size and alpha
    # Labels for top pathways
    geom_text_repel(
      data = top_labels,
      aes(label = pathway),
      size = 3, # Adjusted font size for labels
      box.padding = 0.5,
      point.padding = 0.5,
      segment.color = "grey60",
      segment.size = 0.4,
      min.segment.length = 0.2,
      show.legend = FALSE
    ) +
    # Labels for special pathways
    geom_text_repel(
      data = special_labels,
      aes(label = paste0(pathway, " (CV = ", round(cv, 2), ")")), # More descriptive label
      size = 3,
      color = "black", # Ensure black for special labels
      box.padding = 0.8,
      point.padding = 0.5,
      segment.color = "grey40",
      segment.size = 0.4,
      min.segment.length = 0.2,
      show.legend = FALSE
    ) +
    labs(
      x = "Rank (Descending by Coefficient of Variation)", # More formal x-axis title
      y = "Coefficient of Variation (CV)",
      color = "Data Source", # Clearer legend title
      #title = "Comparison of Metabolic Pathway Variability",
      #caption = "Fig. 1: Variability of metabolic pathways across different datasets. Pathways are ranked by their coefficient of variation (CV) in descending order. The red dashed line indicates a CV of 0.3. Specific pathways are highlighted for their biological relevance." # Placeholder for figure caption
    ) +
    scale_color_manual(values = c("#0072B2", "#D55E00")) + # Example: Colorblind-friendly palette
    # Or use a palette from RColorBrewer:
    # scale_color_brewer(palette = "Set1") + 
    theme_bw(base_size = 14) + # Clean black and white theme, increased base font size
    theme(
      legend.position = c(0.85, 0.9), # Move legend inside the plot
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2), # Add background and border to legend
      legend.title = element_text(face = "bold", size = 12), # Bold and slightly smaller legend title
      legend.text = element_text(size = 11), # Smaller legend text
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16), # Larger, bold, centered title
      plot.caption = element_text(hjust = 0, size = 10, face = "italic"), # Left-aligned, italic caption
      axis.title = element_text(size = 14, face = "bold"), # Bold axis titles
      axis.text = element_text(size = 12), # Axis text size
      panel.grid.major = element_line(linewidth = 0.3, color = "grey90"), # Finer major grid lines
      panel.grid.minor = element_blank(), # Remove minor grid lines
      axis.line = element_line(linewidth = 0.6, color = "black"), # Stronger axis lines
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8) # Stronger panel border
    )
}




draw_cv_sorted_select <- function(df1, name1, df2, name2, top_number, selected_pathways) {
  library(tidyverse)
  library(ggrepel)
  
  # 筛选指定代谢通路 ---------------------------------------------------------------
  filter_pathways <- function(df) {
    df %>% 
      filter(pathway %in% selected_pathways) %>% 
      distinct(pathway, .keep_all = TRUE)
  }
  
  # 预处理数据 --------------------------------------------------------------------
  process_data <- function(df, source_name) {
    df %>% 
      filter_pathways() %>% 
      arrange(desc(cv)) %>% 
      mutate(
        rank = row_number(),
        source = source_name
      )
  }
  
  df1_sorted <- process_data(df1, name1)
  df2_sorted <- process_data(df2, name2)
  
  # 合并数据并检查有效性 -----------------------------------------------------------
  combined_data <- bind_rows(df1_sorted, df2_sorted)
  if (nrow(combined_data) == 0) stop("筛选后无可用数据，请检查selected_pathways参数")
  
  # 计算统计量 --------------------------------------------------------------------
  mean_cv1 <- mean(df1_sorted$cv, na.rm = TRUE)
  mean_cv2 <- mean(df2_sorted$cv, na.rm = TRUE)
  
  # 动态标注设置（修复核心错误）-----------------------------------------------------
  top_labels <- combined_data %>%
    group_by(source) %>%
    group_modify(~ {
      valid_n <- min(top_number, nrow(.x))  # 动态计算有效截断数
      .x %>% slice_head(n = valid_n)        # 安全截断
    }) %>%
    ungroup()
  
  # 关键通路标注（自动匹配selected_pathways中的特殊通路）------------------------------
  key_pathways <- c(
    "Glycolysis / Gluconeogenesis",
    "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",
    "Glutathione metabolism"
  )
  key_labels <- combined_data %>% 
    filter(pathway %in% intersect(key_pathways, selected_pathways))
  
  # 可视化 -----------------------------------------------------------------------
  ggplot(combined_data, aes(x = rank, y = cv, color = source)) +
    geom_hline(
      yintercept = 0.3, 
      color = "red", 
      linetype = "dashed",
      linewidth = 0.5
    ) +
    # geom_hline(
    #   yintercept = c(mean_cv1, mean_cv2),
    #   color = c("#E69F00", "#56B4E9"),
    #   linetype = "dashed",
    #   linewidth = 0.8
    # ) +
    geom_point(
      size = 1.5, 
      alpha = 0.75,
      shape = 16
    ) +
    geom_text_repel(
      data = top_labels,
      aes(label = pathway),
      size = 1,
      box.padding = 0.25,
      segment.color = "grey50",
      segment.size = 0.2,
      max.overlaps = 20,
      show.legend = FALSE
    ) +
    # geom_label_repel(
    #   data = key_labels,
    #   aes(label = paste0(pathway, "\nCV = ", round(cv, 2))),
    #   size = 2,
    #   color = "black",
    #   fill = alpha("white", 0.6),
    #   box.padding = 0.8,
    #   segment.color = "grey40",
    #   segment.size = 0.3,
    #   nudge_x = 0.5,
    #   nudge_y = 0.5,
    #   show.legend = FALSE
    # ) +
    # 特殊通路标签
    geom_text_repel(
      data = key_labels,
      aes(label = paste0(pathway, "\nCV = ", round(cv, 2))),
      size = 2,
      color = "black",
      box.padding = 0.8,
      segment.color = "grey40",
      segment.size = 0.3,
      nudge_x = 0.5,
      nudge_y = 0.5,
      show.legend = FALSE
    ) +
    annotate(
      "text",
      x = -Inf, y = 0.3, 
      label = "CV = 0.3", 
      color = "red",
      size = 3, 
      hjust = -0.1, 
      vjust = -0.5
    ) +
    # annotate(
    #   "text",
    #   x = -Inf, y = mean_cv1,
    #   label = sprintf("%s mean CV: %.2f", name1, mean_cv1),
    #   color = "#E69F00",
    #   size = 3,
    #   hjust = -0.1,
    #   vjust = -0.4
    # ) +
    # annotate(
    #   "text",
    #   x = -Inf, y = mean_cv2,
    #   label = sprintf("%s mean CV: %.2f", name2, mean_cv2),
    #   color = "#56B4E9",
    #   size = 3,
    #   hjust = -0.1,
    #   vjust = -0.4
    # ) +
    # scale_color_manual(values = c("#E69F00", "#56B4E9")) +
    labs(
      x = "Rank (Descending by CV)",
      y = "Coefficient of Variation",
      color = "Dataset",
      title = "Metabolic Pathway Variability Comparison",
      subtitle = paste("Selected pathways:", length(selected_pathways))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = c(0.8, 0.9),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey50"),
      panel.grid.major = element_line(linewidth = 0.2),
      panel.grid.minor = element_blank()
    )
}







draw_cv_sorted_select <- function(df1, name1, df2, name2, top_number, selected_pathways) {
  library(tidyverse)
  library(ggrepel)
  # library(RColorBrewer) # Uncomment if you want to use RColorBrewer palettes
  
  # Filter specified metabolic pathways
  filter_pathways <- function(df) {
    df %>%
      filter(pathway %in% selected_pathways) %>%
      distinct(pathway, .keep_all = TRUE)
  }
  
  # Process data: filter, arrange, and add rank/source
  process_data <- function(df, source_name) {
    df %>%
      filter_pathways() %>%
      arrange(desc(cv)) %>%
      mutate(
        rank = row_number(),
        source = source_name
      )
  }
  
  df1_sorted <- process_data(df1, name1)
  df2_sorted <- process_data(df2, name2)
  
  # Combine data and check for validity
  combined_data <- bind_rows(df1_sorted, df2_sorted)
  if (nrow(combined_data) == 0) {
    stop("No data available after filtering. Please check 'selected_pathways' parameter.")
  }
  
  # Identify top pathways for labeling (dynamic based on available data)
  top_labels <- combined_data %>%
    group_by(source) %>%
    group_modify(~ {
      valid_n <- min(top_number, nrow(.x)) # Ensure 'top_number' doesn't exceed available rows
      .x %>% slice_head(n = valid_n)
    }) %>%
    ungroup()
  
  # Define specific key pathways for highlighting (these will be highlighted if present in selected_pathways)
  key_pathways_to_highlight <- c(
    "Glycolysis / Gluconeogenesis",
    #"Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",
    "Glutathione metabolism"
  )
  
  # Filter for key pathways that are also in the 'selected_pathways' list
  key_labels <- combined_data %>%
    filter(pathway %in% intersect(key_pathways_to_highlight, selected_pathways))
  
  # Visualize
  ggplot(combined_data, aes(x = rank, y = cv, color = source)) +
    # Reference line for CV = 0.3
    geom_hline(
      yintercept = 0.3,
      color = "red",
      linetype = "dashed",
      linewidth = 0.6 # Increased line thickness
    ) +
    # Annotation for the reference line
    annotate(
      "text",
      x = max(combined_data$rank) * 0.1, # Position annotation towards the right
      y = 0.3 + 0.05, # Slightly above the line
      label = "CV = 0.3",
      color = "red",
      size = 3.5, # Slightly larger text
      hjust = 1 # Right-align text
    ) +
    # Data points
    geom_point(
      size = 2, # Increased point size
      alpha = 0.7, # Adjusted alpha for better visibility
      shape = 16 # Solid circles
    ) +
    # Labels for top pathways (dynamic top_number from *selected* pathways)
    geom_text_repel(
      data = top_labels,
      aes(label = pathway),
      size = 3, # Adjusted font size for labels
      box.padding = 0.5,
      point.padding = 0.5,
      segment.color = "grey60",
      segment.size = 0.4,
      min.segment.length = 0.2,
      max.overlaps = Inf, # Allow all labels to be shown, ggrepel handles overlaps
      show.legend = FALSE
    ) +
    # Labels for specific key pathways
    geom_text_repel(
      data = key_labels,
      aes(label = paste0(pathway, " (cv=", round(cv, 2), ")")), # More descriptive label
      size = 3, # Slightly larger for emphasis
      color = "black", # Ensure black for special labels
      #fontface = "bold", # Make key pathway labels bold
      box.padding = 1,
      point.padding = 0.7,
      segment.color = "grey30", # Slightly darker segment for key pathways
      segment.size = 0.5,
      # You might adjust nudge_x/y based on your specific data distribution
      # nudge_x = 0.5,
      # nudge_y = 0.5,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c("#0072B2", "#D55E00")) + # Example: Colorblind-friendly palette
    # Or use a palette from RColorBrewer:
    # scale_color_brewer(palette = "Set1") +
    labs(
      x = "Rank (Descending by Coefficient of Variation)",
      y = "Coefficient of Variation (CV)",
      color = "Data Source",
      #title = "Variability of Selected Metabolic Pathways", # More focused title
      subtitle = paste0("Number of selected pathways: ", length(selected_pathways)), # Informative subtitle
      #caption = "Fig. X: Coefficient of Variation (CV) for a subset of metabolic pathways. Pathways are ranked by CV within each dataset. Highlighted pathways are specific biologically relevant examples." # Placeholder caption
    ) +
    theme_bw(base_size = 14) + # Clean black and white theme, increased base font size
    theme(
      legend.position = c(0.85, 0.9), # Move legend inside the plot
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2), # Add background and border to legend
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 12), # Slightly smaller subtitle
      plot.caption = element_text(hjust = 0, size = 10, face = "italic"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      panel.grid.major = element_line(linewidth = 0.3, color = "grey90"),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = 0.6, color = "black"),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8)
    )
}
