filter_rna <- function(df) {
  df$abundence_adjusted <- df[,2]
  
  # 条件 1: Prob < 0.3 & abundence < 1 → 设为 0
  cond1 <- df$Prob < 0.3 & df[,2] < 1
  df$abundence_adjusted[cond1] <- 0
  
  # 条件 2: Prob < 0.3 & abundence > 3 → 设为 NA
  cond2 <- df$Prob < 0.3 & df[,2] > 3
  df$abundence_adjusted[cond2] <- NA
  
  # 处理 NA 值
  #df[,2]_adjusted[is.na(df$Prob) | is.na(df[,2])] <- NA
  na_count <- sum(is.na(df$abundence_adjusted))
  cat("abundence_adjusted 列中的 NA 数量:", na_count,",",na_count/nrow(df), "\n")
  return(df)
}




handle_epigenetic_transcriptome_conflict <- function(
    df, 
    prob_col = "Prob", 
    abundence_col = "abundence", 
    prob_high_threshold = 0.7, 
    prob_low_threshold = 0.3,
    abundence_high_threshold = NULL,  # 默认自动计算75%分位数
    abundence_low_threshold = NULL,   # 默认自动计算25%分位数
    action = "mark"                   # 可选 "mark", "filter", "integrate"
) {
  # 1. 创建输入数据框的副本，避免修改原始数据
  df_copy <- data.frame(df)
  
  # 2. 检查必要列是否存在
  if (!all(c(prob_col, abundence_col) %in% colnames(df_copy))) {
    stop("输入数据框中必须包含 'Prob' 和 'abundence' 列")
  }
  
  # 3. 自动计算分位数阈值（若未指定）
  if (is.null(abundence_high_threshold) || is.null(abundence_low_threshold)) {
    abundence_quantiles <- quantile(
      df_copy[[abundence_col]], 
      c(0.25, 0.75), 
      na.rm = TRUE
    )
    abundence_low_threshold <- abundence_quantiles[1]
    abundence_high_threshold <- abundence_quantiles[2]
  }
  
  normal_data=normalize_data(df_copy[[abundence_col]])
  # 4. 添加冲突类型和评分列（不修改原数据）
  df_copy <- df_copy %>%
    dplyr::mutate(
      Conflict_Type = dplyr::case_when(
        # Type1: 表观预测激活但表达低
        .data[[prob_col]] >= prob_high_threshold & 
          .data[[abundence_col]] <= abundence_low_threshold ~ "Type1",
        # Type2: 表观预测抑制但表达高
        .data[[prob_col]] <= prob_low_threshold & 
          .data[[abundence_col]] >= abundence_high_threshold ~ "Type2",
        TRUE ~ "No_Conflict"
      ),
      Confidence_Score = abs(normal_data - .data[[prob_col]])
    )
  
  # 5. 根据action参数处理数据
  if (action == "filter") {
    df_copy <- df_copy %>% 
      dplyr::filter(Conflict_Type != "No_Conflict")
  } else if (action == "integrate") {
    df_copy <- df_copy %>%
      dplyr::mutate(
        Integrated_Score = 0.5 * .data[[prob_col]] + 0.5 * .data[[abundence_col]]
      )
  }
  
  # 6. 返回新数据框（原数据未修改）
  return(df_copy)
}




show_co_occurence <- function(type2_symbols_list,name){
# 创建二元矩阵
binary_matrix <- t(sapply(type2_symbols_list, function(x) 
  as.integer(unique(unlist(type2_symbols_list)) %in% x)))
# 添加行名
rownames(binary_matrix) <- names(type2_symbols_list)
# 绘制热图
pheatmap::pheatmap(binary_matrix,
                   cluster_rows = FALSE,cluster_cols = FALSE,
                   color = c("white", "black"),
                   legend_breaks = c(0, 1),
                   legend_labels = c("Absent", "Present"),
                   main = name)
}
