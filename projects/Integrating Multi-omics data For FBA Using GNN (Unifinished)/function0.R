# 输入：counts_df（数据框，第一列为基因名，后续列为样本counts）
# 输出：Log2(CPM + 1) 数据框

calculate_log2cpm <- function(counts_df) {
  # 提取counts矩阵（排除基因名列）
  counts_matrix <- as.matrix(counts_df[, -1])
  
  # 计算CPM = (counts / 样本总counts) * 1e6
  cpm <- sweep(counts_matrix, 2, colSums(counts_matrix), "/") * 1e6
  
  # 计算Log2(CPM + 1)
  log2cpm_plus1 <- log2(cpm + 1)
  
  # 合并基因名并返回
  result <- data.frame(
    Gene = counts_df[, 1],
    log2cpm_plus1,
    check.names = FALSE
  )
  return(result)
}
