fpkm_to_tpm <- function(fpkm) {
  total_fpkm <- sum(fpkm)
  tpm <- fpkm / total_fpkm * 1e6
  return(tpm)
}



calculate_log2tpm <- function(counts_df, gene_length_df) {
  # # 检查必要列是否存在
  # if (!"gene_name" %in% colnames(counts_df)) {
  #   stop("counts_df 必须包含 'gene_name' 列")
  # }
  # if (!all(c("gene_id", "length") %in% colnames(gene_length_df))) {
  #   stop("gene_length_df 必须包含 'gene_id' 和 'length' 列")
  # }
  # 
  # # 检查重复基因
  # if (any(duplicated(counts_df$gene_name))) {
  #   stop("counts_df 中存在重复基因名")
  # }
  # if (any(duplicated(gene_length_df$gene_id))) {
  #   stop("gene_length_df 中存在重复基因ID")
  # }
  # 
  # 合并数据并匹配基因长度
  merged <- merge(
    counts_df,
    gene_length_df,
    by.x = "gene_name",
    by.y = "gene_id",
    all.x = FALSE
  )
  
  # 检查基因是否全部匹配
  # if (nrow(merged) < nrow(counts_df)) {
  #   missing <- setdiff(counts_df$gene_name, merged$gene_name)
  #   stop("以下基因在 gene_length_df 中缺失: ", paste(head(missing), collapse = ", "))
  # }
  
  # 提取数值矩阵和长度
  sample_cols <- setdiff(colnames(merged), c("gene_name", "length"))
  counts_mat <- as.matrix(merged[sample_cols])
  lengths <- merged$length
  
  # 计算 TPM
  rpk <- counts_mat / (lengths / 1000)          # Reads per kilobase
  tpm <- sweep(rpk, 2, colSums(rpk)/1e6, "/")   # Transcripts per million
  
  # 添加 gene_name 并转换为数据框
  result <- data.frame(
    gene_name = merged$gene_name,
    log2(tpm + 1),
    check.names = FALSE
  )
  
  return(result)
}
