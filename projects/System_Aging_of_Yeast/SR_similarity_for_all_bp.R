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

library(data.table)
lifes_set_alpha_bp_df <- fread("./data/lifes_set_alpha_bp_df.csv", sep="," ,header=TRUE, stringsAsFactors = FALSE)
lifes_set_alpha_bp_df=as.data.frame(lifes_set_alpha_bp_df)
gc()

### find root for every gene, every bootstrap sample
# distance_all_bp_df <- data.frame(gene = lifes_set_alpha_bp_df$gene,
#                           eta = NA, epsilon = NA, xc = NA)
# search_interval <- c(-2, 6)
# n = 25
# library(rootSolve)
# for (i in 1:length(lifes_set_alpha_bp_df$gene)) {
#   x0 <- lifes_set_alpha_bp_df$mean_rela_bootstrap[i]
#   y0 <- lifes_set_alpha_bp_df$steepness_rela_bootstrap[i]

#   # 使用tryCatch避免错误中断
#   all_roots_eta <- tryCatch({
#     roots <- uniroot.all(f = eta_eq, interval = search_interval, n = n, x0 = x0, y0 = y0)
#     if(length(roots) == 0) NA else roots
#   }, error = function(e) NA)
  
#   all_roots_epsilon <- tryCatch({
#     roots <- uniroot.all(f = epsilon_eq, interval = search_interval, n = n, x0 = x0, y0 = y0)
#     if(length(roots) == 0) NA else roots
#   }, error = function(e) NA)
  
#   all_roots_xc <- tryCatch({
#     roots <- uniroot.all(f = xc_eq, interval = search_interval, n = n, x0 = x0, y0 = y0)
#     if(length(roots) == 0) NA else roots
#   }, error = function(e) NA)
  
#   eta_dis <- ifelse(!is.na(all_roots_eta), min(sqrt((eta_x(all_roots_eta) - y0)^2 + (all_roots_eta - x0)^2)), NA)
#   epsilon_dis <- ifelse(!is.na(all_roots_epsilon), min(sqrt((epsilon_x(all_roots_epsilon) - y0)^2 + (all_roots_epsilon - x0)^2)), NA)
#   xc_dis <- ifelse(!is.na(all_roots_xc), min(sqrt((xc_x(all_roots_xc) - y0)^2 + (all_roots_xc - x0)^2)), NA)
#   distance_all_bp_df$gene[i] <- lifes_set_alpha_bp_df$gene[i]
#   distance_all_bp_df$eta[i] <- eta_dis
#   distance_all_bp_df$epsilon[i] <- epsilon_dis
#   distance_all_bp_df$xc[i] <- xc_dis
#   # display the progress
#   if (i %% 10 == 0) {
#     print(paste("Processed gene:", lifes_set_alpha_bp_df$gene[i], " ", i, "/", length(lifes_set_alpha_bp_df$gene)))
#   }
# }

# ## 更快的优化版本 - 使用向量化网格搜索 + 缓存
# ## ⚡⚡ 超级快速版本：向量化 + 高效并行 + 自适应网格
# library(rootSolve)
# library(parallel)

# # 预计算导数（只计算一次）
# eta_x_prime <- Deriv(eta_x)
# epsilon_x_prime <- Deriv(epsilon_x)
# xc_x_prime <- Deriv(xc_x)

# # 定义方程（避免重复创建函数对象）
# eta_eq_fast <- function(x, x0, y0) {
#   x - x0 + eta_x_prime(x) * (eta_x(x) - y0)
# }

# epsilon_eq_fast <- function(x, x0, y0) {
#   x - x0 + epsilon_x_prime(x) * (epsilon_x(x) - y0)
# }

# xc_eq_fast <- function(x, x0, y0) {
#   x - x0 + xc_x_prime(x) * (xc_x(x) - y0)
# }

# # 超快速根查找：直接使用 uniroot.all（参数调优）
# find_roots_fast <- function(func, interval, x0, y0, n = 300) {
#   tryCatch({
#     roots <- uniroot.all(f = func, interval = interval, n = n, x0 = x0, y0 = y0)
#     if (length(roots) == 0) return(numeric(0))
#     return(roots)
#   }, error = function(e) numeric(0))
# }

# # 计算距离（矢量化）
# calc_distance <- function(roots, x0, y0, func) {
#   if (length(roots) == 0) return(NA_real_)
#   distances <- sqrt((func(roots) - y0)^2 + (roots - x0)^2)
#   return(min(distances))
# }

# # 处理单行数据（可在并行中使用）
# process_single_row <- function(i, x0_vec, y0_vec, gene_vec) {
#   x0 <- x0_vec[i]
#   y0 <- y0_vec[i]
  
#   roots_eta <- find_roots_fast(eta_eq_fast, c(-2, 6), x0, y0, n = 100)
#   roots_epsilon <- find_roots_fast(epsilon_eq_fast, c(-2, 6), x0, y0, n = 100)
#   roots_xc <- find_roots_fast(xc_eq_fast, c(-2, 6), x0, y0, n = 100)
  
#   return(c(
#     gene = gene_vec[i],
#     eta = calc_distance(roots_eta, x0, y0, eta_x),
#     epsilon = calc_distance(roots_epsilon, x0, y0, epsilon_x),
#     xc = calc_distance(roots_xc, x0, y0, xc_x)
#   ))
# }

# # ==================== 数据准备 ====================
# x0_vec <- lifes_set_alpha_bp_df$mean_rela_bootstrap
# y0_vec <- lifes_set_alpha_bp_df$steepness_rela_bootstrap
# gene_vec <- lifes_set_alpha_bp_df$gene
# total_rows <- length(gene_vec)


# ptm <- proc.time()

# # ==================== 初始化结果向量（关键：使用 NA_real_）====================
# eta_results <- NA_real_[1:total_rows]
# epsilon_results <- NA_real_[1:total_rows]
# xc_results <- NA_real_[1:total_rows]

# # ==================== 主循环：逐行处理 ====================
# for (i in 1:total_rows) {
#   x0 <- x0_vec[i]
#   y0 <- y0_vec[i]
  
#   # 找根并计算距离
#   roots_eta <- find_roots_fast(eta_eq_fast, c(-2, 6), x0, y0, n = 300)
#   if (length(roots_eta) > 0) {
#     eta_results[i] <- min(sqrt((eta_x(roots_eta) - y0)^2 + (roots_eta - x0)^2))
#   }
  
#   roots_epsilon <- find_roots_fast(epsilon_eq_fast, c(-2, 6), x0, y0, n = 300)
#   if (length(roots_epsilon) > 0) {
#     epsilon_results[i] <- min(sqrt((epsilon_x(roots_epsilon) - y0)^2 + (roots_epsilon - x0)^2))
#   }
  
#   roots_xc <- find_roots_fast(xc_eq_fast, c(-2, 6), x0, y0, n = 300)
#   if (length(roots_xc) > 0) {
#     xc_results[i] <- min(sqrt((xc_x(roots_xc) - y0)^2 + (roots_xc - x0)^2))
#   }
  
#   # 进度报告（每1000行一次）
#   if (i %% 1000 == 0) {
#     elapsed <- (proc.time() - ptm)[3]
#     per_row <- elapsed / i
#     remaining <- (total_rows - i) * per_row
#     cat(sprintf("✓ %d/%d (%.1f%%) | %.3fs/行 | 剩余 %.1f分钟\n",
#                 i, total_rows, 100 * i / total_rows, per_row, remaining / 60))
#   }
# }

# total_time <- (proc.time() - ptm)[3]

# # ==================== 构建结果数据框 ====================
# distance_all_bp_df <- data.frame(
#   gene = gene_vec,
#   eta = eta_results,
#   epsilon = epsilon_results,
#   xc = xc_results,
#   stringsAsFactors = FALSE
# )

# # write results to CSV
# write.csv(distance_all_bp_df, "./data/distance_all_bp_df.csv", row.names = FALSE)

# cat(sprintf("\n✨ 完成！总耗时: %.1f秒 (%.4f秒/行)\n", total_time, total_time / total_rows))

# # 统计结果
# non_na_count <- sum(!is.na(distance_all_bp_df$eta))
# cat(sprintf("✓ 找到根的行数: %d/%d (%.1f%%)\n\n", non_na_count, total_rows, 100 * non_na_count / total_rows))

distance_all_bp_df <- fread("./data/distance_all_bp_df.csv", sep="," ,header=TRUE, stringsAsFactors = FALSE)
distance_all_bp_df=as.data.frame(distance_all_bp_df)


# add the distance to scaling line
# distance_all_bp_df=merge(distance_all_bp_df, lifes_set_alpha_bp_df, by="gene", all.x=TRUE)
distance_all_bp_df[,c("mean_rela","steepness_rela")] <- lifes_set_alpha_bp_df[, c("mean_rela_bootstrap","steepness_rela_bootstrap")]
head(distance_all_bp_df)
distance_all_bp_df$scaling=abs(distance_all_bp_df$steepness_rela-1)
distance_all_bp_df$length=sqrt((distance_all_bp_df$steepness_rela-1)^2 + (distance_all_bp_df$mean_rela-1)^2)
# relative distance
distance_all_bp_df$eta_rela=distance_all_bp_df$eta/distance_all_bp_df$length
distance_all_bp_df$epsilon_rela=distance_all_bp_df$epsilon/distance_all_bp_df$length
distance_all_bp_df$xc_rela=distance_all_bp_df$xc/distance_all_bp_df$length
distance_all_bp_df$scaling_rela=distance_all_bp_df$scaling/distance_all_bp_df$length
gc()
head(distance_all_bp_df)
write.csv(distance_all_bp_df, "./data/distance_all_bp_df.csv", row.names = FALSE)
