library(ggplot2)
library(ggrepel)
#BiocManager::install("slingshot")
library(slingshot)
#BiocManager::install("SingleCellExperiment",force = TRUE)
library(SingleCellExperiment)

# 定义突变顺序（从正常到复杂突变）
mutation_order <- c(
  "raw",
  "IDH1", 
  "IDH2",
  "DNMT3A",
  "IDH12",
  "IDH12DNMT3A"
)
library(slingshot)
library(ggplot2)
library(ggrepel)

# 创建SingleCellExperiment对象
sce <- SingleCellExperiment(
  assays = list(logcounts = t(scale(tmp))),  # 注意转置为（基因×细胞）格式
  colData = data.frame(
    mut_state = rownames(tmp),
    stage = factor(rownames(tmp), levels = valid_mutations)  # 保留顺序
  )
)

# 添加PCA结果（直接使用您已计算的pca_result）
reducedDim(sce, "PCA") <- pca_result$x[,1:2]

# 确保mut_state是有序因子（按突变顺序）
colData(sce)$mut_state <- factor(
  colData(sce)$mut_state,
  levels = c("raw", "IDH1", "IDH2", "DNMT3A", "IDH12", "IDH12DNMT3A", "V3","V4","V5","V6","V2")  # 包含所有可能状态
)

# 验证因子水平
levels(colData(sce)$mut_state)  # 应显示预设的顺序

# Slingshot轨迹推断
set.seed(123)
sce <- slingshot(
  sce,
  clusterLabels = "mut_state",
  reducedDim = "PCA",
  start.clus = "raw",
  end.clus = "V2",
  omega = 5,          # 增加曲线平滑度（默认1）
  stretch = 0,        # 禁用曲线拉伸
  thresh = 0.001,     # 降低收敛阈值
  maxit = 50          # 增加最大迭代次数
)
