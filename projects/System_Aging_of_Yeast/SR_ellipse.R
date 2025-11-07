make_ellipse <- function(x, y, level = 0.95, npoints = 100) {
  # 1. 计算均值和协方差
  mu <- colMeans(cbind(x, y))
  sigma <- cov(cbind(x, y))
  
  # 2. 计算卡方分布临界值 (95% CI 椭圆)
  c <- qchisq(level, df = 2)
  
  # 3. 对协方差矩阵做特征分解
  eig <- eigen(sigma)
  order_idx <- order(eig$values, decreasing=TRUE)  # 从大到小排序
  axes <- eig$vectors[, order_idx]                 # 长轴在第一列
  lengths <- sqrt(eig$values[order_idx] * c)
  ellipse=list(axes=axes,lengths=lengths,mu=mu)
  return(ellipse)
}

ellipse_function <- function(ellipse_axes,ellipse_lengths,mu,theta){
  circle <- rbind(cos(theta), sin(theta))  # 单位圆
  ellipse <- t(mu + ellipse_axes %*% diag(ellipse_lengths) %*% circle)
  colnames(ellipse) <- c("x", "y")
  return(as.data.frame(ellipse))
}

cal_rotate <- function(ellipse_axes){
  angle <- atan2(ellipse_axes[2,1], ellipse_axes[1,1])  # 长轴与x轴的夹角
  return(angle)
}



# test
set.seed(1)
x <- rnorm(100)
y <- 1*x + rnorm(100)
ell <- make_ellipse(x, y)
theta_points <- c(0, pi/2, pi, 3*pi/2)  # 0, π/2, π, 3π/2, -angle
points_df <- ellipse_function(ell$axes, ell$lengths, ell$mu, theta_points)
points_df$theta_label <- c("0", "π/2", "π", "3π/2")
theta_seq <- seq(0, 2*pi, length.out = 200)
ellipse_df <- ellipse_function(ell$axes, ell$lengths, ell$mu, theta_seq)
long_right <- as.data.frame(t(ell$mu + ell$axes[,1] * ell$lengths[1]))
long_left  <- as.data.frame(t(ell$mu - ell$axes[,1] * ell$lengths[1]))
short_up   <- as.data.frame(t(ell$mu + ell$axes[,2] * ell$lengths[2]))
short_down <- as.data.frame(t(ell$mu - ell$axes[,2] * ell$lengths[2]))
colnames(long_right) <- colnames(long_left) <- colnames(short_up) <- colnames(short_down) <- c("x","y")
library(ggplot2)
# find cross point of ellipse and line y=tan(θ)*x
# cross should be ...
a <- ell$lengths[1]
b <- ell$lengths[2]
theta_rot <- atan2(ell$axes[2,1], ell$axes[1,1])
theta_line <- pi/6
# alpha*cosθ + beta*sinθ = 0
alpha <- a * sin(theta_rot) - a * cos(theta_rot) * tan(theta_line)
beta  <- b * cos(theta_rot) + b * sin(theta_rot) * tan(theta_line)
theta_cross <- atan2(-alpha, beta)
slope <- tan(theta_line)
intercept <- 0

theta_points <- c(0, pi/2, pi, 3*pi/2, theta_cross, theta_cross+pi)
points_df <- ellipse_function(ell$axes, ell$lengths, ell$mu, theta_points)
points_df$theta_label <- c("0", "π/2", "π", "3π/2", "θ*", "θ*+π")
ggplot() +
  geom_point(aes(x, y), color="grey") +
  geom_path(data=ellipse_df, aes(x, y), color="blue") +
  geom_abline(slope=slope, intercept=intercept, color="black") +
  geom_point(data=points_df, aes(x, y), color="red", size=3) +
  geom_text(data=points_df, aes(x, y, label=theta_label), vjust=-1, color="red") +
  geom_point(data=as.data.frame(t(ell$mu)), aes(x, y), color="green", size=3) +
  geom_text(data=as.data.frame(t(ell$mu)), aes(x, y, label="μ"), vjust=-1, color="green") +
  geom_segment(data=long_right, aes(x=ell$mu[1], y=ell$mu[2], xend=x, yend=y), color="purple", linetype="dashed") +
  geom_segment(data=long_left, aes(x=ell$mu[1], y=ell$mu[2], xend=x, yend=y), color="purple", linetype="dashed") +
  geom_segment(data=short_up, aes(x=ell$mu[1], y=ell$mu[2], xend=x, yend=y), color="orange", linetype="dashed") +
  geom_segment(data=short_down, aes(x=ell$mu[1], y=ell$mu[2], xend=x, yend=y), color="orange", linetype="dashed") +
  coord_equal() +
  theme_minimal()
gc()



library(data.table)
library(readxl)
library(tidyverse)
lifes_set_alpha_bp_df <- fread("./data/lifes_set_alpha_bp_df.csv", sep="," ,header=TRUE, stringsAsFactors = FALSE)
lifes_set_alpha_bp_df=as.data.frame(lifes_set_alpha_bp_df);head(lifes_set_alpha_bp_df)
setDT(lifes_set_alpha_bp_df)
ellipse_df <- lifes_set_alpha_bp_df[, {
  x = mean_rela_bootstrap
  y = steepness_rela_bootstrap
  mu = colMeans(cbind(x, y))
  slope = (mu[2] - 1) / (mu[1] - 1)
  theta_line = atan(slope)
  
  ell = make_ellipse(x, y)
  theta_rot = atan2(ell$axes[2,1], ell$axes[1,1])
  a = ell$lengths[1]
  b = ell$lengths[2]
  alpha = a * sin(theta_rot) - a * cos(theta_rot) * tan(theta_line)
  beta  = b * cos(theta_rot) + b * sin(theta_rot) * tan(theta_line)
  theta_cross = atan2(-alpha, beta)
  # find distance to this direction
  cross=ellipse_function(ell$axes, ell$lengths, ell$mu, c(theta_cross,theta_cross+pi))
  x_cross = cross$x
  y_cross = cross$y
  distance_to_center = sqrt((x_cross - ell$mu[1])^2 + (y_cross - ell$mu[2])^2)[1]
  distance_to_O = sqrt((ell$mu[1] - 1)^2 + (ell$mu[2] - 1)^2)
  mu_x = ell$mu[1]
  mu_y = ell$mu[2]
  .(mu_x=mu_x, mu_y=mu_y, theta_line=theta_line, theta_cross=theta_cross, distance_to_center=distance_to_center, distance_to_O=distance_to_O)
}, by=gene]

# ## add weight to distance_df
# distance_df = read.csv("./data/distance_df.csv", header=TRUE, stringsAsFactors = FALSE)
# distance_weight_df <- distance_df[,paste0(c("eta","epsilon","xc","scaling","xc_epsilon"),"_rela")]
# rownames(distance_weight_df)=distance_df$gene
# 
# # log distance
# 
# gene_weight=(log(ellipse_df$distance_to_O+1)/log(4*ellipse_df$distance_to_center+1+1e-10))^(-1)
# names(gene_weight)=ellipse_df$gene
# # for every gene, multiply the weight to distance_weight_df
# distance_weight_df_weighted=distance_weight_df*gene_weight[rownames(distance_weight_df)]
# distance_df_new=distance_df
# distance_df_new$weight=gene_weight[distance_df_new$gene]
# distance_df_new[,paste0(c("eta","epsilon","xc","scaling","xc_epsilon"),"_rela_new")]=distance_weight_df_weighted
# 
# write.csv(distance_df_new, file="./data/distance_df_new.csv", row.names=FALSE)

