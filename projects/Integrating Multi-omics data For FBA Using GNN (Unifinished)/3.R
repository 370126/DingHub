setwd("C:/Users/86189/Desktop/neo_project")


################ combine through weight mean
w1=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
w2=1-w1

flux_lung=data.frame(reaction=rownames(flux_lung))
for (i in 1:length(w1)){
  comb=as.data.frame(comb_lung$abundence_normalized*w1[i]+comb_lung$Prob*w2[i]);
  rownames(comb)=comb_lung$SYMBOL
  scores<-calculate_reaction_score(comb)
  temp=compute_flux(mras=scores,medium=human_blood)
  flux_lung[[as.character(w1[i])]] = as.vector(temp)
  gc()
}

flux_leukemia=data.frame(reaction=rownames(flux_leukemia))
for (i in 1:length(w1)){
  comb=as.data.frame(comb_leukemia$abundence_normalized*w1[i]+comb_leukemia$Prob*w2[i]);
  rownames(comb)=comb_leukemia$SYMBOL
  scores<-calculate_reaction_score(comb)
  temp=compute_flux(mras=scores,medium=human_blood)
  flux_leukemia[[as.character(w1[i])]] = as.vector(temp)
  gc()
}

rownames(flux_lung)=flux_lung$reaction;flux_lung=flux_lung[,-1]
rownames(flux_leukemia)=flux_leukemia$reaction;flux_leukemia=flux_leukemia[,-1]

## plot of biomass
biomass_lung=flux_lung["biomass_human",];biomass_lung=t(biomass_lung)
biomass_leukemia=flux_leukemia["biomass_human",];biomass_leukemia=t(biomass_leukemia)

ggplot(biomass_lung, aes(x=rownames(biomass_lung),y=biomass_human)) + geom_bar(stat="identity",fill="skyblue2") + ggtitle("Lung") + xlab("Weight of RNA-seq") + ylab("Biomass flux")+theme_minimal()
ggplot(biomass_leukemia, aes(x=rownames(biomass_leukemia),y=biomass_human)) + geom_bar(stat="identity",fill="skyblue2") + ggtitle("Leukemia") + xlab("Weight of RNA-seq") + ylab("Biomass flux")+theme_minimal()

#compute pathway level activity for all samples
pathway_show(human_gem,flux_lung,"lung")

pathway_show(human_gem,flux_leukemia,"leukemia")


################## Bayesian model
#install.packages("rjags")
#install.packages("coda")
library(rjags)  # 用于贝叶斯建模
library(coda)   # 用于MCMC结果分析
model_code <- "
model {
  for (i in 1:N) {
    # 隐变量：z[i] ~ Bernoulli(pi)
    z[i] ~ dbern(pi)
    
    # RNA-seq数据的分布（高斯分布）
    abundence_normalized[i] ~ dnorm(mu[z[i]+1], tau[z[i]+1])
  }
  
  # 超参数的先验分布
  pi ~ dbeta(1, 1)  # 激活比例的先验（均匀分布）
  
  # 激活群的参数
  mu[2] <- 0.8   # 高斯均值
  tau[2] <- 25   # 高斯精度（tau = 1/σ²，σ=0.2）
  
  # 抑制群的参数
  mu[1] <- 0.2
  tau[1] <- 100  # σ=0.1
}
"

# 准备数据
data_list <- list(
  N = nrow(comb_lung),
  abundence_normalized = comb_lung$abundence_normalized
)

# 初始化模型
model <- jags.model(
  textConnection(model_code),
  data = data_list,
  n.chains = 3
)

# 预热（Burn-in）
update(model, 1000)

# 采样
samples <- coda.samples(
  model,
  variable.names = c("z", "pi"),
  n.iter = 5000
)

# 提取隐变量的后验概率
z_post <- summary(samples)$statistics[grep("z", rownames(summary(samples)$statistics)), "Mean"]


z_post=Bayesian(comb_lung)
comb_lung$z_post_prob <- z_post

z_post=Bayesian(comb_leukemia)
comb_leukemia$z_post_prob <- z_post




save(list = ls(), file = "workspace.RData")

load("workspace.RData")
# 查看结果
library(ggplot2)
# 绘制隐变量后验概率 vs RNA-seq数据
ggplot(comb_lung, aes(x = abundence_normalized, y = z_post_prob)) +
  geom_point(size = 1,alpha=0.1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  labs(title = "Posterior Probability of Activation vs RNA-seq Expression",
       x = "Normalized RNA-seq Abundance",
       y = "Posterior Probability of Activation",
       color = "Gene") +
  theme_minimal()



################## Multiply
multi_lung=comb_lung[,c("SYMBOL","abundence_normalized","Prob")]

s1=c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2)
s2=2-s1
flux_multi_lung=data.frame(reaction=rownames(flux_leukemia))
library(METAFlux)
for (i in 1:length(s1)){
  comb=as.data.frame(multi_lung$abundence_normalized^s1[i]*multi_lung$Prob^s2[i]);
  rownames(comb)=comb_lung$SYMBOL;
  scores<-calculate_reaction_score(comb);
  temp=compute_flux(mras=scores,medium=human_blood);
  flux_multi_lung[[as.character(s1[i])]] = as.vector(temp)
  gc()
}

flux_multi_leukemia=data.frame(reaction=rownames(flux_leukemia))
for (i in 1:length(s1)){
  comb=as.data.frame(multi_lung$abundence_normalized^s1[i]*multi_lung$Prob^s2[i]);
  rownames(comb)=comb_lung$SYMBOL;
  scores<-calculate_reaction_score(comb);
  temp=compute_flux(mras=scores,medium=human_blood);
  flux_multi_leukemia[[as.character(s1[i])]] = as.vector(temp)
  gc()
}
