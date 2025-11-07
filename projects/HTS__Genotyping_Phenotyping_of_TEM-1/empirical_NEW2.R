library(dplyr)
library(ggplot2)
library(purrr)
library(ggpubr)
library(minpack.lm)


kinetics_df <- read.csv("./35C_protein_Km_Vmax.tsv", sep = "\t")
# fitness_amp_ex_df <- read.csv("./fitness_amp_ex_df.csv")
fitness_amp_ex_df=read.csv("Dose_response_summary.txt",sep="\t")
colnames(fitness_amp_ex_df)[which(colnames(fitness_amp_ex_df)=="amp_conc")]="Amp_ex"
colnames(fitness_amp_ex_df)[which(colnames(fitness_amp_ex_df)=="mutant")]="Mutant"
# calculate the mean and sd of rep1~3
colnames(fitness_amp_ex_df)
# [1] "Mutant"   "Amp_ex"   "rep1"     "rep2"     "rep3"
# calculate mean and sd of rep1~3
fitness_amp_ex_df <- fitness_amp_ex_df %>%
  rowwise() %>%
  mutate(
    mean = mean(c(rep1, rep2, rep3), na.rm = TRUE),
    sd = sd(c(rep1, rep2, rep3), na.rm = TRUE)
  ) %>%
  ungroup()
fitness_amp_ex_df$CV=fitness_amp_ex_df$sd/fitness_amp_ex_df$mean
# remove rows with CV > 0.3
CV_cutoff=0.3
fitness_amp_ex_df=subset(fitness_amp_ex_df,CV<=CV_cutoff)
# calculate GVR, EGVR
fitness_amp_ex_df$Gvar=log2((fitness_amp_ex_df$mean-0.05)*2/0.02)
# mutate NA and minus infinite to 0
fitness_amp_ex_df$Gvar[is.infinite(fitness_amp_ex_df$Gvar) | is.na(fitness_amp_ex_df$Gvar) | fitness_amp_ex_df$Gvar<0]=0
# calculate EGvar=GVar/GVar(amp_ex=0) for each mutant
fitness_amp_ex_df <- fitness_amp_ex_df %>%
  group_by(Mutant) %>%
  # divide Gvar by the Gvar when Amp_ex=0
  mutate(EGvar = Gvar / Gvar[Amp_ex == 0]) %>%
  ungroup()



wt_km <- kinetics_df %>% filter(Mutant == "WT") %>% pull(P_Km) %>% unique()
wt_vmax <- kinetics_df %>% filter(Mutant == "WT") %>% pull(P_Vmax) %>% unique()
kinetics_df <- kinetics_df %>%
  mutate(
    P_Km_relative = P_Km / wt_km,
    P_Vmax_relative = P_Vmax / wt_vmax,
    log_P_Km_relative = log10(P_Km_relative),
    log_P_Vmax_relative = log10(P_Vmax_relative)
  )


fitness_amp_km_vmax_df <- merge(fitness_amp_ex_df, kinetics_df, by = "Mutant")



dose_response_model <- function(x, E, k, A) {
  B <- 1.0  # 上限固定
  A + (B - A) * exp(-(x / E)^k)
}
dose_para_df=data.frame(
  Mutant=unique(fitness_amp_km_vmax_df$Mutant),
  E=NA,k=NA,A=NA,
  p_E=NA,p_k=NA,p_A=NA,
  r_squared=NA
)
for (i in 1:nrow(dose_para_df)){
  mutant=dose_para_df$Mutant[i]
  df=subset(fitness_amp_km_vmax_df,Mutant==mutant)
  fit <- tryCatch({
    minpack.lm::nlsLM(
      EGvar ~ dose_response_model(Amp_ex, E, k, A),
      data = df,
      start = list(E = median(df$Amp_ex), k = 1,A=0),
      lower = c(1e-6, 0.01, 0),
      # upper = c(10000, 1000)
    )
  }, error = function(e) NULL)
  if (!is.null(fit)){
    dose_para_df$E[i]=coef(fit)["E"]
    dose_para_df$k[i]=coef(fit)["k"]
    dose_para_df$A[i]=coef(fit)["A"]
    dose_para_df$p_E[i]=summary(fit)$coefficients["E","Pr(>|t|)"]
    dose_para_df$p_k[i]=summary(fit)$coefficients["k","Pr(>|t|)"]
    dose_para_df$p_A[i]=summary(fit)$coefficients["A","Pr(>|t|)"]
    # calculate R_squared
    ss_total <- sum((df$EGvar - mean(df$EGvar))^2)
    ss_residual <- sum(residuals(fit)^2)
    r_squared <- 1 - (ss_residual / ss_total)
    dose_para_df$r_squared[i]=r_squared}
}
# visualize the all the fitting results
plot_dose <- function(mutant){
  df=subset(fitness_amp_km_vmax_df,Mutant==mutant)
  E=dose_para_df$E[dose_para_df$Mutant==mutant]
  k=dose_para_df$k[dose_para_df$Mutant==mutant]
  A=dose_para_df$A[dose_para_df$Mutant==mutant]
  p_E=dose_para_df$p_E[dose_para_df$Mutant==mutant]
  p_k=dose_para_df$p_k[dose_para_df$Mutant==mutant]
  p_A=dose_para_df$p_A[dose_para_df$Mutant==mutant]
  r_squared=dose_para_df$r_squared[dose_para_df$Mutant==mutant]
  ggplot(df,aes(x=Amp_ex,y=EGvar))+
    geom_point()+
    stat_function(fun = dose_response_model,args = list(E=E,k=k,A=A),color="red")+
    scale_x_log10()+
    annotation_logticks(sides = "b")+
    # scale_x_continues(limits = c(0,1000))+
    ylim(0,1.05)+
    labs(title=(paste0(mutant,
                       "\nE=",round(E,2),ifelse(p_E<0.05,"*",""),", k=",round(k,2),ifelse(p_k<0.05,"*",""),", A=",round(A,2),ifelse(p_A<0.05,"*",""),
                       "\nR²=",round(r_squared,2))),
         x="log(Amp_ex)",y="EGvar")+
    theme_classic2()
}
plot_list=map(dose_para_df$Mutant,plot_dose)
library(gridExtra)
# p=do.call(grid.arrange,c(plot_list,ncol=5))
# ggsave("figures/dose_response_fitting_exp.png",p,width=50,height=60,units="cm",dpi=500)
# filter out p_E>0.05 or p_k>0.05 or r_squared<0.6
print("Before filtering:")
print(nrow(dose_para_df))
dose_para_df=dose_para_df%>%
  filter(p_E<=0.05 & p_k<=0.05 & r_squared>=0.7)
print("After filtering:")
print(nrow(dose_para_df))

#### E,k and Km,Vmax correlation
dose_para_df$EC99=dose_para_df$E*(99)^(1/dose_para_df$k)
dose_para_df$EC90=dose_para_df$E*(90)^(1/dose_para_df$k)

dose_para_df$EC15=dose_para_df$E*(1/0.85-1)^(1/dose_para_df$k)
dose_para_df$EC10=dose_para_df$E*(1/0.9-1)^(1/dose_para_df$k)
dose_para_df$EC05=dose_para_df$E*(1/0.95-1)^(1/dose_para_df$k)
dose_para_df$EC01=dose_para_df$E*(1/0.99-1)^(1/dose_para_df$k)



## plot E v.s. k
p=ggplot(dose_para_df,aes(x=k,y=E))+
  geom_point(size=2,alpha=0.5)+
  # add mutant name besides points
  geom_text(aes(label=Mutant),hjust=-0.5, vjust=-0.5,size=2)+
  scale_x_log10()+
  scale_y_log10()+
  # geom_smooth(method="loess",se=TRUE,color="red")+
  labs(title="E vs k",
       x="k",y="E")+
  theme_classic2();p








# filter out bad mutants
bad_mutant <- c("A5", "C10", "D9", "D10", "D11", "B1", "D1","C11")
dose_para_df=dose_para_df%>%
  filter(!Mutant %in% bad_mutant)


#### E,k and Km,Vmax correlation
dose_para_df=merge(dose_para_df,kinetics_df,by="Mutant")
dose_para_df$vmax_km_ratio=dose_para_df$P_Vmax/dose_para_df$P_Km
dose_para_df$vmax_km_ratio_relative=dose_para_df$P_Vmax_relative/dose_para_df$P_Km_relative

# # filter outiliers
dose_para_df=dose_para_df%>%
  filter(P_Km<500 & P_Vmax<5 & P_Vmax>0.1)

cor.test(dose_para_df$E,dose_para_df$P_Km,method="spearman")
cor.test(dose_para_df$k,dose_para_df$P_Vmax,method="spearman")
cor.test(dose_para_df$E,dose_para_df$vmax_km_ratio,method="spearman")
cor.test(dose_para_df$k,dose_para_df$vmax_km_ratio,method="spearman")






## use relative Km and Vmax
# plot E vs Km
p1=ggplot(dose_para_df,aes(x=P_Km_relative,y=E))+
  geom_point(size=2,alpha=0.5)+
  # add mutant name besides points
  geom_text(aes(label=Mutant),hjust=-0.5, vjust=-0.5,size=2)+
  geom_smooth(method="loess",se=TRUE,color="red")+
  labs(title="E vs Km",
       x="Relative Km",y="E")+
  theme_classic2();p1
# plot E vs Vmax
p2=ggplot(dose_para_df,aes(x=P_Vmax_relative,y=E))+
  geom_point(size=2,alpha=0.5)+
  geom_text(aes(label=Mutant),hjust=-0.5, vjust=-0.5,size=2)+
  geom_smooth(method="loess",se=TRUE,color="red")+
  labs(title="E vs Vmax",
       x="Relative Vmax",y="E")+
  theme_classic2();p2
# plot k vs Km
p3=ggplot(dose_para_df,aes(x=P_Km_relative,y=k))+
  geom_point(size=2,alpha=0.5)+
  geom_text(aes(label=Mutant),hjust=-0.5, vjust=-0.5,size=2)+
  geom_smooth(method="loess",se=TRUE,color="red")+
  labs(title="k vs Km",
       x="Relative Km",y="k")+
  theme_classic2();p3
# plot k vs Vmax
p4=ggplot(dose_para_df,aes(x=P_Vmax_relative,y=k))+
  geom_point(size=2,alpha=0.5)+
  geom_text(aes(label=Mutant),hjust=-0.5, vjust=-0.5,size=2)+
  geom_smooth(method="lm",se=TRUE,color="red")+
  labs(title="k vs Vmax",
       x="Relative Vmax",y="k")+
  theme_classic2();p4
# plot E vs Vmax/Km
p5=ggplot(dose_para_df,aes(x=vmax_km_ratio_relative,y=E))+
  geom_point(size=2,alpha=0.5)+
  geom_text(aes(label=Mutant),hjust=-0.5, vjust=-0.5,size=2)+
  geom_smooth(method="loess",se=TRUE,color="red")+
  labs(title="E vs Vmax/Km",
       x="Relative Vmax/Km",y="E")+
  theme_classic2();p5
# plot k vs Vmax/Km
p6=ggplot(dose_para_df,aes(x=vmax_km_ratio_relative,y=k))+
  geom_point(size=2,alpha=0.5)+
  geom_text(aes(label=Mutant),hjust=-0.5, vjust=-0.5,size=2)+
  geom_smooth(method="loess",se=TRUE,color="red")+
  labs(title="k vs Vmax/Km",
       x="Relative Vmax/Km",y="k")+
  theme_classic2();p6
library(gridExtra)
#p=grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3)
#ggsave("figures/E_k_vs_Km_Vmax_exp_filtered.png",p,width=50,height=30,units="cm",dpi=500)



## plot EC01~15 vs Km,Vmax,Vmax/Km separately but on the same graph
EC_list=c("EC01","EC05","EC10","EC15","EC90","EC99")
para_list=c("P_Km_relative","P_Vmax_relative","vmax_km_ratio_relative")
plot_EC <- function(EC,para){
  p=ggplot(dose_para_df,aes_string(x=para,y=EC))+
    geom_point(size=2,alpha=0.5)+
    geom_text(aes(label=Mutant),hjust=-0.5, vjust=-0.5,size=2)+
    geom_smooth(method="loess",se=TRUE,color="red")+
    labs(title=paste0(EC," vs ",para),
         x=paste0("Relative ",para),y=EC)+
    theme_classic2()
  return(p)
}
p_list=list()
for (EC in EC_list){
  for (para in para_list){
    p_list[[paste0(EC,"_",para)]]=plot_EC(EC,para)
  }
}
p_all=do.call(grid.arrange,c(p_list,ncol=6))
ggsave("figures/EC_vs_Km_Vmax_exp_filtered.png",p_all,width=100,height=40,units="cm",dpi=500)








## fit E v.s.vmax using log model
E_vmax_x <- function(P_Vmax, b0, b1, b2) {
  b0*log2(P_Vmax+b2)+b1
}
fit_vmax <- tryCatch({
  minpack.lm::nlsLM(
    E ~ E_vmax_x(P_Vmax_relative, b0, b1, b2),
    data = dose_para_df,
    start = list(b0 = 0.3, b1 = 2, b2 = 0.05),
    # lower = c(0, 0, 0),
    # upper = c(1, 10, 1),
    control = minpack.lm::nls.lm.control(maxiter = 1024)
  )
}, error = function(e) NULL)
print(fit_vmax)
if (!is.null(fit_vmax)){
  b0_fit=coef(fit_vmax)["b0"]
  b1_fit=coef(fit_vmax)["b1"]
  b2_fit=coef(fit_vmax)["b2"]
  fit_vmax_df=data.frame(
    x=seq(min(dose_para_df$P_Vmax_relative),max(dose_para_df$P_Vmax_relative),length.out=100),
    y=E_vmax_x(seq(min(dose_para_df$P_Vmax_relative),max(dose_para_df$P_Vmax_relative),length.out=100),b0_fit,b1_fit,b2_fit)
  )
  # calculate R_squared
  ss_total <- sum((dose_para_df$E - mean(dose_para_df$E))^2)
  ss_residual <- sum(residuals(fit_vmax)^2)
  r_squared_vmax <- 1 - (ss_residual / ss_total)
  p_vmax=ggplot(dose_para_df,aes(x=P_Vmax_relative,y=E))+
    geom_point(size=2,alpha=0.5)+
    geom_text(aes(label=Mutant),hjust=-0.5, vjust=-0.5,size=2)+
    geom_line(data=fit_vmax_df,aes(x=x,y=y),color="blue")+
    # add text r_squared
    annotate("text", x = max(dose_para_df$P_Vmax_relative)/2, y = max(dose_para_df$E)*0.9,
             label = paste0("R² = ", round(r_squared_vmax, 3)), size = 6, color = "black")+
    # xlim(0,0.12)+
    labs(title="Fitted equation: E = b0*log2(Vmax + b2) + b1",
         subtitle=paste0("b0=",round(b0_fit,2),", b1=",round(b1_fit,2),", b2=",round(b2_fit,2)),
         x="Relative Vmax",y="E")+
    theme_classic2()
  print(p_vmax)
  # ggsave("figures/E_vs_Vmax_fitting_log.png",width=15,height=10,units="cm",dpi=500)
}



## fit E v.s. Vmax/Km using logarithmic curve
E_km_vmax_x <- function(P_Km, P_Vmax, a0, a1, a2) {
  # seem like a logarithmic curve
  a0*log2((P_Vmax/P_Km)+a2)+a1
}
fit1 <- tryCatch({
  minpack.lm::nlsLM(
    E ~ E_km_vmax_x(P_Km_relative, P_Vmax_relative, a0, a1, a2),
    data = dose_para_df,
    start = list(a0 = 0.3, a1 = 2, a2 = 0.05),
    # lower = c(0, 0, 0),
    # upper = c(1, 10, 1),
    control = minpack.lm::nls.lm.control(maxiter = 1024)
  )
}, error = function(e) NULL)
print(fit1)
if (!is.null(fit1)){
  summary(fit1)
  a0_fit=coef(fit1)["a0"]
  a1_fit=coef(fit1)["a1"]
  a2_fit=coef(fit1)["a2"]
  fit1_df=data.frame(
    x=seq(min(dose_para_df$vmax_km_ratio_relative),max(dose_para_df$vmax_km_ratio_relative),length.out=100),
    y=E_km_vmax_x(1,seq(min(dose_para_df$vmax_km_ratio_relative),max(dose_para_df$vmax_km_ratio_relative),length.out=100),a0_fit,a1_fit,a2_fit)
  )
  # calculate R_squared
  ss_total <- sum((dose_para_df$E - mean(dose_para_df$E))^2)
  ss_residual <- sum(residuals(fit1)^2)
  r_squared1 <- 1 - (ss_residual / ss_total)
  p1=ggplot(dose_para_df,aes(x=vmax_km_ratio_relative,y=E))+
    geom_point(size=2,alpha=0.5)+
    geom_text(aes(label=Mutant),hjust=-0.5, vjust=-0.5,size=2)+
    geom_line(data=fit1_df,aes(x=x,y=y),color="blue")+
    # add text r_squared
    annotate("text", x = max(dose_para_df$vmax_km_ratio_relative)/2, y = max(dose_para_df$E)*0.9,
             label = paste0("R² = ", round(r_squared1, 3)), size = 6, color = "black")+
    # xlim(0,0.12)+
    labs(title="Fitted equation: E = a0*log2(Vmax/Km + a2) + a1",
         subtitle=paste0("a0=",round(a0_fit,2),", a1=",round(a1_fit,2),", a2=",round(a2_fit,2)),
         x="Relative Vmax/Km",y="E")+
    theme_classic2()
  print(p1)
  # ggsave("figures/E_vs_Vmax_Km_fitting_log.png",width=15,height=10,units="cm",dpi=500)
}

# fit E v.s. Vmax/Km using inverse function
E_km_vmax_inv_x <- function(P_Km, P_Vmax, a3, a4, a5) {
  a4 + a5 / (a3 + (P_Vmax / P_Km))
}
fit2 <- tryCatch({
  minpack.lm::nlsLM(
    E ~ E_km_vmax_inv_x(P_Km_relative, P_Vmax_relative, a3, a4, a5),
    data = dose_para_df,
    start = list(a3 = 0.01, a4 = 0.3, a5 = 0.5),
    lower = c(-Inf, -Inf, 0),
    # upper = c(1, 1, 10),
    control = minpack.lm::nls.lm.control(maxiter = 1024)
  )
}, error = function(e) NULL)
print(fit2)
if (!is.null(fit2)){
  summary(fit2)
  a3_fit=coef(fit2)["a3"]
  a4_fit=coef(fit2)["a4"]
  a5_fit=coef(fit2)["a5"]
  fit2_df=data.frame(
    x=seq(min(dose_para_df$vmax_km_ratio_relative),max(dose_para_df$vmax_km_ratio_relative),length.out=100),
    y=E_km_vmax_inv_x(1,seq(min(dose_para_df$vmax_km_ratio_relative),max(dose_para_df$vmax_km_ratio_relative),length.out=100),a3_fit,a4_fit,a5_fit)
  )
  # calculate R_squared
  ss_total <- sum((dose_para_df$E - mean(dose_para_df$E))^2)
  ss_residual <- sum(residuals(fit2)^2)
  r_squared2 <- 1 - (ss_residual / ss_total)
  p2=ggplot(dose_para_df,aes(x=vmax_km_ratio_relative,y=E))+
    geom_point(size=2,alpha=0.5)+
    geom_text(aes(label=Mutant),hjust=-0.5, vjust=-0.5,size=2)+
    geom_line(data=fit2_df,aes(x=x,y=y),color="purple")+
    # add text r_squared
    annotate("text", x = max(dose_para_df$vmax_km_ratio_relative)/2, y = max(dose_para_df$E)*0.9,
             label = paste0("R² = ", round(r_squared2, 3)), size = 6, color = "black")+
    # xlim(0,0.12)+
    labs(title="Fitted equation: E = a4 + a5/(a3 + Vmax/Km)",
         subtitle=paste0("a3=",round(a3_fit,2),", a4=",round(a4_fit,2),", a5=",round(a5_fit,2)),
         x="Relative Vmax/Km",y="E")+
    theme_classic2()
  print(p2)
  # ggsave("figures/E_vs_Vmax_Km_fitting_inverse.png",width=15,height=10,units="cm",dpi=500)
}

# fit E v.s. vmax/km using sigmoid function
E_km_vmax_sig_x <- function(P_Km_relative , P_Vmax_relative, a6, a7, a8) {
  a6/(1+exp(-a7*(P_Vmax_relative/P_Km_relative - a8)))
}
fit3 <- tryCatch({
  minpack.lm::nlsLM(
    E ~ E_km_vmax_sig_x(P_Km_relative, P_Vmax_relative, a6, a7, a8),
    data = dose_para_df,
    start = list(a6 = 2000, a7 = 80, a8= 0.05),
    # lower = c(0, -10000),
    # upper = c(10000, 10000),
    control = minpack.lm::nls.lm.control(maxiter = 1024)
  )
}, error = function(e) NULL)
print(fit3)
if (!is.null(fit3)){
  summary(fit3)
  a6_fit=coef(fit3)["a6"]
  a7_fit=coef(fit3)["a7"]
  a8_fit=coef(fit3)["a8"]
  fit3_df=data.frame(
    x=seq(min(dose_para_df$vmax_km_ratio_relative),max(dose_para_df$vmax_km_ratio_relative),length.out=100),
    y=E_km_vmax_sig_x(1,seq(min(dose_para_df$vmax_km_ratio_relative),max(dose_para_df$vmax_km_ratio_relative),length.out=100),a6_fit,a7_fit,a8_fit)
  )
  # calculate R_squared
  ss_total <- sum((dose_para_df$E - mean(dose_para_df$E))^2)
  ss_residual <- sum(residuals(fit3)^2)
  r_squared3 <- 1 - (ss_residual / ss_total)
  p3=ggplot(dose_para_df,aes(x=vmax_km_ratio_relative,y=E))+
    geom_point(size=2,alpha=0.5)+
    geom_text(aes(label=Mutant),hjust=-0.5, vjust=-0.5,size=2)+
    geom_line(data=fit3_df,aes(x=x,y=y),color="green")+
    # add text r_squared
    annotate("text", x = max(dose_para_df$vmax_km_ratio_relative)/2, y = max(dose_para_df$E)*0.9,
             label = paste0("R² = ", round(r_squared3, 3)), size = 6, color = "black")+
    # xlim(0,0.12)+
    labs(title="Fitted equation: E = a6/(1+exp(-a7*(Vmax/Km - a8)))",
         subtitle=paste0("a6=",round(a6_fit,2),", a7=",round(a7_fit,2),", a8=",round(a8_fit,2)),
         x="Relative Vmax/Km",y="E")+
    theme_classic2()
  print(p3)
  # ggsave("figures/E_vs_Vmax_Km_fitting_sigmoid.png",width=15,height=10,units="cm",dpi=500)
}
# plot all three fitting results
p_all=grid.arrange(p1,p2,p3,ncol=3)
ggsave("figures/E_vs_Vmax_Km_fitting_exp.png",p_all,width=45,height=15,units="cm",dpi=500)





