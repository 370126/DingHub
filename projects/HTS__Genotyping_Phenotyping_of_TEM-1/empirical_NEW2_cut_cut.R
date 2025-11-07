library(dplyr)
library(ggplot2)
library(purrr)
library(ggpubr)
library(minpack.lm)
library(gridExtra)


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

# write.csv(fitness_amp_ex_df,"Dose_minilib_summary.csv",row.names = FALSE,quote=FALSE)


# remove rows with CV > 0.3
CV_cutoff=0.35
fitness_amp_ex_df=subset(fitness_amp_ex_df,CV<=CV_cutoff)
wt_km <- kinetics_df %>% filter(Mutant == "WT") %>% pull(P_Km) %>% unique()
wt_vmax <- kinetics_df %>% filter(Mutant == "WT") %>% pull(P_Vmax) %>% unique()
kinetics_df <- kinetics_df %>%
  mutate(
    P_Km_relative = P_Km / wt_km,
    P_Vmax_relative = P_Vmax / wt_vmax,
    log_P_Km_relative = log10(P_Km_relative),
    log_P_Vmax_relative = log10(P_Vmax_relative)
  )
kinetics_df$vmax_km_ratio=kinetics_df$P_Vmax/kinetics_df$P_Km

# rename Mutant "WT" in kinetics_df to "B12"
kinetics_df$Mutant[kinetics_df$Mutant=="WT"]="B12"



# merge fitness_amp_ex_df and kinetics_df by Mutant, but preserve all rows in fitness_amp_ex_df
fitness_amp_km_vmax_df <- merge(fitness_amp_ex_df, kinetics_df, by = "Mutant", all.x = TRUE)
# fitness_amp_km_vmax_df <- merge(fitness_amp_ex_df, kinetics_df, by = "Mutant")
###### NEW: CUTOFF AMP CONCENTRATION
# amp_cutoff=2049
# fitness_amp_km_vmax_df=fitness_amp_km_vmax_df%>%
#   filter(Amp_ex<=amp_cutoff)
EGvar_cutoff=0.05;
fitness_amp_km_vmax_df=fitness_amp_km_vmax_df%>%
  filter(EGvar>=EGvar_cutoff)





cutoff_mutant <- c(
  "A1"=0.65,  "A2"=0.62,  "A3"=0.60,  "A4"=0.60,  "A5"=0.60,  "A6"=0.61,  
  "A7"=0.61,  "A8"=0.65,  "A9"=0.50,  "A10"=NA,  "A11"=0.50,
  "B1"=NA,  "B2"=0.40,  "B3"=0.60,  "B4"=0.62,  "B5"=NA,  "B6"=0.50,
  "B7"=0.40,  "B8"=NA,  "B9"=0.62,  "B10"=0.60,  "B11"=0.60,"B12"=0.70,
  "C1"=0.63,  "C2"=NA,  "C3"=0.35,  "C4"=0.50,  "C5"=0.50,  "C6"=NA,
  "C7"=0.50,  "C8"=0.60,  "C9"=0.50,  "C10"=NA,  "C11"=0.45,
  "D1"=0.10,  "D2"=0.35,  "D3"=NA,  "D4"=0.17,  "D5"=0.50,  "D6"=0.40,
  "D7"=0.30,  "D8"=0.40,  "D9"=NA,  "D11"=0.50
)

cutoff_half_mutant <- c(
  "A1"=0.83,  "A2"=0.80,  "A3"=0.82,  "A4"=0.80,  "A5"=0.85,  "A6"=0.88,  
  "A7"=0.82,  "A8"=0.75,  "A9"=0.70,  "A10"=NA,  "A11"=0.80,
  "B1"=NA,  "B2"=0.60,  "B3"=0.80,  "B4"=0.71,  "B5"=NA,  "B6"=0.80,
  "B7"=0.60,  "B8"=NA,  "B9"=0.80,  "B10"=0.75,  "B11"=0.80,"B12"=0.75,
  "C1"=0.80,  "C2"=NA,  "C3"=0.50,  "C4"=0.70,  "C5"=0.65,  "C6"=NA,
  "C7"=0.75,  "C8"=0.70,  "C9"=0.70,  "C10"=NA,  "C11"=0.70,
  "D1"=0.50,  "D2"=0.70,  "D3"=NA,  "D4"=0.40,  "D5"=0.65,  "D6"=0.50,
  "D7"=0.40,  "D8"=0.70,  "D9"=NA,  "D11"=0.75
)

# add cutoff column to fitness_amp_km_vmax_df
fitness_amp_km_vmax_df$cutoff <- cutoff_mutant[fitness_amp_km_vmax_df$Mutant]
fitness_amp_km_vmax_df$cutoff_half <- cutoff_half_mutant[fitness_amp_km_vmax_df$Mutant]






# df=fitness_amp_km_vmax_df[,c("Mutant","Amp_ex","EGvar","cutoff","cutoff_half")]
# write.csv(df,"summary_phaseI_phaseIA_cutoff.csv",row.names = FALSE,quote=FALSE)
# 

fitness_amp_km_vmax_df=na.omit(fitness_amp_km_vmax_df)

# filter out rows with EGvar < cutoff
fitness_amp_km_vmax_cut_cut_df <- subset(fitness_amp_km_vmax_df, EGvar >= cutoff_half)


df=fitness_amp_km_vmax_cut_cut_df[,c("Mutant","Amp_ex","EGvar","cutoff","cutoff_half")]



# plot dose response curve for each mutant after cutoff, with error bars, with automatic x-axis and y-axis
p = ggplot(fitness_amp_km_vmax_cut_cut_df, aes(x = Amp_ex, y = EGvar, color = Mutant)) +
  geom_point() +
  geom_errorbar(aes(ymin = EGvar - sd, ymax = EGvar + sd), width = 0.1) +
  geom_line() +
  # add horizontal line at cutoff, for each mutant
  # geom_hline(aes(yintercept = cutoff), linetype = "dashed", color = "red") +
  # annotate("text", x = 10, y = cutoff+0.02, label = "phase I cutoff", color = "red",size=3) +
  # plot every mutant separately
  # facet_wrap(~ Mutant, ncol = 5, scales = "free") + # free scale
  # scale_x_log10() +
  # annotation_logticks(sides = "b") +
  # coord_cartesian(ylim = c(0, 1.1)) +
  labs(title = "Dose Response Curve after Cutoff", x = "Ampicillin Concentration (Amp_ex)", y = "EGvar") +
  theme_classic()+
  # remove legend
  theme(legend.position = "none")
print(p)
# # ggsave("figures/dose_response_curve_nolog_separate_phaseI_half.png", p, width = 40, height = 60, units = "cm", dpi = 500)


dose_para_cut_cut_df=data.frame(
  Mutant=unique(fitness_amp_km_vmax_cut_cut_df$Mutant),
  E=NA,k=NA,A=dose_para_cutoff_df$A,
  p_E=NA,p_k=NA,p_A=NA,
  r_squared=NA
)


dose_response_half_model <- function(x, E, k, A=0.5) {
  A+(1-A)/(1+(x/E)^k)
}

for (i in 1:nrow(dose_para_cut_cut_df)){
  mutant=dose_para_cut_cut_df$Mutant[i]
  A_temp=dose_para_cut_cut_df$A[i]
  df=subset(fitness_amp_km_vmax_cut_cut_df,Mutant==mutant)
  fit <- tryCatch({
    minpack.lm::nlsLM(
      EGvar ~ dose_response_half_model(Amp_ex, E, k, A=A_temp),
      data = df,
      start = list(E = median(df$EGvar), k = 4),
      lower = c(1e-6, 0.01),
      upper = c(100000, 100000)
    )
  }, error = function(e) NULL)
  if (!is.null(fit)){
    dose_para_cut_cut_df$E[i]=coef(fit)["E"]
    dose_para_cut_cut_df$k[i]=coef(fit)["k"]
    # dose_para_cut_cut_df$A[i]=coef(fit)["A"]
    dose_para_cut_cut_df$p_E[i]=summary(fit)$coefficients["E","Pr(>|t|)"]
    dose_para_cut_cut_df$p_k[i]=summary(fit)$coefficients["k","Pr(>|t|)"]
    # dose_para_cut_cut_df$p_A[i]=summary(fit)$coefficients["A","Pr(>|t|)"]
    # calculate R_squared
    ss_total <- sum((df$EGvar - mean(df$EGvar))^2)
    ss_residual <- sum(residuals(fit)^2)
    r_squared <- 1 - (ss_residual / ss_total)
    dose_para_cut_cut_df$r_squared[i]=r_squared}
}
# visualize the all the fitting results
plot_dose_cutoff <- function(mutant){
  df=subset(fitness_amp_km_vmax_cut_cut_df,Mutant==mutant)
  E=dose_para_cut_cut_df$E[dose_para_cut_cut_df$Mutant==mutant]
  k=dose_para_cut_cut_df$k[dose_para_cut_cut_df$Mutant==mutant]
  A=dose_para_cut_cut_df$A[dose_para_cut_cut_df$Mutant==mutant]
  p_E=dose_para_cut_cut_df$p_E[dose_para_cut_cut_df$Mutant==mutant]
  p_k=dose_para_cut_cut_df$p_k[dose_para_cut_cut_df$Mutant==mutant]
  # p_A=dose_para_cut_cut_df$p_A[dose_para_cut_cut_df$Mutant==mutant]
  r_squared=dose_para_cut_cut_df$r_squared[dose_para_cut_cut_df$Mutant==mutant]
  ggplot(df,aes(x=Amp_ex,y=EGvar))+
    geom_point()+
    geom_errorbar(aes(ymin = EGvar - sd, ymax = EGvar + sd), width = 0.1) +
    stat_function(fun = dose_response_half_model,args = list(E=E,k=k,A=A),color="red")+
    # scale_x_log10()+
    # annotation_logticks(sides = "b")+
    # scale_x_continues(limits = c(0,1000))+
    coord_cartesian(ylim = c(0,1.1))+
    labs(title=(paste0(mutant,
                       "\nE=",round(E,3),ifelse(p_E<0.05,"*",""),", k=",round(k,3),ifelse(p_k<0.05,"*",""),
                       "\nR²=",round(r_squared,4))),
         x="log(Amp_ex)",y="EGvar")+
    # make title smaller
    theme_classic2()+
    theme(plot.title = element_text(size = 12))
}
#plot_list_cutoff=map(dose_para_cut_cut_df$Mutant,plot_dose_cutoff)
#p=do.call(grid.arrange,c(plot_list_cutoff,ncol=5))
#ggsave("figures/dose_response_curve_fitting_nolog_phaseI_half.png", p, width = 40, height = 60, units = "cm", dpi = 500)


## calculate ECs
dose_para_cut_cut_df_filtered=dose_para_cut_cut_df%>%
  filter(p_E<0.05 & p_k<0.05 & r_squared>0.6)
# merge with kinetics_df
dose_para_cut_cut_df_filtered <- merge(dose_para_cut_cut_df_filtered, kinetics_df, by = "Mutant")
dose_para_cut_cut_df_filtered$vmax_km_ratio=dose_para_cut_cut_df_filtered$P_Vmax/dose_para_cut_cut_df_filtered$P_Km
dose_para_cut_cut_df_filtered$vmax_km_ratio_relative=dose_para_cut_cut_df_filtered$P_Vmax_relative/dose_para_cut_cut_df_filtered$P_Km_relative

# filter out bad mutants
bad_mutant <- c("A5", "C10", "D9", "D10", "D11", "B1", "D1","C11","D6")
dose_para_cut_cut_df_filtered=dose_para_cut_cut_df_filtered%>%
  filter(!Mutant %in% bad_mutant)

# filter outiliers
# dose_para_cut_cut_df_filtered=dose_para_cut_cut_df_filtered%>%
#   filter(P_Km<500 & P_Vmax>0.1 & P_Vmax<5)


# calculate EC values
dose_para_cut_cut_df_filtered$EC01=0.01^(1/dose_para_cut_cut_df_filtered$k)*dose_para_cut_cut_df_filtered$E
dose_para_cut_cut_df_filtered$EC05=0.05^(1/dose_para_cut_cut_df_filtered$k)*dose_para_cut_cut_df_filtered$E
dose_para_cut_cut_df_filtered$EC10=0.1^(1/dose_para_cut_cut_df_filtered$k)*dose_para_cut_cut_df_filtered$E
dose_para_cut_cut_df_filtered$EC15=0.15^(1/dose_para_cut_cut_df_filtered$k)*dose_para_cut_cut_df_filtered$E
dose_para_cut_cut_df_filtered$EC90=0.9^(1/dose_para_cut_cut_df_filtered$k)*dose_para_cut_cut_df_filtered$E
dose_para_cut_cut_df_filtered$EC99=0.99^(1/dose_para_cut_cut_df_filtered$k)*dose_para_cut_cut_df_filtered$E

# write.csv(dose_para_cut_cut_df_filtered,"summary_phaseI_half_fitting_kinetics_filtered.csv",row.names = FALSE,quote=FALSE)



## plot E/k vs Km,Vmax,Vmax/Km separately but on the same graph
shape_list=c("E","k")
para_list=c("P_Km_relative","P_Vmax_relative","vmax_km_ratio_relative")
plot_enzyme <- function(shape,para){
  ggplot(dose_para_cut_cut_df_filtered, aes_string(x = para, y = shape)) +
    geom_point() +
    geom_text(aes(label=Mutant),hjust=0, vjust=0, size=3)+
    geom_smooth(method = "loess", se = TRUE, color = "red") +
    labs(title = paste0(shape, " vs ", para), x = para, y = shape) +
    theme_classic()
}
plot_list_shape=list()
for (shape in shape_list){
  for (para in para_list){
    plot_list_shape[[paste0(shape,"_",para)]]=plot_enzyme(shape,para)
  }
}
p=do.call(grid.arrange,c(plot_list_shape,ncol=3))
ggsave("figures/E_k_vs_Km_Vmax_phaseI_half.png", p, width = 30, height = 20, units = "cm", dpi = 500)


EC_list=c("EC01","EC05","EC10","EC15","EC90","EC99")
plot_list_EC=list()
for (EC in EC_list){
  for (para in para_list){
    plot_list_EC[[paste0(EC,"_",para)]]=plot_enzyme(EC,para)
  }
}
p=do.call(grid.arrange,c(plot_list_EC,ncol=6))
ggsave("figures/EC_vs_Km_Vmax_phaseI_half.png", p, width = 80, height = 40, units = "cm", dpi = 500)




bad_mutant2="B2" # remove B2 because it is an outlier
dose_para_cut_cut_df_filtered=dose_para_cut_cut_df_filtered%>%
  filter(Mutant!=bad_mutant2)
## fitting EC vs Vmax
EC_lm_results <- data.frame(EC = EC_list, slope = NA, intercept = NA, r_squared = NA, p_value = NA, rho_pearson=NA, rho_speaerman=NA)
for (i in 1:nrow(EC_lm_results)){
  EC <- EC_lm_results$EC[i]
  df <- dose_para_cut_cut_df_filtered[, c(EC, "P_Vmax_relative")]
  colnames(df) <- c("EC", "P_Vmax_relative")
  fit <- lm(EC ~ P_Vmax_relative, data = df)
  summary_fit <- summary(fit)
  EC_lm_results$slope[i] <- summary_fit$coefficients["P_Vmax_relative", "Estimate"]
  EC_lm_results$intercept[i] <- summary_fit$coefficients["(Intercept)", "Estimate"]
  EC_lm_results$r_squared[i] <- summary_fit$r.squared
  EC_lm_results$p_value[i] <- summary_fit$coefficients["P_Vmax_relative", "Pr(>|t|)"]
  # calculate pearson correlation
  cor_pearson <- cor.test(df$EC, df$P_Vmax_relative, method = "pearson")
  EC_lm_results$rho_pearson[i] <- cor_pearson$estimate
  # calculate spearman correlation
  cor_spearman <- cor.test(df$EC, df$P_Vmax_relative, method = "spearman")
  EC_lm_results$rho_speaerman[i] <- cor_spearman$estimate
}

# plot all EC fittings(EC v.s. vmax) on same graph, separately with r_squared and rhos
plot_list_EC_lm <- list()
for (i in 1:nrow(EC_lm_results)){
  EC <- EC_lm_results$EC[i]
  slope <- round(EC_lm_results$slope[i], 3)
  intercept <- round(EC_lm_results$intercept[i], 3)
  r_squared <- round(EC_lm_results$r_squared[i], 3)
  p_value <- signif(EC_lm_results$p_value[i], 3)
  rho_pearson <- round(EC_lm_results$rho_pearson[i], 3)
  rho_speaerman <- round(EC_lm_results$rho_speaerman[i], 3)
  
  plot_list_EC_lm[[EC]] <- ggplot(dose_para_cut_cut_df_filtered, aes_string(x = "P_Vmax_relative", y = EC)) +
    geom_point() +
    geom_text(aes(label=Mutant),hjust=0, vjust=0, size=3)+
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(title = paste0(EC,
                        "\nSlope=", slope, ", Intercept=", intercept,
                        ", R²=", r_squared, ifelse(p_value<0.05,"*",""),
                        "\nρ_pearson=", rho_pearson, ifelse(cor_pearson$p.value<0.05,"*",""),
                        " | ρ_spearman=", rho_speaerman, ifelse(cor_spearman$p.value<0.05,"*","")
    ),
    x = "P_Vmax_relative", y = EC) +
    theme_classic()+
    theme(plot.title = element_text(size = 16))
}
p=do.call(gridExtra::grid.arrange,c(plot_list_EC_lm,ncol=3))
ggsave("figures/EC_vs_Vmax_phaseI_half.png", p, width = 60, height = 40, units = "cm", dpi = 500)







