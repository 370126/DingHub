library(data.table)
library(ggplot2)
library(tidyverse)
library(pheatmap)

###################### read data ####################
# use data.table to read in big data
dose_all_df <- fread("Amp_2_2048_median_EGvar.stat")
dose_all_df <- as.data.frame(dose_all_df)
colnames(dose_all_df) <- c("site","Amp_2","Amp_4","Amp_8","Amp_16","Amp_32","Amp_64","Amp_128","Amp_256","Amp_512","Amp_1024","Amp_2048")

dose_all_0_df <- fread("Amp_0_median_EGvar.stat")
dose_all_0_df <- as.data.frame(dose_all_0_df)
colnames(dose_all_0_df) <- c("site","Amp_0")
dose_all_df=merge(dose_all_df,dose_all_0_df,by="site",all.x=TRUE)
rm(dose_all_0_df)
# reorder columns
dose_all_df <- dose_all_df[, c("site","Amp_0","Amp_2","Amp_4","Amp_8","Amp_16","Amp_32","Amp_64","Amp_128","Amp_256","Amp_512","Amp_1024","Amp_2048")]


# # set.seed(123)
# sample_n <- 1000
# dose_sample <- dose_all_df[sample(1:nrow(dose_all_df), sample_n), ]
# my_col <- colorRampPalette(c("white", "red"))(100)
# 
# p=pheatmap(dose_sample[,-1],
#          cluster_rows = TRUE,
#          cluster_cols = FALSE,
#          show_rownames = FALSE,
#          show_colnames = TRUE,
#          color = my_col,
#          main = paste("Heatmap of EGvar for", sample_n, "Randomly Sampled Mutants"),
#          border_color = NA
#          
#          )
# print(p)
# # ggsave("figures/story/dose_response_EGvar_heatmap.png", p, width = 10, height = 8, dpi = 300)
# 
# 
# # plot lines
# dose_long <- dose_sample %>%
#   pivot_longer(cols = starts_with("Amp_"), names_to = "Dose", values_to = "EGvar")
# 
# # filter out EEGvar <= 0.05
# # dose_long <- dose_long[dose_long$EEGvar > 0.05, ]
# 
# dose_long$Dose <- as.numeric(gsub("Amp_", "", dose_long$Dose))
# p=ggplot(dose_long, aes(x = Dose, y = EGvar, group = site)) +
#   geom_hline(yintercept = 0.00, linetype = "dashed", color = "grey") +
#   geom_hline(yintercept = 1.00, linetype = "dashed", color = "grey") +
#   geom_line(alpha = 0.1,color="red") +
#   # geom_point(alpha = 1,color="blue", size=1) +
#   # violin plot to show density
#   # set no edge color and light blue fill color, bigger width
#   geom_violin(aes(x = Dose, y = EGvar,group = Dose),fill = "blue", alpha = 0.8, color = NA,width = 0.6)+
#   # coord_cartesian(xlim = c(2, 200))+
#   scale_x_log10(breaks = c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)) +
#   labs(title = "Dose-Response Curves for Randomly Sampled Mutants",
#        x = "Dose (log scale)",
#        y = "EGvar") +
#   theme_classic()
# # ggsave("figures/story/dose_response_EGvar_lines.png", p, width = 10, height = 8, dpi = 300)   



###################### compare mini library solo v.s. mix ####################
# read mini library
mini_mutant=read.csv("minilib_BC_PLATE_mutant.stat",sep="\t",header=FALSE)
colnames(mini_mutant)=c("seq","Mutant","PCR","site")
dose_mini_solo_df=read.csv("Dose_minilib_summary.csv",sep=",",header=TRUE)
# use only AMP_ex <8's data
# dose_mini_solo_df=subset(dose_mini_solo_df, Amp_ex <8)
# with more Amp_ex >=8 data
dose_mini_solo_more_df = read.csv("35C_BC_genotype_TEM1minilib_Gvar.tsv", sep="\t", header=TRUE)
dose_mini_solo_more_df=dose_mini_solo_more_df[,c("Mutant","Amp","EGvar")]
colnames(dose_mini_solo_more_df)[colnames(dose_mini_solo_more_df)=="Amp"]="Amp_ex"
colnames(dose_mini_solo_more_df)[colnames(dose_mini_solo_more_df)=="EGvar"]="EGvar_mini_mix"
# set colname Amp to Amp_ex
# colnames(dose_mini_solo_more_df)[colnames(dose_mini_solo_more_df)=="Amp"]="Amp_ex"
# get together for 0 ~ 2048
# dose_mini_solo_more_df=subset(dose_mini_solo_more_df, Amp_ex >=8)
# dose_mini_solo_more_df=dose_mini_solo_more_df[,c("Mutant","Amp_ex","EGvar")]
dose_mini_solo_df=dose_mini_solo_df[,c("Mutant","Amp_ex","EGvar")]
# dose_mini_solo_df=rbind(dose_mini_solo_df,dose_mini_solo_more_df)
# rm(dose_mini_solo_more_df)

mini_mutant=merge(mini_mutant,dose_all_df,by="site",all.x=TRUE)
dose_mini_mix_df=mini_mutant[,c("Mutant","site","seq","Amp_0","Amp_2","Amp_4","Amp_8","Amp_16","Amp_32","Amp_64","Amp_128","Amp_256","Amp_512","Amp_1024","Amp_2048")]
dose_mini_mix_df=na.omit(dose_mini_mix_df)
# remove replicates
dose_mini_mix_df=dose_mini_mix_df[!duplicated(dose_mini_mix_df$Mutant), ]
# convert into long format
dose_mini_mix_df <- dose_mini_mix_df %>%
  pivot_longer(cols = starts_with("Amp_"), names_to = "Amp_ex", values_to = "EGvar")
dose_mini_mix_df$Amp_ex <- as.numeric(gsub("Amp_", "", dose_mini_mix_df$Amp_ex))


## merge solo and mix data by both Mutant and Amp_ex, preserved all rows in mix data
# rename EGvar to EGvar_solo/EGvar_mix
colnames(dose_mini_solo_df)[colnames(dose_mini_solo_df)=="EGvar"]="EGvar_solo"
colnames(dose_mini_mix_df)[colnames(dose_mini_mix_df)=="EGvar"]="EGvar_mix"
# merge
growth_df=merge(dose_mini_mix_df,dose_mini_solo_df,by=c("Mutant","Amp_ex"),all.x=TRUE,all.y=TRUE)
# growth_df=na.omit(growth_df)
# calculate ratio
growth_df$ratio=growth_df$EGvar_mix/growth_df$EGvar_solo
# check how many mutants & which mutants could be used at different doses
growth_n_df=growth_df %>%
  group_by(Amp_ex) %>%
  summarise(
    total_mutants = n(),
    mutants_with_solo_data = sum(!is.na(EGvar_solo)),
    percent_with_solo_data = mutants_with_solo_data / total_mutants * 100
  )
dose_mini_mix_df=as.data.frame(dose_mini_mix_df)
dose_mini_solo_df=as.data.frame(dose_mini_solo_df)
p=ggplot()+
  geom_hline(yintercept = 1.00, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  geom_line(data=dose_mini_mix_df, aes(x=Amp_ex, y=EGvar_mix, group=Mutant, color=Mutant), alpha=0.75,linetype="dashed")+
  # geom_point(data=dose_mini_mix_df, aes(x=Amp_ex, y=EGvar_mix, group=Mutant), alpha=0.5, size=1)+
  geom_line(data=dose_mini_solo_df, aes(x=Amp_ex, y=EGvar_solo, group=Mutant, color=Mutant), alpha=0.75)+
  # geom_point(data=dose_mini_solo_df, aes(x=Amp_ex, y=EGvar_solo, group=Mutant), alpha=0.5, size=1)+
  # geom_text(data=growth_n_df, aes(x=Amp_ex, y=0.5, label=paste0(mutants_with_solo_data,"/",total_mutants)), color="black", size=3)+
  # scale_x_log10(breaks = c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)) +
  scale_color_viridis_d() +
  facet_wrap(~ Mutant, ncol=6, scales = "free_y") +
  labs(title = "Dose-Response Curves for Mini-Library Mutants (mix & solo)",
       x = "Dose (log scale)",
       y = "EGvar") +
  theme_classic()+
  theme(legend.position = "none")
# print(p)
# ggsave("figures/story/dose_response_EGvar_mix_and_solo_all_mutants_nolog.png", p, width = 16, height = 16, dpi = 500)



calib_mutant <- function(df){
  # df: growth_df for one mutant
  max_EGvar=max(df$EGvar, na.rm=TRUE)
  max_dose=df$Amp_ex[which.max(df$EGvar)]
  # set all EGvar at dose <= max_dose to max_EGvar
  df$EGvar[df$Amp_ex <= max_dose]=max_EGvar
  # then normalize to <=1
  EGvar_calib=df$EGvar/max_EGvar
  return(EGvar_calib)
}
dose_mini_mix_df$EGvar_mix_calib=NA
for (mutant in unique(dose_mini_mix_df$Mutant)){
  df=subset(dose_mini_mix_df, Mutant==mutant)
  colnames(df)[colnames(df)=="EGvar_mix"]="EGvar"
  # calibrate
  df$EGvar_mix_calib=calib_mutant(df)
  dose_mini_mix_df$EGvar_mix_calib[dose_mini_mix_df$Mutant==mutant]=df$EGvar_mix_calib
}
# compare EGvar_mix_calib & EGvar_mix & EGvar_solo
p=ggplot()+
  geom_hline(yintercept = 1.00, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  geom_line(data=dose_mini_solo_df, aes(x=Amp_ex, y=EGvar_solo, group=Mutant, color=Mutant), alpha=0.25)+
  geom_point(data=dose_mini_solo_df, aes(x=Amp_ex, y=EGvar_solo, group=Mutant, color=Mutant),size=1, shape=1)+
  geom_line(data=dose_mini_mix_df, aes(x=Amp_ex, y=EGvar_mix, group=Mutant, color=Mutant), alpha=0.75,linetype="dashed")+
  geom_point(data=dose_mini_mix_df, aes(x=Amp_ex, y=EGvar_mix, group=Mutant, color=Mutant), size=1, shape=2)+
  geom_line(data=dose_mini_mix_df, aes(x=Amp_ex, y=EGvar_mix_calib, group=Mutant, color=Mutant), alpha=0.95,linetype="dotted")+
  geom_point(data=dose_mini_mix_df, aes(x=Amp_ex, y=EGvar_mix_calib, group=Mutant, color=Mutant),size=1, shape=3)+
  # scale_x_log10(breaks = c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)) +
  scale_color_viridis_d() +
  facet_wrap(~ Mutant, ncol=6, scales = "free_y") +
  labs(title = "Dose-Response Curves for Mini-Library Mutants (mix & solo)",
       x = "Dose (log scale)",
       y = "EGvar") +
  theme_classic()+
  theme(legend.position = "none")+
  # add legend for linetype
  guides(linetype = guide_legend(override.aes = list(color = "black")))
# print(p)
# ggsave("figures/story/dose_response_EGvar_mix_and_solo_all_mutants_calibrated_nolog.png", p, width = 16, height = 16, dpi = 500)




calib_mutant_wide <- function(df, dose_cols) {
  mat <- as.matrix(df[, dose_cols])
  colnames(mat) <- dose_cols
  diff_2_4 <- df$Amp_4 - df$Amp_2
  
  mat_new <- t(sapply(seq_len(nrow(mat)), function(i) {
    x <- mat[i, ]
    if (diff_2_4[i] > 0) {
      # high activity
      j <- which.max(x)
      x[1:j] <- max(x, na.rm = TRUE)
      x / max(x, na.rm = TRUE)
    } else {
      # low activity
      max_after_4 <- max(x[names(x) %in% dose_cols[dose_cols != "Amp_0" & dose_cols != "Amp_2"]], na.rm = TRUE)
      x["Amp_0"] <- max_after_4
      x["Amp_2"] <- max_after_4
      x / max_after_4
    }
  }))
  colnames(mat_new) <- dose_cols
  df[, dose_cols] <- as.data.frame(mat_new)
  return(df)
}
dose_cols <- grep("^Amp_", names(dose_all_df), value = TRUE)
dose_all_df_calib <- calib_mutant_wide(dose_all_df, dose_cols)

# round to 4 decimal places
dose_all_df_calib[, dose_cols] <- round(dose_all_df_calib[, dose_cols], 4)


dose_all_df_calib_diff <- data.frame(site = dose_all_df_calib$site)
for (i in 1:(length(dose_cols)-1)){
  col1 <- dose_cols[i]
  col2 <- dose_cols[i+1]
  new_col <- paste0(col1, "_to_", col2, "_diff")
  dose_all_df_calib_diff[[new_col]] <- dose_all_df_calib[[col2]] - dose_all_df_calib[[col1]]
}
# filter out mutants with sharp increase (diff > 0.3) at any dose interval
sharp_increase_sites <- unique(unlist(lapply(dose_all_df_calib_diff[, -1], function(x) dose_all_df_calib_diff$site[which(x > 0.1)])))
dose_all_df_calib_filter=dose_all_df_calib[!dose_all_df_calib$site %in% sharp_increase_sites, ]
print(paste("Number of mutants filtered out due to sharp increase:", length(sharp_increase_sites)))




# randomly check
set.seed(123)
sample_n <- 1000
dose_sample <- dose_all_df[sample(1:nrow(dose_all_df), sample_n), ]
dose_sample_calib <- dose_all_df_calib[sample(1:nrow(dose_all_df_calib), sample_n), ]
dose_sample_calib_filter <- dose_all_df_calib_filter[sample(1:nrow(dose_all_df_calib_filter), sample_n), ]

dose_sample_long <- dose_sample_calib_filter %>%
  pivot_longer(cols = starts_with("Amp_"), names_to = "Dose", values_to = "EGvar")
dose_sample_long$Dose <- as.numeric(gsub("Amp_", "", dose_sample_long$Dose))
dose_sample_calib_long <- dose_sample_calib %>%
  pivot_longer(cols = starts_with("Amp_"), names_to = "Dose", values_to = "EGvar")
dose_sample_calib_long$Dose <- as.numeric(gsub("Amp_", "", dose_sample_calib_long$Dose))
p1=ggplot(dose_sample_long, aes(x = Dose, y = EGvar, group = site)) +
  geom_hline(yintercept = 0.00, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 1.00, linetype = "dashed", color = "grey") +
  geom_line(alpha = 0.1,color="red") +
  scale_x_log10(breaks = c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)) +
  labs(title = "Dose-Response Curves for Randomly Sampled Mutants",
       x = "Dose (log scale)",
       y = "EGvar") +
  theme_classic()
print(p1)

library(gridExtra)
# p=grid.arrange(p1, p2, ncol=1)
# ggsave("figures/story/dose_response_EGvar_lines_calibrated_comparison.png", p, width = 12, height = 16, dpi = 300)



judge_activity <- function(df) {
  df$activity <- ifelse(df$Amp_4 > df$Amp_2, "high", "low")
  return(df$activity)
}
dose_all_df_calib$activity <- judge_activity(dose_all_df)
# check how many high/low activity mutants
table(dose_all_df_calib$activity)

# 
# 
# # plot EGvar_mix & EGvar_solo v.s. Amp_ex, with different colors for different mutants, in one plot
# mutant_ROI=c("D11","B1","D9","D10","A11")
# growth_df_ROI=subset(growth_df, Mutant %in% mutant_ROI)
# p=ggplot(growth_df, aes(x = Amp_ex)) +
#   geom_hline(yintercept = 1.00, linetype = "dashed", color = "grey") +
#   geom_line(aes(y = EGvar_mix, color = paste(Mutant, "mix")), alpha = 0.7) +
#   geom_point(aes(y = EGvar_mix, color = paste(Mutant, "mix")), alpha = 1, size=1) +
#   geom_line(aes(y = EGvar_solo, color = paste(Mutant, "solo")), alpha = 0.7, linetype = "dashed") +
#   geom_point(aes(y = EGvar_solo, color = paste(Mutant, "solo")), alpha = 1, size=1, shape=17) +
#   scale_x_log10(breaks = c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)) +
#   scale_color_viridis_d()+
#   labs(title = "EGvar (mix & solo) vs Dose for Selected Mini-Library Mutants",
#        x = "Dose (log scale)",
#        y = "EGvar",
#        color = "Mutant and Condition") +
#   theme_bw()
# print(p)
# ggsave("figures/story/dose_response_EGvar_mix_and_solo_selected_mutants.png", p, width = 10, height = 6, dpi = 300)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# growth_df=subset(growth_df, Mutant !="D11")
# # plot ratio v.s. Amp_ex, with different colors for different mutants
# p=ggplot(growth_df, aes(x = Amp_ex, y = ratio, group = Mutant, color = Mutant)) +
#   geom_hline(yintercept = 1.00, linetype = "dashed", color = "grey") +
#   geom_line(alpha = 0.7) +
#   geom_point(alpha = 1, size=1) +
#   # violin plot to show density
#   # set no edge color and light blue fill color, bigger width
#   geom_violin(aes(x = Amp_ex, y = ratio, group = Amp_ex),fill = "blue", alpha = 0.75, color = NA,width = 0.5)+
#   # coord_cartesian(xlim = c(2, 200))+
#   scale_x_log10(breaks = c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)) +
#   scale_color_viridis_d() +
#   labs(title = "Ratio of EGvar (mix/solo) for Mini-Library Mutants",
#        x = "Dose (log scale)",
#        y = "EGvar ratio (mix/solo)") +
#   theme_classic()
# print(p)
# # ggsave("figures/story/dose_response_EGvar_ratio.png", p, width = 10, height = 8, dpi = 300)
# 
# # plot ratio distribution at different doses
# p=ggplot(growth_df, aes(x = ratio, group = Amp_ex, fill = as.factor(Amp_ex))) +
#   geom_histogram(binwidth = 0.2, position = "identity", alpha = 0.6) +
#   # geom_density(aes(y = ..count.. ), alpha = 0.2, color = "black", linewidth=0.5) +
#   scale_fill_viridis_d() +
#   facet_wrap(~ Amp_ex, ncol=4, scales = "free_y") +
#   labs(title = "Distribution of EGvar Ratio (mix/solo) at Different Doses",
#        x = "EGvar ratio (mix/solo)",
#        y = "Count",
#        fill = "Dose") +
#   theme_classic()+
#   theme(legend.position = "none")
# print(p)
# ggsave("figures/story/dose_response_EGvar_ratio_distribution.png", p, width = 12, height = 8, dpi = 300)
# 
# 
# # plot ratio v.s. EGvar_solo for different doses, with different colors for different doses, separate plots for different mutants
# p=ggplot(growth_df, aes(x = EGvar_solo, y = EGvar_mix, color = Mutant)) +
#   # geom_hline(yintercept = 1.00, linetype = "dashed", color = "grey") +
#   geom_point(alpha = 0.85, size=1) +
#   stat_smooth(method = "loess", se = TRUE, color = "lightblue",linewidth=0.75, alpha=0.2) +
#   # scale_x_log10() +
#   scale_color_viridis_d() +
#   labs(title = "EGvar_mix vs EGvar_solo for Mini-Library Mutants",
#        x = "EGvar_solo",
#        y = "EGvar_mix") +
#   theme_classic()# +
#   # facet_wrap(~ Amp_ex, scales = "free")
# print(p)
# # ggsave("figures/story/dose_response_EGvar_mix_vs_solo_together.png", p, width = 10, height = 6, dpi = 300)
# # # ggsave("figures/story/dose_response_EGvar_mix_vs_solo_saperate.png", p, width = 16, height = 12, dpi = 500)
# 
# 
# p=ggplot(growth_df, aes(x = EGvar_solo, y = ratio, color = Mutant)) +
#   # geom_hline(yintercept = 1.00, linetype = "dashed", color = "grey") +
#   geom_point(alpha = 0.85, size=1) +
#   stat_smooth(method = "loess", se = TRUE, color = "lightblue",linewidth=0.75, alpha=0.2) +
#   # scale_x_log10() +
#   # xlim(0.1,0.6)+
#   scale_color_viridis_d() +
#   labs(title = "ratio vs EGvar_solo for Mini-Library Mutants",
#        x = "EGvar_solo",
#        y = "ratio") +
#   theme_classic()+
#   facet_wrap(~ Amp_ex, scales = "free")
# print(p)
# # ggsave("figures/story/dose_response_EGvar_ratio_vs_1_together.png", p, width = 10, height = 6, dpi = 300)
# # # ggsave("figures/story/dose_response_EGvar_ratio_vs_1_saperate.png", p, width = 16, height = 12, dpi = 500)
# 
# ## use trimmed mean to calculate the average ratio at different doses
# ratio_mean_df=growth_df %>%
#   group_by(Amp_ex) %>%
#   summarise(
#     mean_ratio = mean(ratio, na.rm = TRUE),
#     trimmed_mean_ratio = mean(ratio, trim = 0.2, na.rm = TRUE),
#     median_ratio = median(ratio, na.rm = TRUE),
#     sd_ratio = sd(ratio, na.rm = TRUE),
#     n = n()
#   )
# # plot mean ratio v.s. Amp_ex
# p=ggplot(ratio_mean_df, aes(x = Amp_ex)) +
#   geom_line(data=growth_df, aes(x = Amp_ex, y = ratio, group = Mutant), color = "lightgrey", alpha = 0.5) +
#   geom_hline(yintercept = 1.00, linetype = "dashed", color = "grey") +
#   geom_line(aes(y = mean_ratio, color = "Mean"), size=1) +
#   geom_point(aes(y = mean_ratio, color = "Mean"), size=2) +
#   geom_line(aes(y = trimmed_mean_ratio, color = "Trimmed Mean"), size=1) +
#   geom_point(aes(y = trimmed_mean_ratio, color = "Trimmed Mean"), size=2) +
#   geom_line(aes(y = median_ratio, color = "Median"), size=1) +
#   geom_point(aes(y = median_ratio, color = "Median"), size=2) +
#   scale_x_log10(breaks = c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)) +
#   scale_color_viridis_d() +
#   labs(title = "Average EGvar Ratio (mix/solo) at Different Doses",
#        x = "Dose (log scale)",
#        y = "Average EGvar ratio (mix/solo)",
#        color = "Statistic") +
#   theme_classic()
# print(p)
# ggsave("figures/story/dose_response_EGvar_ratio_average.png", p, width = 10, height = 6, dpi = 300)
# 
# 
# ## use trimmed mean to calibrate the mix data to solo data
# # calculate calibration factor at each dose
# ratio_mean_df$calibration_factor=1/ratio_mean_df$median_ratio
# # merge calibration factor to growth_df
# growth_df=merge(growth_df,ratio_mean_df[,c("Amp_ex","calibration_factor")],by="Amp_ex",all.x=TRUE)
# # calculate calibrated EGvar_mix
# growth_df$EGvar_mix_calibrated=growth_df$EGvar_mix*growth_df$calibration_factor
# # plot calibrated EGvar_mix v.s. Amp_ex, with different colors for different mutants
# p=ggplot(growth_df, aes(x = Amp_ex, y = EGvar_mix_calibrated, group = Mutant, color = Mutant)) +
#   geom_hline(yintercept = 1.00, linetype = "dashed", color = "grey") +
#   geom_line(alpha = 0.7) +
#   geom_point(alpha = 1, size=1) +
#   # violin plot to show density
#   # set no edge color and light blue fill color, bigger width
#   geom_violin(aes(x = Amp_ex, y = EGvar_mix_calibrated, group = Amp_ex),fill = "blue", alpha = 0.75, color = NA,width = 0.5)+
#   # coord_cartesian(xlim = c(2, 200))+
#   scale_x_log10(breaks = c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)) +
#   scale_color_viridis_d() +
#   labs(title = "Calibrated EGvar_mix vs Dose for Mini-Library Mutants",
#        x = "Dose (log scale)",
#        y = "Calibrated EGvar_mix") +
#   theme_classic()
# print(p)





################ clustering ####################
################ clustering ####################
# # set.seed(123)
sample_n <- 2000
dose_sample <- dose_all_df_calib_filter[sample(1:nrow(dose_all_df_calib_filter), sample_n), ]
rownames(dose_sample)=dose_sample$site
dose_sample=subset(dose_sample, select = -c(site))
dose_sample <- na.omit(dose_sample)
mat <- as.matrix(dose_sample)
# (1) 欧氏距离
# dist_euc <- dist(mat, method = "euclidean")
# (2) 相关距离 (推荐)
cor_dist <- as.dist(1 - cor(t(mat), method = "pearson"))
hc <- hclust(cor_dist, method = "ward.D2")  # 也可以试 "average" 或 "complete"
clusters <- cutree(hc, k = 12)   # 例如切成4类
table(clusters)

# line plot for each cluster
dose_sample$site <- rownames(dose_sample)
dose_sample$cluster <- as.factor(clusters)

dose_long <- dose_sample %>%
  pivot_longer(cols = starts_with("Amp_"), names_to = "Dose", values_to = "EGvar")
dose_long$Dose <- as.numeric(gsub("Amp_", "", dose_long$Dose))
p=ggplot(dose_long, aes(x = Dose, y = EGvar, group = site, color = cluster)) +
  geom_hline(yintercept = 0.00, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 25.00, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 125.00, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 1.00, linetype = "dashed", color = "grey") +
  # smooth line for each mutant
  # stat_smooth(method = "loess", se = FALSE, alpha = 0.1, linewidth=0.1, span=0.3) +
  geom_line(alpha = 0.1) +
  # geom_point(alpha = 1,color="blue", size=1) +
  # violin plot to show density
  # set no edge color and light blue fill color, bigger width
  # geom_violin(aes(x = Dose, y = EGvar,group = Dose),fill = "blue", alpha = 0.8, color = NA,width = 0.6)+
  coord_cartesian(xlim = c(0, 200))+
  # scale_x_log10(breaks = c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)) +
  scale_color_viridis_d() +
  facet_wrap(~ cluster, ncol=4, scales = "free_y") +
  labs(title = "Dose-Response Curves for Randomly Sampled Mutants",
       x = "Dose (log scale)",
       y = "EGvar") +
  # remove legend
  theme_classic()+
  theme(legend.position = "none")
print(p)

























###################### Manually grab the pattern (fixed summary function) ####################
library(data.table)
# ---------- 参数 (保持不变) ----------
control_threshold <- 0.10
low_drop_thr <- 0.50
low_min_low <- 0.05
low_min_high <- 0.50
mid_cv_thr <- 0.5
plateau_close_thr <- 0.20

n_cand_plot <- 1000
n_neg_sample_plot <- 3000
set.seed(421)

# ---------- 数据准备 (保持不变) ----------
dt <- as.data.table(dose_all_df_calib_filter)
low_cols  <- c("Amp_0","Amp_2","Amp_4","Amp_8", "Amp_16")
mid_cols  <- c("Amp_16","Amp_32","Amp_64","Amp_128")
high_cols <- c("Amp_32","Amp_64","Amp_128","Amp_256","Amp_512","Amp_1024","Amp_2048")
all_cols  <- setdiff(names(dt), "site")

# ---------- summary 指标 (计算 low_drop 和 low_steep_rise) ----------
dt[, low_min := do.call(pmin, .SD), .SDcols = low_cols]
dt[, low_max := do.call(pmax, .SD), .SDcols = low_cols]

# 1. low_drop 使用 Amp_0 和 Amp_16 的差值计算
dt[, low_drop := (Amp_0 - Amp_16) ]

dt[, mid_sd   := apply(.SD, 1, sd,  na.rm = TRUE), .SDcols = mid_cols]
dt[, mid_mean := apply(.SD, 1, mean, na.rm = TRUE), .SDcols = mid_cols]
dt[, mid_cv   := fifelse(mid_mean > 0, mid_sd / mid_mean, NA_real_)]

dt[, high_max := do.call(pmax, .SD), .SDcols = high_cols]

# # 2. 新增低浓度段陡升检查：所有相邻点增量都 <= 0.1
# dt[, low_steep_rise :=
#      (Amp_2 - Amp_0 <= 0.1) &
#      (Amp_4 - Amp_2 <= 0.1) &
#      (Amp_8 - Amp_4 <= 0.1) &
#      (Amp_16 - Amp_8 <= 0.1)]


# ---------- negative control (新增 low_steep_rise 约束) ----------
neg_ctrl_dt <- dt[
  (low_min < control_threshold) & 
    (high_max < control_threshold),# & 
    # (low_steep_rise == TRUE), # 排除低浓度段陡升的样本
  .(site, low_min, high_max)
]

# ---------- candidates (保持不变) ----------
candidates_dt <- dt[
  low_drop > low_drop_thr &
    low_min > low_min_low & low_min < low_min_high &
    mid_cv < mid_cv_thr &
    abs(mid_mean - low_min) < plateau_close_thr,# &
    # low_steep_rise == TRUE, 
  .(site, low_min, low_max, low_drop, mid_mean, mid_sd, mid_cv)
]
cat("候选 mutants 数量:", nrow(candidates_dt), "\n")
cat("negative controls 数量:", nrow(neg_ctrl_dt), "\n")


head(candidates_dt)
head(neg_ctrl_dt)


# ---------- 抽样准备 ----------
cand_plot_sites <- head(candidates_dt$site, n_cand_plot)
neg_plot_sites  <- if (nrow(neg_ctrl_dt) <= n_neg_sample_plot) neg_ctrl_dt$site else sample(neg_ctrl_dt$site, n_neg_sample_plot)

# ---------- 稳健的 summary (median ± sd) 函数 ----------
make_summary_long_sd <- function(dt_subset, label, base_dt, cols) {
  if (is.null(dt_subset) || nrow(dt_subset) == 0) return(NULL)
  
  # 计算 median / sd
  med_dt <- base_dt[site %in% dt_subset$site, lapply(.SD, median, na.rm = TRUE), .SDcols = cols]
  sd_dt  <- base_dt[site %in% dt_subset$site, lapply(.SD, sd, na.rm = TRUE), .SDcols = cols]
  
  med_vec <- as.numeric(med_dt[1,])
  sd_vec  <- as.numeric(sd_dt[1,])
  dose_names <- cols
  
  res <- data.table(
    Dose = as.numeric(sub("^Amp_", "", dose_names)),
    med  = med_vec,
    sd   = sd_vec,
    type = label
  )
  setorder(res, Dose)
  return(res)
}

# ---------- 计算 neg / cand 的 summary long (median ± SD) ----------
neg_summary_long  <- make_summary_long_sd(neg_ctrl_dt,  "Neg",  dt, all_cols)
cand_summary_long <- make_summary_long_sd(candidates_dt,"Cand", dt, all_cols)



to_melt_sites <- c(cand_plot_sites, neg_plot_sites)
melt_dt <- data.table::melt(dt[site %in% to_melt_sites, c("site", all_cols), with = FALSE], 
                            id.vars = "site", variable.name = "Dose", value.name = "EGvar") 
melt_dt[, Dose := as.numeric(sub("Amp_", "", Dose))] 
melt_dt[, type := ifelse(site %in% cand_plot_sites, "Candidate", "Neg_sample")]
# ---------- 绘图 (median ± SD with error bars) ----------
melt_dt[, Dose_f := factor(Dose, levels = sort(unique(Dose)))]

p <- ggplot() +
  # sampled negative individuals
  geom_line(data = melt_dt[type == "Neg_sample"], aes(x = Dose_f, y = EGvar, group = site),
            color = "grey70", linetype = "dashed", alpha = 0.1) +
  # candidate individual curves
  geom_line(data = melt_dt[type == "Candidate"], aes(x = Dose_f, y = EGvar, color = site, group = site),
            alpha = 0.5) +
  # median lines
  stat_summary(data = melt_dt[type == "Candidate"], aes(x = as.numeric(Dose_f), y = EGvar),
               fun = median, geom = "line", color = "red", size = 1.2) +
  stat_summary(data = melt_dt[type == "Neg_sample"], aes(x = as.numeric(Dose_f), y = EGvar),
               fun = median, geom = "line", color = "black", size = 1.2, linetype = "dotted") +
  # violin plot
  geom_violin(data = melt_dt[type == "Candidate"], aes(x = Dose_f, y = EGvar),
              fill = "red", alpha = 0.5, color = NA) +
  geom_violin(data = melt_dt[type == "Neg_sample"], aes(x = Dose_f, y = EGvar),
              fill = "black", alpha = 0.5, color = NA) +
  scale_color_viridis_d() +
  theme_classic() +
  labs(title = "Representative candidates vs negative controls",
       subtitle = sprintf("candidates: %d; negs: %d; median trend with violin distribution",
                          nrow(candidates_dt), nrow(neg_ctrl_dt)),
       x = "Ampicillin (concentration)", y = "EGvar") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        # make legend smaller
        legend.text = element_text(size = 5),
        # make space between legend rows smaller
        legend.spacing.y = unit(0.025, "cm"))

print(p)
# ggsave(p, file = "figures/story/dose_response_curve_candidates_and_nc.png", width = 12, height = 8, dpi = 500)


##################### mutation sites visualization ####################

# 氨基酸分类
aa_groups <- list(
  Nonpolar = c("A","V","L","I","M","F","W","P","G"),
  Polar    = c("S","T","C","N","Q","Y"),
  Positive = c("K","R","H"),
  Negative = c("D","E")
)
aa_group_colors <- c(
  Nonpolar = "#FDB863",   # 橙
  Polar    = "#80CDC1",   # 青绿
  Positive = "#B2182B",   # 红
  Negative = "#2166AC"   # 蓝
)
# set factor levels
aa_group_levels <- names(aa_groups)
aa_group_colors <- aa_group_colors[aa_group_levels]



aa2group <- function(aa) {
  for (grp in names(aa_groups)) {
    if (aa %in% aa_groups[[grp]]) return(grp)
  }
  return(NA_character_)
}

analyze_mut_sites <- function(site_vec, label="Group", wt_seq = NULL, label_every = 5) {
  # parse inputs
  dt_sites <- data.table(site = site_vec)
  # explode comma-separated mutation lists into rows
  dt_mut <- dt_sites[, .(mutation = unlist(strsplit(site, ","))), by = .(site)]
  # split mutation tokens into WT, Pos, Mut
  dt_mut[, c("WT","Pos","Mut") := tstrsplit(mutation, "_")]
  dt_mut[, Pos := as.integer(Pos)]
  # add AA group labels
  dt_mut[, WT_group  := vapply(WT, aa2group, character(1))]
  dt_mut[, Mut_group := vapply(Mut, aa2group, character(1))]
  dt_mut[, label := label]
  
  # Determine sequence length L from wt_seq or from observed positions
  if (!is.null(wt_seq)) {
    if (length(wt_seq) == 1) wt_seq <- strsplit(wt_seq, "")[[1]]
    L <- length(wt_seq)
    wt_vec <- wt_seq
  } else {
    L <- max(dt_mut$Pos, na.rm = TRUE)
    wt_vec <- rep(NA_character_, L)
  }
  
  # build counts per position (include zero counts)
  counts_dt <- data.table(Pos = seq_len(L))
  pos_counts <- dt_mut[, .N, by = Pos]
  counts_dt <- merge(counts_dt, pos_counts, by = "Pos", all.x = TRUE)
  counts_dt[is.na(N), N := 0]
  # 把 N 转换为频率
  counts_dt[, freq := N / sum(N)]

  # chunk into two panels (top = chunk 1, bottom = chunk 2)
  split_pt <- ceiling(L/2)
  counts_dt[, chunk := ifelse(Pos <= split_pt, 1L, 2L)]
  
  # build seq_df for labels (every label_every positions)
  seq_df <- data.table(Pos = seq_len(L), aa = wt_vec)
  seq_df[, chunk := ifelse(Pos <= split_pt, 1L, 2L)]
  seq_df[, label_flag := ( (Pos - 1) %% label_every == 0 )]  # label every label_every positions
  
  # P1: per-position barplot (binwidth=1) split into two rows
  ymax <- max(counts_dt$N, 1)
  active_sites <- c(68, 164)
  # 把 WT_group 加到 seq_df 里
  seq_df[, WT_group := vapply(aa, aa2group, character(1))]
  
  # 合并 counts_dt 和 WT_group 信息，方便给柱子上色
  plot_dt <- merge(counts_dt, seq_df[, .(Pos, WT_group, aa)], by = "Pos", all.x = TRUE)
  
  # 只在每隔 5 个位置显示数字
  seq_df[, show_label := ifelse(Pos %% 5 == 0, Pos, "")]
  
  
  # 更新 ymax 基于频率
  ymax <- max(counts_dt$freq, 1e-6)
  
  p1 <- ggplot(plot_dt, aes(x = Pos, y = freq, fill = WT_group)) +
    geom_col(width = 0.9, color = NA) +
    facet_wrap(~chunk, ncol = 1, scales = "free_x") +
    geom_text(data = seq_df, aes(x = Pos, y = -0.08 * ymax, label = aa, color = WT_group),
              inherit.aes = FALSE, size = 3, vjust = 1) +
    geom_text(data = seq_df, aes(x = Pos, y = -0.18 * ymax, label = show_label),
              inherit.aes = FALSE, size = 2.5, vjust = 1) +
    geom_point(data = plot_dt[Pos %in% active_sites],
               aes(x = Pos, y = freq + 0.05 * ymax),
               shape = 24, size = 1.5, fill = "red", color = "white", inherit.aes = FALSE) +
    scale_fill_manual(values = aa_group_colors, name = "AA group") +
    scale_color_manual(values = aa_group_colors, guide = "none") +
    scale_y_continuous(limits = c(-0.25 * ymax, ymax * 1.25), expand = c(0,0)) +
    labs(title = paste("Mutation position frequency (with WT sequence) -", label),
         x = "Residue position", y = "Frequency") +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = unit(c(1,1,3,1), "lines"),
      legend.position = "bottom"
    ) +
    coord_cartesian(clip = "off")
  
  
  # P2: top substitution counts (top 20)
  repl_counts <- dt_mut[, .N, by = .(WT, Mut)][order(-N)]
  repl_counts[, sub := paste0(WT, "→", Mut)]
  top_repl <- repl_counts[1:min(50, .N)]
  p2 <- ggplot(top_repl, aes(x = reorder(sub, N), y = N)) +
    geom_col(fill = "darkorange") + coord_flip() +
    labs(title = paste("Top AA substitutions -", label), x = "Substitution", y = "Count") +
    theme_bw()
  
  # P3: property transitions
  prop_dt <- dt_mut[, .N, by = .(WT_group, Mut_group)]
  # filter out rows with NA groups
  prop_dt <- prop_dt[!is.na(WT_group) & !is.na(Mut_group)]
  p3 <- ggplot(prop_dt, aes(x = WT_group, y = N, fill = Mut_group)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = paste("AA property transitions -", label),
         x = "WT group", y = "Count", fill = "Mut group") +
    scale_fill_manual(values = aa_group_colors) +
    theme_bw()
  
  # return parsed data and plots
  return(list(
    parsed = dt_mut,
    counts_per_pos = counts_dt,
    seq_map = seq_df,
    plots = list(site_dist = p1, substitution = p2, property = p3)
  ))
}

wt_seq="MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW"
# candidates
cand_sites <- candidates_dt$site
cand_res <- analyze_mut_sites(cand_sites, label="Candidates", wt_seq=wt_seq)
# negative controls
neg_sites <- neg_ctrl_dt$site
neg_res <- analyze_mut_sites(neg_sites, label="Neg controls", wt_seq=wt_seq)

# all
all_sites <- dt$site
all_res <- analyze_mut_sites(all_sites, label="All mutants", wt_seq=wt_seq)



# low activity
low_activity_sites <- dose_all_df_calib[dose_all_df_calib$activity == "low", "site"]
low_activity_res <- analyze_mut_sites(low_activity_sites, label="Low activity mutants", wt_seq=wt_seq)


# high activity
high_activity_sites <- dose_all_df_calib[dose_all_df_calib$activity == "high", "site"]
high_activity_res <- analyze_mut_sites(high_activity_sites, label="High activity mutants", wt_seq=wt_seq)
print(high_activity_res$plots$site_dist)
# # active-site mutated mutants
# acgive_sites: 68, 164
active_sites <- dose_all_df_calib[grep("_68|_164", dose_all_df_calib$site), "site"]
active_res <- analyze_mut_sites(active_sites, label="Active-site mutated mutants", wt_seq=wt_seq)


# # save plots
# ggsave(cand_res$plots$site_dist, file = "figures/mutation/mutation_position_counts_candidates.png", width = 14, height = 6, dpi = 500)
# ggsave(cand_res$plots$substitution, file = "figures/mutation/top_aa_substitutions_candidates.png", width = 6, height = 8, dpi = 500)
# ggsave(cand_res$plots$property, file = "figures/mutation/aa_property_transitions_candidates.png", width = 6, height = 6, dpi = 500)
# 
# ggsave(neg_res$plots$site_dist, file = "figures/mutation/mutation_position_counts_neg_controls.png", width = 14, height = 6, dpi = 500)
# ggsave(neg_res$plots$substitution, file = "figures/mutation/top_aa_substitutions_neg_controls.png", width = 6, height = 8, dpi = 500)
# ggsave(neg_res$plots$property, file = "figures/mutation/aa_property_transitions_neg_controls.png", width = 6, height = 6, dpi = 500)
# 
# ggsave(all_res$plots$site_dist, file = "figures/mutation/mutation_position_counts_all_mutants.png", width = 14, height = 6, dpi = 500)
# ggsave(all_res$plots$substitution, file = "figures/mutation/top_aa_substitutions_all_mutants.png", width = 6, height = 8, dpi = 500)
# ggsave(all_res$plots$property, file = "figures/mutation/aa_property_transitions_all_mutants.png", width = 6, height = 6, dpi = 500)
# 
# ggsave(low_activity_res$plots$site_dist, file = "figures/mutation/mutation_position_counts_low_activity_mutants.png", width = 14, height = 6, dpi = 500)
# ggsave(low_activity_res$plots$substitution, file = "figures/mutation/top_aa_substitutions_low_activity_mutants.png", width = 6, height = 8, dpi = 500)
# ggsave(low_activity_res$plots$property, file = "figures/mutation/aa_property_transitions_low_activity_mutants.png", width = 6, height = 6, dpi = 500)
# 
# ggsave(high_activity_res$plots$site_dist, file = "figures/mutation/mutation_position_counts_high_activity_mutants.png", width = 14, height = 6, dpi = 500)
# ggsave(high_activity_res$plots$substitution, file = "figures/mutation/top_aa_substitutions_high_activity_mutants.png", width = 6, height = 8, dpi = 500)
# ggsave(high_activity_res$plots$property, file = "figures/mutation/aa_property_transitions_high_activity_mutants.png", width = 6, height = 6, dpi = 500)
# 
# ggsave(active_res$plots$site_dist, file = "figures/mutation/mutation_position_counts_active_site_mutated_mutants.png", width = 14, height = 6, dpi = 500)
# ggsave(active_res$plots$substitution, file = "figures/mutation/top_aa_substitutions_active_site_mutated_mutants.png", width = 6, height = 8, dpi = 500)
# ggsave(active_res$plots$property, file = "figures/mutation/aa_property_transitions_active_site_mutated_mutants.png", width = 6, height = 6, dpi = 500)
# 

