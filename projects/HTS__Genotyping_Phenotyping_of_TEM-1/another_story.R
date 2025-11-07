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
CV_cutoff=0.35
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





# select several mutants as ROI
mutant_ROI=c("D9","B1","D10","D11")
df=subset(fitness_amp_ex_df, Mutant %in% mutant_ROI)
# set plot order: D9,B1,D10,D11
df$Mutant=factor(df$Mutant, levels = mutant_ROI)

# plot dose response curves for ROI


# plot dose response curve for each mutant, with error bars
p=ggplot(df, aes(x = Amp_ex, y = EGvar, color = Mutant)) +
  geom_point() +
  geom_errorbar(aes(ymin = EGvar - sd, ymax = EGvar + sd), width = 0.1) +
  geom_line() +
  geom_hline(yintercept = 0.35, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 0.37, label = "0.35", color = "red",size=3) +
  geom_hline(yintercept = 0.20, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 0.22, label = "0.20", color = "red",size=3) +
  geom_hline(yintercept = 0.0, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 0.02, label = "0.00", color = "red",size=3) +
  geom_hline(yintercept = 0.55, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 0.57, label = "0.55", color = "red",size=3) +
  scale_color_brewer(palette = "Set2") +
  # scale_x_log10() +
  # annotation_logticks(sides = "b") +
  # plot every mutant separately
  # show x-axis for each facet
  # facet_wrap(~ Mutant, ncol = 2, scales = "free_x") + 
  # show y 0~1.1
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Dose Response Curve (really low enzyme activity)", x = "Ampicillin Concentration (Amp_ex)", y = "EGvar") +
  theme_classic() 
  # remove legend
  # theme(legend.position = "none")
print(p)
# ggsave("figures/story/dose_response_curve_ROI.png", p, width = 20, height = 15, units = "cm", dpi = 500)
# ggsave("figures/story/dose_response_curve_ROI_separate.png", p, width = 20, height = 15, units = "cm", dpi = 500)
# ggsave("figures/story/dose_response_curve_ROI_separate_cut_log.png", p, width = 20, height = 15, units = "cm", dpi = 500)





