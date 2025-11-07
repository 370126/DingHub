## load tsv files
fitness_amp_df <- read.csv("./35C_BC_genotype_TEM1minilib_Gvar.tsv", sep = "\t")
kinetics_df <- read.csv("./35C_protein_Km_Vmax.tsv", sep = "\t")

library(ggplot2)
library(tidyverse)
library(minpack.lm)
x=subset(fitness_amp_df, Mutant == "D10")$Amp
y=subset(fitness_amp_df, Mutant == "D10")$EGvar

# fitting fitness v.s. amp_p
fit <- nlsLM(
  y ~ 1/(1 + (x/E)^k),
  start = list(E = 60,k=1.5),
  control = nls.lm.control(maxiter = 1024)
)
summary(fit)
temp=summary(fit)$coefficients[,"Pr(>|t|)"][1]
# print R-square
rss <- sum(residuals(fit)^2);tss <- sum((y - mean(y))^2);
rsq <- 1 - rss/tss
coef_fit <- coef(fit)
E_fit <- coef_fit["E"]; k_fit <- coef_fit["k"]
fitness_amp_x <- function(x) {
  1/(1 + (x/E_fit)^k_fit)
  # (1/(1 + (x/E_fit)^k_fit))
}
ggplot(subset(fitness_amp_df, Mutant == "D10"), aes(x = Amp, y = EGvar)) +
  geom_point() +
  stat_function(fun = fitness_amp_x, color = "red", size = 1) +
  theme_classic() +
  labs(title = paste("Fitness vs Amp_p | R-squared: ", round(rsq, 3)), x = "Amp_p", y = "Fitness (EGvar)") +
  xlim(0, 2100) +
  ylim(0, 1)



# guess parameters
fitness_km_vmax_x <- function(x,p2,p4) {
  x1=(x-p4+sqrt((p4-x)^2+4*p2*x))/2
  y=1/(1+(x1/E_fit)^k_fit)
  return(y)
}
p2_try=12
p4_try=2600

amp_p_amp_ex_x <- function(x,p2,p4){
  x1=(x-p4+sqrt((p4-x)^2+4*p2*x))/2
}

diff_fitness_x <- function(x,p2=p2_try,p4=p4_try){
  fitness_km_vmax_x(x,p2,p4)-fitness_amp_x(x)
}

x=seq(1,2200,by=1)
ggplot() +
  geom_line(data=data.frame(x = x), aes(x = x, y = fitness_km_vmax_x(x,p2_try,p4_try)), color = "blue") +
  xlim(0, 2200) +
  theme_classic() +
  labs(title = "Fitness diff between Mutant and NC", x = "Amp_ex", y = "fitness_mutant - fitness_NC")




# fitting fitness v.s. Km, Vmax
mutant_list=list()
for (mutant in unique(fitness_amp_df$Mutant)) {
  df=subset(fitness_amp_df, Mutant == mutant)
  df=df[,c("Amp","EGvar")]
  colnames(df)=c("Amp_ex","EGvar")
  mutant_list[[mutant]]=df
}
# minus D10 EGvar
for (mutant in names(mutant_list)) {
  mutant_list[[mutant]]$EGvar_relative=mutant_list[[mutant]]$EGvar - mutant_list[["D10"]]$EGvar
}
# call list and keep mutant name
fitness_amp_ex_df=do.call(rbind, lapply(names(mutant_list), function(name) {
  df <- mutant_list[[name]]
  df$Mutant <- name
  return(df)
}))
# save new fitness_amp_ex_df
write.csv(fitness_amp_ex_df, "./fitness_amp_ex_df.csv", row.names = FALSE)


# plot all mutants together, with different colors and lines
p=ggplot() +
  # guess fit line
  geom_point(data = fitness_amp_ex_df, aes(x = Amp_ex, y = EGvar_relative, color = Mutant)) +
  # add lines for each mutant
  geom_line(data = fitness_amp_ex_df, aes(x = Amp_ex, y = EGvar_relative, color = Mutant)) +
  theme_classic() +
  labs(title = "Fitness-Fitness_NC vs Amp_ex for all mutants", x = "Amp_ex", y = "Fitness (EGvar)") +
  xlim(0, 2200) +
  # ylim(0, 1)+
  # log-scale x axis
  # scale_x_log10() +
  theme(legend.position = "none") # hide legend for clarity
x=seq(1,2200,by=1)
p=p+geom_line(data=data.frame(x = x), aes(x = x, y = diff_fitness_x(x)), color = "black", size=1)
print(p)

p=ggplot() +
  # guess fit line
  geom_point(data = fitness_amp_ex_df, aes(x = Amp_ex, y = EGvar, color = Mutant)) +
  # add lines for each mutant
  geom_line(data = fitness_amp_ex_df, aes(x = Amp_ex, y = EGvar, color = Mutant)) +
  theme_classic() +
  labs(title = "Fitness-Fitness_NC vs Amp_ex for all mutants", x = "Amp_ex", y = "Fitness (EGvar)") +
  xlim(0, 2200) +
  # ylim(0, 1)+
  # log-scale x axis
  # scale_x_log10() +
  theme(legend.position = "none") # hide legend for clarity
p=p+geom_line(data=data.frame(x = x), aes(x = x, y = fitness_km_vmax_x(x,p2_try,p4_try)), color = "black", size=1)
print(p)

# p2_start=mean(kinetics_df$P_Km,na.rm=TRUE)
# p4_start=5000
p2_sample <- seq(-6, 6, length.out=104)
p4_sample <- seq(-6, 6, length.out=104)
# param_grid <- expand.grid(p2=p2_range, p4=p4_range)
cal_parameter <- function(x,y,p2_range=p2_sample,p4_range=p4_sample){
  param_grid <- expand.grid(p2=p2_range, p4=p4_range)
  results=data.frame(p2=param_grid$p2, p4=param_grid$p4, R_squared=NA,
                    p_p2=NA, p_p4=NA,p_normal=NA,cv_p2=NA,cv_p4=NA,
                    Km_pred=NA, Vmax_D_pred=NA)
  for (i in 1:nrow(param_grid)) {
    # if (i %% 1000 == 0) {
    #   cat("Fitting parameter set", i, "of", nrow(param_grid), "\n")
    # }
    p2_start=param_grid$p2[i]
    p4_start=param_grid$p4[i]
    try({
      fit <- nlsLM(
        y ~ 1/(1+( (x - 10^p2 - 10^p4 + sqrt((10^p2 + 10^p4 - x)^2 + 4 * 10^p2 * x)) / 2 / E_fit )^k_fit) - 1/(1+(x/E_fit)^k_fit),
        # y  ~ 1/(1+(x-p4+sqrt((p4-x)^2+4*p2*x)/2*E_fit)^k_fit),
        # y ~ p2*x+p4,
        # E_fit and k_fit are from previous fitting
        # p2 is km, p4 = p2 + vmax/D, D is diffusion rate, here we fit p4 directly
        start = list(p2 = p2_start, p4=p4_start),
        control = nls.lm.control(maxiter = 1024),
        # set limit for p2 and p4
        #lower = c(0, 0),
        #upper = c(1000, Inf)
      )
      # ## filtering criteria
      # p_critical_para=0.10
      # p_critical_norm=0.10
      # std_ratio=1
      # # if the parameters are not significant, skip
      # if (summary(fit)$coefficients[,"Pr(>|t|)"][1] > p_critical_para | summary(fit)$coefficients[,"Pr(>|t|)"][2] > p_critical_para) {
      #   next
      # }
      # # if parameters standard error is too large, skip
      # if (summary(fit)$coefficients[,"Std. Error"][1] > std_ratio*abs(coef(fit)[1]) | summary(fit)$coefficients[,"Std. Error"][2] > std_ratio*abs(coef(fit)[2])) {
      #   next
      # }
      # # if the distribution of residuals is not normal, skip
      # if (shapiro.test(residuals(fit))$p.value < p_critical_norm) {
      #   next
      # }
      rss <- sum(residuals(fit)^2);tss <- sum((y - mean(y))^2);
      rsq <- 1 - rss/tss
      coef_fit <- coef(fit)
      p2_fit <- coef_fit["p2"]; p4_fit <- coef_fit["p4"]
      std_fit <- summary(fit)$coefficients[,"Std. Error"]
      p_fit <- summary(fit)$coefficients[,"Pr(>|t|)"]
      std_p2 <- std_fit["p2"]; std_p4 <- std_fit["p4"]
      cv_p2 <- std_p2/abs(p2_fit); cv_p4 <- std_p4/abs(p4_fit)
      p_p2 <- p_fit["p2"]; p_p4 <- p_fit["p4"]
      p_normal <- shapiro.test(residuals(fit))$p.value
      results[i,c("R_squared","p_p2","p_p4","p_normal","cv_p2","cv_p4","Km_pred","Vmax_D_pred")]=c(rsq,p_p2,p_p4,p_normal,cv_p2,cv_p4,10^p2_fit,10^p4_fit)
      # results[i,c("R_squared","Km_pred","Vmax_D_pred")]=c(rsq,p2_fit,p4_fit-p2_fit)
    }, silent=TRUE)
  }
  return(results)
  # fitness_km_vmax_x <- function(x) {
  #   x1=(x-p4_fit+sqrt((p4_fit-x)^2+4*p2_fit*x))/2
  #   y=1/(1+(x1/E_fit)^k_fit)
  #   return(y)
  #   # p2_fit*x+p4_fit
  # }
  # ggplot(mutant_list[[1]], aes(x = Amp_ex, y = EGvar)) +
  #   geom_point() +
  #   stat_function(fun = fitness_km_vmax_x, color = "red", size = 1) +
  #   theme_classic() +
  #   labs(title = paste("Fitness vs Amp_ex | R-squared: ", round(rsq, 3)), x = "Amp_ex", y = "Fitness (EGvar)") +
  #   xlim(0, 2100) +
  #   ylim(0, 1)+
  #   scale_x_log10()
}
# test for one mutant
names(mutant_list)
x=mutant_list[[45]]$Amp_ex
y=mutant_list[[45]]$EGvar_relative
# plot(x, y)
# lines(x, diff_fitness_x(x, p2 = p2_try, p4 = p4_try), col = "red")
results=cal_parameter(x,y)
results=na.omit(results);
# filter results
p_critical_para=0.25
r_squared_critical=0.7
cv_critical=5

# see the histogram of R-squared, Km_pred, Vmax_D_pred on one plot
ggplot(results, aes(x = Km_pred)) +
  xlim(0, 1000) +
  geom_histogram(binwidth = 10, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Km_predict", x = "R-squared", y = "Count") +
  theme_classic()

p_critical_para_try=0.25
r_squared_critical_try=0.5
cv_critical_try=5
min_count_try=100
# final parameters: most frequent Km_pred and Vmax_D_pred
decide_final_params <- function(results,p_critical_para=p_critical_para_try,
                                r_squared_critical=r_squared_critical_try,
                                cv_critical=cv_critical_try,
                                min_count=min_count_try){
  results=na.omit(results)
  results <- results %>%
    filter(R_squared > r_squared_critical & p_p2 < p_critical_para & p_p4 < p_critical_para & cv_p2 < cv_critical & cv_p4 < cv_critical & Km_pred < 1000)
  if (nrow(results) == 0) {
    return(data.frame(Km_pred=NA, Vmax_D_pred=NA, count=0))
  }
  # round to 1 decimal place
  results$Km_pred <- round(results$Km_pred, 1)
  results$Vmax_D_pred <- round(results$Vmax_D_pred, 1)
  final_params <- results %>%
    group_by(Km_pred, Vmax_D_pred) %>%
    summarise(count = n()) %>%
    arrange(desc(count)) %>%
    slice(1)
  if (final_params$count[1] < min_count) {
    return(data.frame(Km_pred=NA, Vmax_D_pred=NA, count=0))
  }
  if (nrow(final_params) > 1) {
    final_params <- final_params[1, ]
  }
  return(final_params)
}

temp=decide_final_params(results)
print(temp)


## loop over all mutants
mutant_parameters <- data.frame(Mutant=names(mutant_list),
                                Km_pred=NA,
                                Vmax_D_pred=NA,
                                # R_squared=NA,
                                count=NA,
                                stringsAsFactors=FALSE)
for (i in 1:length(mutant_list)) {
  mutant_name=names(mutant_list)[i]
  cat("Processing mutant", mutant_name, "(", i, "of", length(mutant_list), ")\n")
  df=mutant_list[[i]]
  x=df$Amp_ex
  y=df$EGvar_relative
  results=cal_parameter(x,y)
  params=decide_final_params(results)
  mutant_parameters[i,"Km_pred"]=params$Km_pred
  mutant_parameters[i,"Vmax_D_pred"]=params$Vmax_D_pred
  mutant_parameters[i,"count"]=params$count
}

# merge with kinetics_df
check_df <- merge(kinetics_df, mutant_parameters, by = "Mutant")
check_df <- na.omit(check_df)
# plot Km_pred vs Km
ggplot(check_df, aes(x = P_Km, y = Km_pred)) +
  geom_point() +
  # geom_smooth(method = "lm", se = FALSE, color = "blue") +
  # add y=x line
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  # point out WT
  geom_point(data = subset(check_df, Mutant == "WT"), aes(x = P_Km, y = Km_pred), color = "red", size = 2) +
  annotate("text", x = subset(check_df, Mutant == "WT")$P_Km, y = subset(check_df, Mutant == "WT")$Km_pred, label = "WT", vjust = -1, color = "red") +
  # label with mutant names
  geom_text(aes(label = Mutant), hjust = 0, vjust = 1.5, size = 3) +
  coord_fixed(ratio = 1) +
  labs(title = "Predicted Km vs Measured Km", x = "Measured Km", y = "Predicted Km") +
  xlim(0, 1.5*max(check_df$Km, check_df$Km_pred)) +
  ylim(0, 1.5*max(check_df$Km, check_df$Km_pred))+
  # log-scale both axes
  # scale_x_log10() +
  # scale_y_log10()+
  theme_classic()











