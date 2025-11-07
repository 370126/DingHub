## load tsv files
fitness_amp_ex_df <- read.csv("./fitness_amp_ex_df.csv")
# bad_mutant=c("A5", "C10", "D9", "D10", "D11","B1","D1")
bad_amp_ex=c(256)
# fitness_amp_ex_df <- fitness_amp_ex_df[!fitness_amp_ex_df$Mutant %in% bad_mutant, ]
fitness_amp_ex_df <- fitness_amp_ex_df[!fitness_amp_ex_df$Amp_ex %in% bad_amp_ex, ]

kinetics_df <- read.csv("./35C_protein_Km_Vmax.tsv", sep = "\t")
kinetics_df$P_vm_km <- kinetics_df$P_Vmax / kinetics_df$P_Km
# log
kinetics_df$log_P_Km <- log10(kinetics_df$P_Km)
kinetics_df$log_P_Vmax <- log10(kinetics_df$P_Vmax)

library(ggplot2)
library(tidyverse)
library(minpack.lm)

# split data into two parts
# sort by P_vm_km, then split into two halves, training is the odd rows, testing is the even rows
kinetics_df <- kinetics_df[order(kinetics_df$P_vm_km), ]
kinetics_train_df <- kinetics_df[seq(1, nrow(kinetics_df), by = 2), ]
kinetics_test_df <- kinetics_df[seq(2, nrow(kinetics_df), by = 2), ]
fitness_amp_train_df <- subset(fitness_amp_ex_df, Mutant %in% kinetics_train_df$Mutant)
fitness_amp_test_df <- subset(fitness_amp_ex_df, Mutant %in% kinetics_test_df$Mutant)


amp_ex_amp_p_x <- function(x,km,vmax,D) {
  # x is Amp_ex, x1 is Amp_p
  x1=(x-km-vmax/D+sqrt((km+vmax/D-x)^2+4*km*x))/2
  return(x1)
}

fitness_relative_amp_ex_x <- function(x,km,vmax,D,E,k) {
  # x is Amp_ex, x1 is Amp_p
  # y is fitness relative to NC
  # for NC, Amp_ex=Amp_p
  x1=(x-km-vmax/D+sqrt((km+vmax/D-x)^2+4*km*x))/2
  y=1/(1+(x1/E)^k)-1/(1+(x/E)^k)
  return(y)
}
# fit D,E,k using training data
fitness_amp_km_vmax_df=merge(fitness_amp_train_df, kinetics_train_df[,c("Mutant","P_Km","P_Vmax")], by="Mutant")
# with known Km and Vmax for each mutant, fit D,E,k
fit <- nlsLM(EGvar_relative ~ fitness_relative_amp_ex_x(Amp_ex, P_Km, P_Vmax, D, E, k),
             data = fitness_amp_km_vmax_df,
             start = list(D = 10, E = 10, k = 1),
             control = nls.lm.control(maxiter = 1000))
# check fitting summary
summary(fit)
# calculate r-squared
fitness_amp_km_vmax_df$Pred <- predict(fit)
SSE <- sum((fitness_amp_km_vmax_df$EGvar_relative - fitness_amp_km_vmax_df$Pred)^2)
SST <- sum((fitness_amp_km_vmax_df$EGvar_relative - mean(fitness_amp_km_vmax_df$EGvar_relative))^2)
R2 <- 1 - SSE/SST
print(paste("R-squared: ", R2))

# get fitted parameters
params <- coef(fit)
D_fit <- params["D"]
E_fit <- params["E"]
k_fit <- params["k"]
# plot fitting curve for each mutant separately on same plot
curve_df <- fitness_amp_km_vmax_df %>%
  group_by(Mutant) %>%
  mutate(Pred = fitness_relative_amp_ex_x(Amp_ex, P_Km, P_Vmax,
                                          D_fit, E_fit, k_fit))
ggplot(curve_df, aes(x = Amp_ex, color = Mutant)) +
  geom_point(aes(y = EGvar_relative)) +
  geom_line(aes(y = Pred)) +
  scale_color_viridis_d() +
  scale_x_log10() +
  theme_classic() +
  labs(title = "Training Data: Fitting Curve for Each Mutant",
       x = "Amp_ex",
       y = "EGvar_relative to NC") +
  theme(legend.position = "right")




#### now fit vmax & km for test data
mutants_names=unique(fitness_amp_test_df$Mutant)


km_sample <- seq(-6, 6, length.out=104)
vmax_sample <- seq(-6, 6, length.out=104)
cal_parameter_split <- function(x,y,km_range=km_sample,vmax_range=vmax_sample){
  param_grid <- expand.grid(km=km_range, vmax=vmax_range)
  results=data.frame(km=param_grid$km, vmax=param_grid$vmax, R_squared=NA,
                    p_km=NA, p_vmax=NA,p_normal=NA,cv_km=NA,cv_vmax=NA,
                    Km_pred=NA, Vmax_pred=NA)
  for (i in 1:nrow(param_grid)) {
    # if (i %% 1000 == 0) {
    #   cat("Fitting parameter set", i, "of", nrow(param_grid), "\n")
    # }
    km_start=param_grid$km[i]
    vmax_start=param_grid$vmax[i]
    try({
      fit <- nlsLM(
        y ~ 1/(1+( (x - 10^km - 10^vmax/D_fit + sqrt((10^km + 10^vmax/D_fit - x)^2 + 4 * 10^km * x)) / 2 / E_fit )^k_fit) - 1/(1+(x/E_fit)^k_fit),
        # E_fit and k_fit are from previous fitting
        # km is km, vmax = km + vmax/D, D is diffusion rate, here we fit vmax directly
        start = list(km = km_start, vmax=vmax_start),
        control = nls.lm.control(maxiter = 1024),
      )
      rss <- sum(residuals(fit)^2);tss <- sum((y - mean(y))^2);
      rsq <- 1 - rss/tss
      coef_fit <- coef(fit)
      km_fit <- coef_fit["km"]; vmax_fit <- coef_fit["vmax"]
      std_fit <- summary(fit)$coefficients[,"Std. Error"]
      p_fit <- summary(fit)$coefficients[,"Pr(>|t|)"]
      std_km <- std_fit["km"]; std_vmax <- std_fit["vmax"]
      cv_km <- std_km/abs(km_fit); cv_vmax <- std_vmax/abs(vmax_fit)
      p_km <- p_fit["km"]; p_vmax <- p_fit["vmax"]
      p_normal <- shapiro.test(residuals(fit))$p.value
      results[i,c("R_squared","p_km","p_vmax","p_normal","cv_km","cv_vmax","Km_pred","Vmax_pred")]=c(rsq,p_km,p_vmax,p_normal,cv_km,cv_vmax,10^km_fit,10^vmax_fit)
      # results[i,c("R_squared","Km_pred","Vmax_pred")]=c(rsq,km_fit,vmax_fit-km_fit)
    }, silent=TRUE)
  }
  return(results)
}
# test for one mutant
x=subset(fitness_amp_test_df, Mutant=="B6")$Amp_ex
y=subset(fitness_amp_test_df, Mutant=="B6")$EGvar_relative
plot(x, y)
results=cal_parameter_split(x,y)
results=na.omit(results);

# filter results
p_critical_para_try=0.25
r_squared_critical_try=0.5
cv_critical_try=5
min_count_try=100
decide_final_params_split <- function(results,p_critical_para=p_critical_para_try,
                                r_squared_critical=r_squared_critical_try,
                                cv_critical=cv_critical_try,
                                min_count=min_count_try){
  results=na.omit(results)
  results <- results %>%
    filter(R_squared > r_squared_critical & p_km < p_critical_para & p_vmax < p_critical_para & cv_km < cv_critical & cv_vmax < cv_critical & Km_pred < 1000)
  if (nrow(results) == 0) {
    return(data.frame(Km_pred=NA, Vmax_pred=NA, count=0))
  }
  # round to 1 decimal place
  results$Km_pred_round <- round(results$Km_pred, 1)
  results$Vmax_pred_round <- round(results$Vmax_pred, 1)
  # find the most frequent parameter set by round values, and return the median of original values
  final_params <- results %>%
    group_by(Km_pred_round, Vmax_pred_round) %>%
    summarise(count = n()) %>%
    arrange(desc(count)) %>%
    slice(1)
  if (final_params$count[1] < min_count) {
    return(data.frame(Km_pred=NA, Vmax_pred=NA, count=0))
  }
  if (nrow(final_params) > 1) {
    final_params <- final_params[1, ]
  }
  median_Km <- median(results$Km_pred[results$Km_pred_round == final_params$Km_pred_round])
  median_Vmax <- median(results$Vmax_pred[results$Vmax_pred_round == final_params$Vmax_pred_round])
  final_params$Km_pred <- median_Km
  final_params$Vmax_pred <- median_Vmax
  return(final_params)
}
# decide_final_params_split(results)



## loop over all mutants
mutant_parameters <- data.frame(Mutant=mutants_names,
                                Km_pred=NA,
                                Vmax_pred=NA,
                                # R_squared=NA,
                                count=NA,
                                stringsAsFactors=FALSE)
for (i in 1:length(mutants_names)) {
  mutant_name=mutants_names[i]
  cat("Processing mutant", mutant_name, "(", i, "of", length(mutants_names), ")\n")
  df=subset(fitness_amp_test_df, Mutant==mutant_name)
  x=df$Amp_ex
  y=df$EGvar_relative
  results=cal_parameter_split(x,y)
  params=decide_final_params_split(results)
  mutant_parameters[i,"Km_pred"]=params$Km_pred
  mutant_parameters[i,"Vmax_pred"]=params$Vmax_pred
  mutant_parameters[i,"count"]=params$count
}

mutant_parameters$log_Km_pred <- log10(mutant_parameters$Km_pred)
mutant_parameters$log_Vmax_pred <- log10(mutant_parameters$Vmax_pred)

# merge with kinetics_df
check_df <- merge(kinetics_test_df , mutant_parameters, by = "Mutant")
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
  # coord_fixed(ratio = 1) +
  labs(title = "Predicted Km vs Measured Km", x = "Measured Km", y = "Predicted Km") +
  #xlim(0, 1.5*max(check_df$Km, check_df$Km_pred)) +
  #ylim(0, 1.5*max(check_df$Km, check_df$Km_pred))+
  # log-scale both axes
  # scale_x_log10() +
  # scale_y_log10()+
  theme_classic()

# plot Vmax_pred vs Vmax
ggplot(check_df, aes(x = P_Vmax, y = Vmax_pred)) +
  geom_point() +
  # geom_smooth(method = "lm", se = FALSE, color = "blue") +
  # add y=x line
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  # point out WT
  geom_point(data = subset(check_df, Mutant == "WT"), aes(x = P_Vmax, y = Vmax_pred), color = "red", size = 2) +
  annotate("text", x = subset(check_df, Mutant == "WT")$P_Vmax, y = subset(check_df, Mutant == "WT")$Vmax_pred, label = "WT", vjust = -1, color = "red") +
  # label with mutant names
  geom_text(aes(label = Mutant), hjust = 0, vjust = 1.5, size = 3) +
  # coord_fixed(ratio = 1) +
  labs(title = "Predicted Vmax vs Measured Vmax", x = "Measured Vmax", y = "Predicted Vmax") +
  #xlim(0, 1.5*max(check_df$Vmax, check_df$Vmax_pred)) +
  #ylim(0, 1.5*max(check_df$Vmax, check_df$Vmax_pred))+
  # log-scale both axes
  # scale_x_log10() +
  # scale_y_log10()+
  theme_classic()

# plot log_Km_pred vs log_Km
ggplot(check_df, aes(x = log_P_Km, y = log_Km_pred)) +
  geom_point() +
  # geom_smooth(method = "lm", se = FALSE, color = "blue") +
  # add y=x line
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  # point out WT
  geom_point(data = subset(check_df, Mutant == "WT"), aes(x = log_P_Km, y = log_Km_pred), color = "red", size = 2) +
  annotate("text", x = subset(check_df, Mutant == "WT")$log_P_Km, y = subset(check_df, Mutant == "WT")$log_Km_pred, label = "WT", vjust = -1, color = "red") +
  # label with mutant names
  geom_text(aes(label = Mutant), hjust = 0, vjust = 1.5, size = 3) +
  # coord_fixed(ratio = 1) +
  labs(title = "Predicted log(Km) vs Measured log(Km)", x = "Measured log(Km)", y = "Predicted log(Km)") +
  xlim(min(check_df$log_P_Km, check_df$log_Km_pred), max(check_df$log_P_Km, check_df$log_Km_pred)) +
  ylim(min(check_df$log_P_Km, check_df$log_Km_pred), max(check_df$log_P_Km, check_df$log_Km_pred))+
  theme_classic()

# plot log_Vmax_pred vs log_Vmax
ggplot(check_df, aes(x = log_P_Vmax, y = log_Vmax_pred)) +
  geom_point() +
  # geom_smooth(method = "lm", se = FALSE, color = "blue") +
  # add y=x line
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  # point out WT
  geom_point(data = subset(check_df, Mutant == "WT"), aes(x = log_P_Vmax, y = log_Vmax_pred), color = "red", size = 2) +
  annotate("text", x = subset(check_df, Mutant == "WT")$log_P_Vmax, y = subset(check_df, Mutant == "WT")$log_Vmax_pred, label = "WT", vjust = -1, color = "red") +
  # label with mutant names
  geom_text(aes(label = Mutant), hjust = 0, vjust = 1.5, size = 3) +
  # coord_fixed(ratio = 1) +
  labs(title = "Predicted log(Vmax) vs Measured log(Vmax)", x = "Measured log(Vmax)", y = "Predicted log(Vmax)") +
  xlim(min(check_df$log_P_Vmax, check_df$log_Vmax_pred), max(check_df$log_P_Vmax, check_df$log_Vmax_pred)) +
  ylim(min(check_df$log_P_Vmax, check_df$log_Vmax_pred), max(check_df$log_P_Vmax, check_df$log_Vmax_pred))+
  theme_classic()




