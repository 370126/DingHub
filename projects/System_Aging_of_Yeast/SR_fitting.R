# fitting functions for the relationship between steepness and lifespan_rela under different parameters



# setwd("C:/Users/20145/Desktop/westlake_summer/")

simu_long=read.csv("./data/simu_long.csv")

### fit simu_long data

x=simu_long$lifespan_rela[simu_long$para=="eta"]

y=simu_long$steepness_rela[simu_long$para=="eta"]

library(minpack.lm)

fit <- nlsLM(

  y ~ a_fit * exp(-b_fit * x) + c_fit * exp(-d_fit * x),

  start = list(a_fit = max(y), b_fit = 0.1, c_fit = min(y), d_fit = 0.01),

  control = nls.lm.control(maxiter = 10000)

)

summary(fit)

coef_fit <- coef(fit)

a_fit <- coef_fit["a_fit"];b_fit <- coef_fit["b_fit"];c_fit <- coef_fit["c_fit"];d_fit <- coef_fit["d_fit"]

eta_x <- function(x) {

  a_fit * exp(-b_fit * x) + c_fit * exp(-d_fit * x)

}

# plot(x,y)
# 
# lines(x, eta_x(x), col = "red")
# 
# title("eta fitting")



x=simu_long$lifespan_rela[simu_long$para=="epsilon"]

y=simu_long$steepness_rela[simu_long$para=="epsilon"]

fit <- nlsLM(

  y ~ e_fit * exp(f_fit * x) + g_fit * exp(h_fit * x),

  start = list(e_fit = max(y), f_fit = 0.1, g_fit = min(y), h_fit = 1),

  control = nls.lm.control(maxiter = 10000)

)

summary(fit)

coef_fit <- coef(fit)

e_fit <- coef_fit["e_fit"];f_fit <- coef_fit["f_fit"];g_fit <- coef_fit["g_fit"];h_fit <- coef_fit["h_fit"];

epsilon_x <- function(x) {

  e_fit * exp(f_fit * x) + g_fit * exp(h_fit * x)

}

# plot(x,y)
# 
# lines(x, epsilon_x(x), col = "red")
# 
# title("epsilon fitting")



x=simu_long$lifespan_rela[simu_long$para=="X_c"]

y=simu_long$steepness_rela[simu_long$para=="X_c"]

fit <- nlsLM(

  y ~ i_fit * x^j_fit + k_fit,

  start = list(i_fit = max(y), j_fit = 0.1, k_fit = min(y)),

  control = nls.lm.control(maxiter = 10000)

)

summary(fit)

coef_fit <- coef(fit)

i_fit <- coef_fit["i_fit"];j_fit <- coef_fit["j_fit"];k_fit <- coef_fit["k_fit"];

xc_x <- function(x) {

  i_fit * x^j_fit + k_fit

}

# plot(x,y)
# 
# lines(x, xc_x(x), col = "red")
# 
# title("X_c fitting")



rm(fit, coef_fit,x,y)

