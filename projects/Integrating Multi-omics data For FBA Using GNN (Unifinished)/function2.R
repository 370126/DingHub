Bayesian = function(comb){

library(rjags)
library(coda)

model_code <- "
model {
  for (i in 1:N) {
    z[i] ~ dbern(Prob[i])
    y[i] ~ dnorm(mu[z[i]+1], tau[z[i]+1])
  }
  mu[2] <- 0.8
  tau[2] <- 25
  mu[1] <- 0.2
  tau[1] <- 100
}
"

data_list <- list(
  N = nrow(comb),
  y = comb$abundence,
  Prob = comb$Prob
)

model <- jags.model(
  textConnection(model_code),
  data = data_list,
  n.chains = 3
)

update(model, 1000)
samples <- coda.samples(
  model,
  variable.names = "z",  # 不再需要pi
  n.iter = 1000
)

z_post <- summary(samples)$statistics[grep("^z\\[", rownames(summary(samples)$statistics)), "Mean"]

return(z_post)
}

