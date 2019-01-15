library(winference)
rm(list = ls())
setmytheme()
set.seed(11)
prefix = ""
target <- get_autoregressive()
# number of observations
nobservations <- 10000
# parameter of data-generating process
true_theta <- c(0.7, 0.9)
obs <- target$robservation(nobservations, true_theta,
                           target$parameters, target$generate_randomness(nobservations))

save(true_theta, obs, file = paste0(prefix, "ar1data.RData"))
# plot the observations
plot(obs, type = "l")
acf(obs)
