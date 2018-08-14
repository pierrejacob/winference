library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()

prefix = ""

set.seed(11)
target <- get_toggleswitch()

# number of observations
nobservations <- 2000
# parameter of data-generating process
# true_theta <- c(alpha_1, alpha_2, beta_1, beta_2, mu, sigma, gamma)
true_theta <- c(22, 12, 4, 4.5, 325, 0.25, 0.15)

obs <- target$robservation(nobservations, true_theta,
                           target$parameters, target$generate_randomness(nobservations))
save(true_theta, obs, file = paste0(prefix,"toggleswitchdata.RData"))

hist(obs)
# plot(obs[1:1000], type = "l")
# mean(obs)
# sd(obs)

