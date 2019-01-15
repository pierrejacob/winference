library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()

set.seed(11)

target <- get_levydriven()

prefix <- ""

# number of observations
nobservations <- 50000
# parameter of data-generating process
true_theta <- c(0, 0, 0.5, 0.0625, 0.01)
obs <- target$robservation(nobservations, true_theta,
                           target$parameters, target$generate_randomness(nobservations))

save(true_theta, obs, file = paste0(prefix,"levydrivendata.RData"))

# hist(obs)
plot(obs, type = "l")
# plot(obs[1:250], type = "l")
# acf(obs)
# mean(obs)
# sd(obs)
# approx equal to exp(0.9)/sqrt(1-0.7^2)

library(microbenchmark)
microbenchmark(
obs <- target$robservation(1e4, true_theta,
                           target$parameters, target$generate_randomness(1e4)),
times = 100)


