library(winference)
registerDoParallel(cores = 6)
rm(list = ls())
setmytheme()

set.seed(11)

target <- get_levydriven()

# number of observations
nobservations <- 50000
# parameter of data-generating process
true_theta <- c(0, 0, 0.5, 0.0625, 0.01)
obs <- target$robservation(nobservations, true_theta,
                           target$parameters, target$generate_randomness(nobservations))

save(true_theta, obs, file = "levydrivendata.RData")

# hist(obs)
plot(obs, type = "l")
# plot(obs[1:250], type = "l")
# acf(obs)
# mean(obs)
# sd(obs)
# approx equal to exp(0.9)/sqrt(1-0.7^2)
