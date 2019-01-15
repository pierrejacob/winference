library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()
set.seed(11)

target <- get_cosine()

prefix = ""

# number of observations
nobservations <- 10000
# parameter of data-generating process
true_theta <- c(1/80, pi/4, 0, log(2))
obs <- target$robservation(nobservations, true_theta,
                           target$parameters, target$generate_randomness(nobservations))

save(true_theta, obs, file = paste0(prefix,"cosinedata.RData"))

plot(obs[1:100], type = "l")
# hist(obs)
# acf(obs, 100)
# mean(obs)
# sd(obs)
