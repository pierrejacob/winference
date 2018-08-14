#+ presets, echo = FALSE, warning = FALSE, message = FALSE
library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()

set.seed(11)

prefix = ""

target <- get_gandk()

# number of observations
nobservations <- 10000
# parameter of data-generating process
true_theta <- c(3, 1, 2, 0.5)
obs <- target$robservation(nobservations, true_theta,
                           target$parameters, target$generate_randomness(nobservations))

save(true_theta, obs, file = paste0(prefix,"gandkdata.Rdata"))

hist(obs)
plot(obs[1:1000], type = "l")
mean(obs)
sd(obs)

