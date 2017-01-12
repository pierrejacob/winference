#+ presets, echo = FALSE, warning = FALSE, message = FALSE
library(winference)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(foreach)
library(doMC)
library(doRNG)
library(reshape2)
registerDoMC(cores = 10)
rm(list = ls())
setmytheme()

set.seed(11)

target <- get_toggleswitch()

# number of observations
nobservations <- 10000
# parameter of data-generating process
# true_theta <- c(alpha_1, alpha_2, beta_1, beta_2, mu, sigma, gamma)
true_theta <- c(22, 12, 4, 4.5, 325, 0.25, 0.15)

obs <- target$robservation(nobservations, true_theta,
                           target$parameters, target$generate_randomness(nobservations))
save(true_theta, obs, file = "~/Dropbox/ABCD/Results/data/toggleswitchdata.RData")

hist(obs)
plot(obs[1:1000], type = "l")
mean(obs)
sd(obs)
