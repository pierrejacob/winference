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

target <- get_gandk()

# number of observations
nobservations <- 10000
# parameter of data-generating process
true_theta <- c(3, 1, 2, 0.5)
obs <- target$robservation(nobservations, true_theta,
                           target$parameters, target$generate_randomness(nobservations))

save(true_theta, obs, file = "~/Dropbox/ABCD/Results/data/gandkdata.RData")

hist(obs)
plot(obs[1:1000], type = "l")
mean(obs)
sd(obs)
