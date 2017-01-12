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

target <- get_normal()

# number of observations
nobservations <- 10000
# parameter of data-generating process
true_theta <- c(0.4, 0.9)
obs <- target$robservation(nobservations, true_theta,
                                  target$parameters, target$generate_randomness(nobservations))

save(true_theta, obs, file = "~/Dropbox/ABCD/Results/data/normaldata.RData")

hist(obs)
mean(obs)
sd(obs)
