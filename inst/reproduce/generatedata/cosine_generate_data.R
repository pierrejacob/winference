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

target <- get_cosine()

# number of observations
nobservations <- 10000
# parameter of data-generating process
true_theta <- c(1/80, pi/4, 0, log(2))
obs <- target$robservation(nobservations, true_theta,
                           target$parameters, target$generate_randomness(nobservations))

save(true_theta, obs, file = "~/Dropbox/ABCD/Results/data/cosinedata.RData")

hist(obs)
plot(obs[1:1000], type = "l")
acf(obs, 100)
mean(obs)
sd(obs)
# approx equal to exp(0.9)/sqrt(1-0.7^2)
