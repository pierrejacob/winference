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

target <- get_garch()

# number of observations
nobservations <- 10000
# parameter of data-generating process
true_theta <- c(0.1, 0.15, 0.6)
nobservations <- 50000
obs <- target$robservation(nobservations, true_theta, list(), list())

save(true_theta, obs, file = "~/Dropbox/ABCD/Results/data/garchdata.RData")

hist(obs)
# plot(obs, type = "l")
plot(obs[1:250], type = "l")
acf(obs)
mean(obs)
sd(obs)
# approx equal to exp(0.9)/sqrt(1-0.7^2)
