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

# number of observations
nobservations <- 10000
# parameter of data-generating process
obs <- rgamma(nobservations, shape = 10, rate = 5)
save(obs, file = "~/Dropbox/ABCD/Results/data/gammadata.RData")

hist(obs)
mean(obs)
sd(obs)
