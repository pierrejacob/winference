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
nobservations <- 10
load(file = "~/Dropbox/ABCD/Results/data/normaldata.RData")
obs <- obs[1:nobservations]
obs_sorted <- sort(obs)

# function to compute distance between observed data and data generated given theta
compute_d <- function(theta, metric = metricL2){
  fake_rand <- target$generate_randomness(nobservations)
  fake_obs <- target$robservation(nobservations, theta, target$parameters, fake_rand)
  fake_obs_sorted <- sort(fake_obs)
  return(metric(obs_sorted, fake_obs_sorted))
}

######################
proposal <- independent_proposal()

param_algo <- list(nthetas = 1024, nmoves = 25, proposal = proposal,
                   nsteps = 10, minimum_diversity = 0.5)

filename <- paste0("~/Dropbox/ABCD/Results/univariatenormal/normaldata.n",
                   nobservations, ".wsmc.RData")
wsmcresults <- wsmc(compute_d, target, param_algo)
wsmc.df <- wsmc_to_dataframe(wsmcresults, target$parameter_names)
nsteps <- max(wsmc.df$step)
save(wsmc.df, wsmcresults, nsteps, file = filename)

names(wsmcresults)
gthreshold <- qplot(x = 1:(length(wsmcresults$threshold_history)), y = wsmcresults$threshold_history, geom = "line")
gthreshold <- gthreshold + xlab("step") + ylab("threshold")
gthreshold
tail(wsmcresults$threshold_history)


head(wsmc.df)

g <- ggplot(wsmc.df, aes(x = mu, y = sigma, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
g

g <- ggplot(wsmc.df, aes(x = mu, colour = factor(step), group = step)) + geom_density(aes(y = ..density..))
g <- g + geom_vline(xintercept = true_theta[1]) + theme(legend.position = "none")
g

# wsmcresults_2 <- wsmc_fixedthresholds(wsmcresults$threshold_history, compute_d, target, param_algo)
# df <- wsmc_to_dataframe(wsmcresults_2, target$parameter_names)
# nsteps <- max(df$step)
# g <- ggplot(df, aes(x = mu, colour = factor(step), group = step)) + geom_density(aes(y = ..density..))
# g <- g + geom_vline(xintercept = true_theta[1])  + theme(legend.position = "none")
# g
#
# thresholds_ <- c(0.5, 0.4, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.05, 0.05, 0.05, 0.05, 0.05)
# wsmcresults_3 <- wsmc_fixedthresholds(thresholds_, compute_d, target, param_algo)
# df <- wsmc_to_dataframe(wsmcresults_3, target$parameter_names)
# nsteps <- max(df$step)
# g <- ggplot(df, aes(x = mu, colour = factor(step), group = step)) + geom_density(aes(y = ..density..))
# g <- g + geom_vline(xintercept = true_theta[1])
# g
#
