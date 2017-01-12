#+ presets, echo = FALSE, warning = FALSE, message = FALSE
library(winference)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(foreach)
library(doMC)
library(doRNG)
library(reshape2)
library(tidyr)
registerDoMC(cores = 6)
rm(list = ls())
setmytheme()
set.seed(11)
target <- get_toggleswitch()
# number of observations
nobservations <- 2000
load(file = "~/Dropbox/ABCD/Results/data/toggleswitchdata.RData")
obs <- obs[1:nobservations]
obs_sorted <- sort(obs)

# function to compute distance between observed data and data generated given theta
compute_d <- function(theta, metric = metricL2){
  fake_rand <- target$generate_randomness(nobservations)
  fake_obs <- target$robservation(nobservations, theta, target$parameters, fake_rand)
  fake_obs_sorted <- sort(fake_obs)
  return(metric(obs_sorted, fake_obs_sorted))
}
# compute_d(true_theta)

proposal <- mixture_proposal()

param_algo <- list(nthetas = 2048, nmoves = 1, proposal = proposal,
                   nsteps = 50, minimum_diversity = 0.5, R = 2, maxtrials = 1000)

filename <- paste0("~/Dropbox/ABCD/Results/toggleswitch/toggleswitchdata.n",
                   nobservations, ".wsmc_rhit.RData")

# results <- wsmc_rhit(compute_d, target, param_algo, savefile = filename)
# wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
# nsteps <- max(wsmc.df$step)
# save(wsmc.df, results, nsteps, file = filename)
#
load(filename)
wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps <- max(wsmc.df$step)
target$parameter_names

qplot(x = 1:nsteps, y = results$threshold_history, geom ="line") + scale_y_log10()

g <- ggplot(wsmc.df, aes(x = alpha_1, y = alpha_2, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(alpha[1])) + ylab(expression(alpha[2]))
g + geom_point(data=NULL, aes(x = true_theta[1], y = true_theta[2], colour = NULL, group = NULL),
               size = 5)

g <- ggplot(wsmc.df, aes(x = alpha_1, group = step))
g <- g + geom_density(aes(y = ..density..))
g
g <- ggplot(wsmc.df, aes(x = alpha_2, group = step))
g <- g + geom_density(aes(y = ..density..))
g

g <- ggplot(wsmc.df, aes(x = beta_1, y = beta_2, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(beta[1])) + ylab(expression(beta[2]))
g + geom_point(data=NULL, aes(x = true_theta[3], y = true_theta[4], colour = NULL, group = NULL),
               size = 5)

g <- ggplot(wsmc.df, aes(x = mu, y = sigma, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(mu)) + ylab(expression(sigma))
g + geom_point(data=NULL, aes(x = true_theta[5], y = true_theta[6], colour = NULL, group = NULL),
               size = 5)

g <- ggplot(wsmc.df, aes(x = sigma, y = gamma, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(sigma)) + ylab(expression(gamma))
g + geom_point(data=NULL, aes(x = true_theta[6], y = true_theta[7], colour = NULL, group = NULL),
               size = 5)

g <- ggplot(wsmc.df %>% filter(step > 20), aes(x = gamma, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none")
g <- g + xlab(expression(gamma))
g
