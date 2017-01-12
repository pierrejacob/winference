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
registerDoMC(cores = 2)
rm(list = ls())
setmytheme()
set.seed(11)
target <- get_gandk()
# number of observations
nobservations <- 1000
load(file = "~/Dropbox/ABCD/Results/data/gandkdata.RData")
obs <- obs[1:nobservations]
obs_sorted <- sort(obs)

# function to compute distance between observed data and data generated given theta
compute_d <- function(theta, metric = metricL2){
  fake_rand <- target$generate_randomness(nobservations)
  fake_obs <- target$robservation(nobservations, theta, target$parameters, fake_rand)
  fake_obs_sorted <- sort(fake_obs)
  return(metric(obs_sorted, fake_obs_sorted))
}


proposal <- mixture_proposal()

param_algo <- list(nthetas = 1024, nmoves = 5, proposal = proposal,
                   nsteps = 25, minimum_diversity = 0.5, R = 2, maxtrials = 1000)

filename <- paste0("~/Dropbox/ABCD/Results/gandk/gandkdata.n",
                   nobservations, ".wsmc_rhit.RData")
# results <- wsmc_rhit(compute_d, target, param_algo, savefile = filename)
# wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
# nsteps <- max(wsmc.df$step)
# save(wsmc.df, results, nsteps, file = filename)
#
load(filename)
wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps <- max(wsmc.df$step)

g <- ggplot(wsmc.df, aes(x = A, y = B, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(A)) + ylab(expression(B))
g

g <- ggplot(wsmc.df, aes(x = g, y = k, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(g)) + ylab(expression(k))
g


g <- ggplot(wsmc.df %>% filter(step > 20), aes(x = g, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none")
g <- g + xlab(expression(g))
g
