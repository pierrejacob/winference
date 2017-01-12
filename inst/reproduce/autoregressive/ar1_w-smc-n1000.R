#+ presets, echo = FALSE, warning = FALSE, message = FALSE
library(winference)
library(doMC)
library(doRNG)
library(dplyr)
library(ggthemes)

registerDoMC(cores = 6)

rm(list = ls())
set.seed(11)
setmytheme()
#
# model
target <- get_autoregressive()
#
# number of observations
nobservations <- 1000
load(file = "~/Dropbox/ABCD/Results/data/ar1data.RData")
obs <- obs[1:nobservations]
obs_sorted <- sort(obs)

compute_d <- function(theta, metric = metricL1){
  fake_rand <- target$generate_randomness(nobservations)
  fake_obs <- target$robservation(nobservations, theta, target$parameters, fake_rand)
  fake_ordered_obs <- sort(fake_obs)
  return(metric(obs_sorted, fake_ordered_obs))
}


proposal <- mixture_proposal()

param_algo <- list(nthetas = 1024, nmoves = 5, proposal = proposal,
                   nsteps = 10, minimum_diversity = 0.5, R = 2, maxtrials = 1000)

filename <- paste0("~/Dropbox/ABCD/Results/autoregressive/ar1data.n",
                   nobservations, ".wsmc.RData")
# results <- wsmc(compute_d, target, param_algo, savefile = filename)
# wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
# nsteps <- max(wsmc.df$step)
# save(wsmc.df, results, nsteps, file = filename)
load(filename)
wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps <- max(wsmc.df$step)

g <- ggplot(wsmc.df, aes(x = rho, group = step)) + geom_density(aes(y = ..density..), colour = "darkgrey")
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(rho))
g

g <- ggplot(wsmc.df, aes(x = logsigma, group = step)) + geom_density(aes(y = ..density..), colour = "darkgrey")
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(log(sigma)))
g

g <- ggplot(wsmc.df, aes(x = rho, y = logsigma, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(rho)) + ylab(expression(log(sigma)))
g

