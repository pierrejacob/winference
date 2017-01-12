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
registerDoMC(cores = 5)
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


proposal <- independent_proposal()

param_algo <- list(nthetas = 1024, nmoves = 1, proposal = proposal,
                   nsteps = 10, minimum_diversity = 0.5, R = 2, maxtrials = 1000)

filename <- paste0("~/Dropbox/ABCD/Results/univariatenormal/normaldata.n",
nobservations, ".wsmc_rhit.RData")
results <- wsmc_rhit(compute_d, target, param_algo, savefile = filename)
wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps <- max(wsmc.df$step)
save(wsmc.df, results, nsteps, file = filename)

# load(filename)
# results <- wsmc_rhit(compute_d, target, param_algo)
# wsmccheck <- wsmc_fixedthresholds(results$threshold_history, compute_d, target, param_algo)
#
# wsmcrhit.df <- wsmc_to_dataframe(results, target$parameter_names)
# nsteps <- max(wsmcrhit.df$step)
# wsmccheck.df <- wsmc_to_dataframe(wsmccheck, target$parameter_names)
#
# summary(wsmccheck.df$step)
# summary(wsmcrhit.df$step)
#
# g <- ggplot(wsmcrhit.df, aes(x = mu, group = step)) + geom_density(aes(y = ..density..), colour = "darkgrey")
# g <- g +  theme(legend.position = "none")
# g <- g + xlab(expression(mu))
# g <- g + geom_density(data=wsmccheck.df  %>% filter(step <= nsteps), colour = "orange")
# g
#

# wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)

# gthreshold <- qplot(x = 1:(length(results$threshold_history)), y = results$threshold_history, geom = "line")
# gthreshold <- gthreshold + xlab("step") + ylab("threshold")
# gthreshold
# tail(results$threshold_history)
#
#
# head(wsmc.df)
#
# g <- ggplot(wsmc.df, aes(x = mu, y = sigma, colour = step, group = step))
# g <- g + geom_point(alpha = 0.5)
# g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
# g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
# g
#
# g <- ggplot(wsmc.df, aes(x = mu, colour = factor(step), group = step)) + geom_density(aes(y = ..density..))
# g <- g + geom_vline(xintercept = true_theta[1]) + theme(legend.position = "none")
# g
#
# ## load results from standard MCMC using the exact likelihood
# filename <- paste0("~/Dropbox/ABCD/Results/univariatenormal/normaldata.n", nobservations, ".metropolis.RData")
# load(filename)
# chainlist_to_dataframe <- function(chains_list){
#   nchains <- length(chains_list)
#   niterations <- nrow(chains_list[[1]])
#   chaindf <- foreach (i = 1:nchains, .combine = rbind) %do% {
#     data.frame(ichain = rep(i, niterations), iteration = 1:niterations, X = chains_list[[i]])
#   }
#   return(chaindf)
# }
# chaindf <- chainlist_to_dataframe(mh$chains)
# # plot
# chaindf.melt <- melt(chaindf, id.vars = c("ichain", "iteration"))
# mcmc.df <- chaindf.melt %>% filter(iteration > 2000) %>% spread(variable, value)
# #############
# g <- ggplot(wsmc.df, aes(x = mu, group = step)) + geom_density(aes(y = ..density..), colour = "darkgrey")
# g <- g +  theme(legend.position = "none")
# g <- g + geom_density(data = mcmc.df, aes(x = X.1, group = NULL, colour = NULL), colour = "black", linetype = 2)
# g <- g + xlab(expression(mu))
# g
#
# g <- ggplot(wsmc.df, aes(x = sigma, group = step)) + geom_density(aes(y = ..density..), colour = "darkgrey")
# g <- g + theme(legend.position = "none")
# g <- g + geom_density(data = mcmc.df, aes(x = X.2, group = NULL, colour = NULL), colour = "black", linetype = 2)
# g <- g + xlab(expression(sigma))
# g
#



