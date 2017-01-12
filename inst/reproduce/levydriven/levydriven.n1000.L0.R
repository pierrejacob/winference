library(winference)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(foreach)
library(doMC)
library(doRNG)
library(reshape2)
registerDoMC(cores = 6)
rm(list = ls())
setmytheme()

set.seed(13)

target <- get_levydriven()

# number of observations
load(file = "~/Dropbox/ABCD/Results/data/levydrivendata.RData")

nobservations <- 1000
obs <- obs[1:nobservations]
plot(obs, type = "l")

lagvalue <- 0
#
# lag_obs <- create_lagmatrix(matrix(obs, nrow = 1), lagvalue)
# lag_obs <- lag_obs[,(lagvalue+1):ncol(lag_obs)]
# order_obs <- hilbert_order(lag_obs)
# orderded_obs <- lag_obs[,order_obs]

# compute_d <- function(theta){
#   fake_rand <- target$generate_randomness(nobservations)
#   fake_obs <- target$robservation(nobservations, theta, target$parameters, fake_rand)
#   fake_obs <- create_lagmatrix(matrix(fake_obs, nrow = 1), lagvalue)
#   fake_obs <- fake_obs[,(lagvalue+1):ncol(fake_obs)]
#   order_fake <- hilbert_order(fake_obs)
#   distance <- nrow(fake_obs) * mean(abs(orderded_obs - fake_obs[,order_fake]))
#   return(distance)
# }


# ##### function to compute distance between observed data and data generated given theta
obs_sorted <- sort(obs)
compute_d <- function(theta, metric = metricL2){
  fake_rand <- target$generate_randomness(nobservations)
  fake_obs <- target$robservation(nobservations, theta, target$parameters, fake_rand)
  fake_obs_sorted <- sort(fake_obs)
  return(metric(obs_sorted, fake_obs_sorted))
}

compute_d(true_theta)

proposal <- mixture_proposal()

param_algo <- list(nthetas = 1024, nmoves = 1, proposal = proposal,
                   nsteps = 20, minimum_diversity = 0.5, R = 2, maxtrials = 1000)

filename <- paste0("~/Dropbox/ABCD/Results/levydriven/levydriven.n",
                   nobservations, ".L", lagvalue, ".wsmc_rhit-hilbert.RData")
# results <- wsmc_rhit(compute_d, target, param_algo, savefile = filename)
# wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
# nsteps <- max(wsmc.df$step)
# save(wsmc.df, results, nsteps, param_algo, file = filename)

load(filename)
wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps <- max(wsmc.df$step)
#
results$threshold_history
target$parameter_names

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = mu, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(mu)) + geom_vline(xintercept = true_theta[1])
g

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = beta, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(beta)) + geom_vline(xintercept = true_theta[2])
g

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = xi, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(xi)) + geom_vline(xintercept = true_theta[3])
g

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = omega2, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(omega^2)) + geom_vline(xintercept = true_theta[4])
g

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = lambda, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(lambda)) + geom_vline(xintercept = true_theta[5])
g

##

g <- ggplot(wsmc.df, aes(x = mu, y = beta, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(mu)) + ylab(expression(beta))
g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
g

g <- ggplot(wsmc.df, aes(x = xi, y = omega2, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(xi)) + ylab(expression(omega^2))
g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
g

g <- ggplot(wsmc.df, aes(x = omega2, y = lambda, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + ylab(expression(lambda)) + xlab(expression(omega^2))
g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
g



dist.df <- foreach(irep = 1:nrow(results$thetas_history[[nsteps]]), .combine = rbind) %dorng%{
  c(compute_d(results$thetas_history[[nsteps]][irep,]), compute_d(true_theta))
}

dist.df <- melt(data.frame(dist.df))
dist.df %>% head

ggplot(dist.df, aes(x = value, group = variable, fill = variable, colour = variable)) + geom_density(aes(y = ..density..), alpha = 0.5)
