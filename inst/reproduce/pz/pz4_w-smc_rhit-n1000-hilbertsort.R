#+ presets, echo = FALSE, warning = FALSE, message = FALSE
library(winference)
library(doMC)
library(doRNG)
library(dplyr)
library(ggthemes)

registerDoMC(cores = 4)

rm(list = ls())
set.seed(11)
setmytheme()
#
# model
target <- get_pz_4param()
#
# number of observations
nobservations <- 1000
load(file = "~/Dropbox/ABCD/Results/data/pzdata.RData")
obs <- obs[1:nobservations]
# plot(obs, type = "l")
# now using the embeddings, and Hilbert sort
lagvalue <- 1
lag_obs <- create_lagmatrix(matrix(obs, nrow = 1), lagvalue)
lag_obs <- lag_obs[,(lagvalue+1):ncol(lag_obs)]
order_obs <- hilbert_order(lag_obs)
orderded_obs <- lag_obs[,order_obs]

compute_d_hilbert <- function(theta){
  fake_rand <- target$generate_randomness(nobservations)
  fake_obs <- target$robservation(nobservations, theta, target$parameters, fake_rand)
  fake_obs <- create_lagmatrix(matrix(fake_obs, nrow = 1), lagvalue)
  fake_obs <- fake_obs[,(lagvalue+1):ncol(fake_obs)]
  order_fake <- hilbert_order(fake_obs)
  distance <- nrow(fake_obs) * mean(abs(orderded_obs - fake_obs[,order_fake]))
  return(distance)
}

proposal <- mixture_proposal()

param_algo <- list(nthetas = 1024, nmoves = 1, proposal = proposal,
                   nsteps = 40, minimum_diversity = 0.5, R = 2, maxtrials = 1000)


filename <- paste0("~/Dropbox/ABCD/Results/pz/pzdata.n",
                   nobservations, ".L", lagvalue, ".wsmc_rhit-hilbert.RData")
# results <- wsmc_rhit(compute_d_hilbert, target, param_algo, savefile = filename)
# wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
# nsteps <- max(wsmc.df$step)
# save(wsmc.df, results, nsteps, file = filename)
load(filename)
wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps <- max(wsmc.df$step)

target$parameter_names

# g <- ggplot(wsmc.df, aes(x = omega, group = step)) + geom_density(aes(y = ..density..), colour = "darkgrey")
# g <- g +  theme(legend.position = "none")
# g <- g + xlab(expression(omega))
# g
#
# g <- ggplot(wsmc.df, aes(x = phi, group = step)) + geom_density(aes(y = ..density..), colour = "darkgrey")
# g <- g +  theme(legend.position = "none")
# g <- g + xlab(expression(phi))
# g

g <- ggplot(wsmc.df, aes(x = mu_alpha, y = sigma_alpha, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(mu[alpha])) + ylab(expression(sigma[alpha]))
g

g <- ggplot(wsmc.df, aes(x = c, y = e, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(c)) + ylab(expression(e))
g

dist.df <- foreach(irep = 1:nrow(results$thetas_history[[nsteps]]), .combine = rbind) %dorng%{
  c(compute_d_hilbert(results$thetas_history[[nsteps]][irep,]), compute_d_hilbert(true_theta))
}

dist.df <- melt(data.frame(dist.df))
ggplot(dist.df, aes(x = value, group = variable, fill = variable, colour = variable)) + geom_density(aes(y = ..density..), alpha = 0.5)

