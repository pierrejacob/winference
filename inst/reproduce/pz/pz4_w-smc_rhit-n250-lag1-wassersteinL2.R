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
target <- get_pz_4param()
#
# number of observations
nobservations <- 250
load(file = "~/Dropbox/ABCD/Results/data/pzdata.RData")
obs <- obs[1:nobservations]
# plot(obs, type = "l")
# now using the embeddings, and Hilbert sort
lagvalue <- 1
lag_obs <- create_lagmatrix(matrix(obs, nrow = 1), lagvalue)
lag_obs <- lag_obs[,(lagvalue+1):ncol(lag_obs)]

eps <- 0.05
# compute regularized transport distance between delayed embeddings of time series
compute_d_wasserstein <- function(theta, transportiterations = 100){
  fake_rand <- target$generate_randomness(nobservations)
  fake_obs <- target$robservation(nobservations, theta, target$parameters, fake_rand)
  lag_fake_obs <- create_lagmatrix(matrix(fake_obs, nrow = target$ydim), lagvalue)
  C <- cost_matrix_L2(lag_obs, lag_fake_obs[,(lagvalue+1):nobservations,drop=FALSE])
  epsilon <- eps * median(C)
  equalw <- rep(1/(nobservations-lagvalue), (nobservations-lagvalue))
  wass <- wasserstein(equalw, equalw, C, epsilon, transportiterations)
  return(as.numeric(wass$distances))
}

proposal <- mixture_proposal()

param_algo <- list(nthetas = 1024, nmoves = 1, proposal = proposal,
                   nsteps = 40, minimum_diversity = 0.5, R = 2, maxtrials = 1000)


filename <- paste0("~/Dropbox/ABCD/Results/pz/pzdata.n",
                   nobservations, ".L", lagvalue, ".wsmc_rhit-wassersteinL2.RData")
# results <- wsmc_rhit(compute_d_wasserstein, target, param_algo, savefile = filename)
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
g <- g + geom_hline(yintercept = true_theta[2]) + geom_vline(xintercept = true_theta[1])
g

g <- ggplot(wsmc.df, aes(x = c, y = e, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(c)) + ylab(expression(e))
g <- g + geom_hline(yintercept = true_theta[4]) + geom_vline(xintercept = true_theta[3])
g

g + geom_rug(alpha = 0.2)

dist.df <- foreach(irep = 1:nrow(results$thetas_history[[nsteps]]), .combine = rbind) %dorng%{
  c(compute_d_wasserstein(results$thetas_history[[nsteps]][irep,]), compute_d_wasserstein(true_theta))
}

dist.df <- melt(data.frame(dist.df))
ggplot(dist.df, aes(x = value, group = variable, fill = variable, colour = variable)) + geom_density(aes(y = ..density..), alpha = 0.5)


