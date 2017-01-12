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
target <- get_cosine()
#
# number of observations
nobservations <- 500
load(file = "~/Dropbox/ABCD/Results/data/cosinedata.RData")
obs <- obs[1:nobservations]
# plot(obs, type = "l")
# now using the embeddings, and Hilbert sort
lagvalue <- 3
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

param_algo <- list(nthetas = 1024, nmoves = 5, proposal = proposal,
                   nsteps = 20, minimum_diversity = 0.5, R = 2, maxtrials = 1000)


filename <- paste0("~/Dropbox/ABCD/Results/cosine/cosinedata.n",
                   nobservations, ".L", lagvalue, ".wsmc_rhit-wassersteinL2.RData")
load(filename)
wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps <- max(wsmc.df$step)

g <- ggplot(wsmc.df, aes(x = omega, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(omega))
g <- g + geom_vline(xintercept = true_theta[1])
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
g
ggsave(filename = "~/Dropbox/ABCD/draft5/cosine_withlag3_omega.pdf", plot = g, width = 5, height = 5 )

g <- ggplot(wsmc.df, aes(x = phi, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(phi))
g <- g + geom_vline(xintercept = true_theta[2])
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
g
ggsave(filename = "~/Dropbox/ABCD/draft5/cosine_withlag3_phi.pdf", plot = g, width = 5, height = 5 )

g <- ggplot(wsmc.df, aes(x = logsigma, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(log(sigma)))
g <- g + geom_vline(xintercept = true_theta[3])
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
g
ggsave(filename = "~/Dropbox/ABCD/draft5/cosine_withlag3_logsigma.pdf", plot = g, width = 5, height = 5 )

g <- ggplot(wsmc.df, aes(x = logA, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(log(A)))
g <- g + geom_vline(xintercept = true_theta[4])
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
g
ggsave(filename = "~/Dropbox/ABCD/draft5/cosine_withlag3_logA.pdf", plot = g, width = 5, height = 5 )

#

g <- ggplot(wsmc.df, aes(x = omega, y = logsigma, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(omega)) + ylab(expression(log(sigma)))
g

g <- ggplot(wsmc.df, aes(x = phi, y = logA, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(phi)) + ylab(expression(log(A)))
g

