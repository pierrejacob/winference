#+ presets, echo = FALSE, warning = FALSE, message = FALSE
library(winference)
library(doMC)
library(doRNG)
library(dplyr)
library(ggthemes)

registerDoMC(cores = 8)

rm(list = ls())
set.seed(13)
setmytheme()
#
# model
target <- get_pz_4param()
target$robservation <- function(nobservations, theta, parameters, randomness){
  # untr_theta <- untransform_theta(theta)
  state <- matrix(exp(log(2) + randomness$x_0), nrow = 2)
  states <- matrix(nrow = 2, ncol = nobservations+1)
  states[,1] <- state
  log_obs <- rep(0, nobservations)
  for (t in 1:nobservations){
    alpha <- theta[1] #theta[2] * randomness$x[t] + theta[1]
    state <- pz_transition(state, alpha, t-1, c(theta[3:4], 0.1, 0.1))
    states[,t+1] <- state
    log_obs[t] <- 0.2*randomness$y[t] + log(state[1,1])
  }
  return(log_obs)
}

# plot(obs, type = "l")
reconstitute <- function(observations, theta){
  nobservations <- length(observations)
  state <- matrix(exp(log(2)), nrow = 2) # matrix(exp(log(2) + rnorm(2, 0, 1)), nrow = 2)
  states <- matrix(nrow = 2, ncol = nobservations+1)
  states[,1] <- state
  residuals <- rep(0, nobservations)
  for (t in 1:nobservations){
    alpha <- theta[1] # theta[2] * rnorm(1) + theta[1]
    state <- pz_transition(state, alpha, t-1, c(theta[3:4], 0.1, 0.1))
    states[,t+1] <- state
    residuals[t] <- (observations[t] - log(state[1,1]))/0.2
  }
  return(residuals)
}


#
# number of observations
nobservations <- 500
load(file = "~/Dropbox/ABCD/Results/data/pzdata.RData")
obs <- obs[1:nobservations]
# plot(obs, type = "l")

# now using the embeddings, and Hilbert sort
lagvalue <- 10
tau <- 2
create_lagmatrix <- function(timeseries, k, tau){
  res <- matrix(NA, nrow = k+1, ncol = ncol(timeseries))
  res[1,] <- timeseries
  for (lagvalue in 1:k){
    res[lagvalue+1,] <- lag(timeseries[1,], lagvalue*tau)
  }
  return(res)
}
# lag_obs <- create_lagmatrix(matrix(obs, nrow = 1), lagvalue, tau)
# lag_obs <- lag_obs[,(lagvalue*tau+1):ncol(lag_obs)]
# order_obs <- hilbert_order(lag_obs)
# orderded_obs <- lag_obs[,order_obs]

compute_d <- function(theta){
  eps <- reconstitute(obs, theta)
  lag_obs <- create_lagmatrix(matrix(eps, nrow = 1), lagvalue, tau)
  lag_obs <- lag_obs[,(lagvalue*tau+1):ncol(lag_obs)]
  order_obs <- hilbert_order(lag_obs)
  orderded_obs <- lag_obs[,order_obs]
  fake_eps <- rnorm(length(obs))
  lag_fake_obs <- create_lagmatrix(matrix(fake_eps, nrow = 1), lagvalue, tau)
  lag_fake_obs <- lag_fake_obs[,(lagvalue*tau+1):ncol(lag_fake_obs)]
  order_fake_obs <- hilbert_order(lag_fake_obs)
  orderded_fake_obs <- lag_fake_obs[,order_fake_obs]
  distance <- nrow(orderded_obs) * mean(abs(orderded_obs - orderded_fake_obs))
  return(distance)
}

compute_d(true_theta)

proposal <- mixture_proposal()

param_algo <- list(nthetas = 1024, nmoves = 1, proposal = proposal,
                   nsteps = 40, minimum_diversity = 0.5, R = 2, maxtrials = 1000)


filename <- paste0("~/Dropbox/ABCD/Results/pz/pzdata.n",
                   nobservations, ".L", lagvalue, ".tau", tau, ".wsmc_rhit-residuals-hilbert.RData")
# results <- wsmc_rhit(compute_d, target, param_algo, savefile = filename)
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

dist.df <- foreach(irep = 1:40, .combine = rbind) %dorng%{
  c(compute_d(results$thetas_history[[nsteps]][irep,]), compute_d(true_theta))
}

dist.df <- melt(data.frame(dist.df))
ggplot(dist.df, aes(x = value, group = variable, fill = variable, colour = variable)) + geom_density(aes(y = ..density..), alpha = 0.5)

