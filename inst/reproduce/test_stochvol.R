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

set.seed(11)

get_stochvol <- function(){
  rprior <- function(nparticles, ...){
    phis <- runif(nparticles, min = -1, max = 1)
    logsigmas <- rnorm(nparticles, 0, 1)
    return(cbind(phis, logsigmas))
  }
  dprior <- function(thetas, ...){
    if (is.null(dim(thetas))) thetas <- matrix(thetas, nrow = 1)
    logdensities <- dnorm(thetas[,2], 0, 1, log = TRUE)
    logdensities <- logdensities + dunif(thetas[,1], min = -1, max = 1, log = TRUE)
    return(logdensities)
  }
  #
  generate_randomness <- function(nobservations){
    return(list())
  }
  robservation <- function(nobservations, theta, parameters, randomness){
    obs <- rep(0, nobservations)
    phi <- theta[1]
    sigma <- exp(theta[2])
    state <- rnorm(1, mean = 0, sd = sigma / sqrt(1- (phi^2)))
    for (t in 1:nobservations){
      state <- phi * state + rnorm(1, mean = 0, sd = sigma)
      obs[t] <- exp(state) * rnorm(1, mean = 0, sd = 1)
    }
    return(obs)
  }
  #
  model <- list(rprior = rprior,
                dprior = dprior,
                generate_randomness = generate_randomness,
                robservation = robservation,
                parameter_names = c("phi", "log_sigma"),
                thetadim = 2, ydim = 1)
  return(model)
}

target <- get_stochvol()

# number of observations
nobservations <- 1000
# parameter of data-generating process
true_theta <- c(0.9, log(0.5))
obs <- target$robservation(nobservations, true_theta,
                           target$parameters, target$generate_randomness(nobservations))
plot(obs, type = "l")

lagvalue <- 1
#
lag_obs <- create_lagmatrix(matrix(obs, nrow = 1), lagvalue)
lag_obs <- lag_obs[,(lagvalue+1):ncol(lag_obs)]
order_obs <- hilbert_order(lag_obs)
orderded_obs <- lag_obs[,order_obs]

compute_d <- function(theta){
  fake_rand <- target$generate_randomness(nobservations)
  fake_obs <- target$robservation(nobservations, theta, target$parameters, fake_rand)
  fake_obs <- create_lagmatrix(matrix(fake_obs, nrow = 1), lagvalue)
  fake_obs <- fake_obs[,(lagvalue+1):ncol(fake_obs)]
  order_fake <- hilbert_order(fake_obs)
  distance <- nrow(fake_obs) * mean(abs(orderded_obs - fake_obs[,order_fake]))
  return(distance)
}


# obs_sorted <- sort(obs)
# # function to compute distance between observed data and data generated given theta
# compute_d <- function(theta, metric = metricL2){
#   fake_rand <- target$generate_randomness(nobservations)
#   fake_obs <- target$robservation(nobservations, theta, target$parameters, fake_rand)
#   fake_obs_sorted <- sort(fake_obs)
#   return(metric(obs_sorted, fake_obs_sorted))
# }
compute_d(true_theta)
proposal <- mixture_proposal()

param_algo <- list(nthetas = 1024, nmoves = 1, proposal = proposal,
                   nsteps = 20, minimum_diversity = 0.5, R = 2, maxtrials = 1000)

filename <- paste0("~/Dropbox/ABCD/Results/stochvol/stochvol.n",
                   nobservations, ".L", lagvalue, ".wsmc_rhit-hilbert.RData")
# results <- wsmc_rhit(compute_d, target, param_algo, savefile = filename)
# wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
# nsteps <- max(wsmc.df$step)
# save(wsmc.df, results, nsteps, file = filename)

load(filename)
wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps <- max(wsmc.df$step)
#
target$parameter_names

# g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = phi, group = step, colour = step)) + geom_density(aes(y = ..density..))
# g <- g +  theme(legend.position = "none")
# g <- g + xlab(expression(phi)) + geom_vline(xintercept = true_theta[1])
# g
#
g <- ggplot(wsmc.df, aes(x = phi, group = step)) + geom_density(aes(y = ..density..), colour = "darkgrey")
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(phi)) + geom_vline(xintercept = true_theta[1])
g

g <- ggplot(wsmc.df, aes(x = log_sigma, group = step)) + geom_density(aes(y = ..density..), colour = "darkgrey")
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(log(sigma))) + geom_vline(xintercept = true_theta[2])
g


g <- ggplot(wsmc.df, aes(x = phi, y = log_sigma, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(phi)) + ylab(expression(log(sigma)))
g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
g


