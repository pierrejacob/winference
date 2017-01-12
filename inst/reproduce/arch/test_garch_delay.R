library(winference)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(foreach)
library(doMC)
library(doRNG)
library(reshape2)
registerDoMC(cores = 4)
rm(list = ls())
setmytheme()

set.seed(11)

rprior <- function(nparticles, ...){
  dummy <- rexp(nparticles, rate = 1)
  alpha1_unnormalized <- rexp(nparticles, rate = 1)
  beta1_unnormalized <- rexp(nparticles, rate = 1)
  alpha0 <- rexp(nparticles, rate = 1)
  alpha1 <- alpha1_unnormalized / (dummy + alpha1_unnormalized + beta1_unnormalized)
  beta1 <- beta1_unnormalized / (dummy + alpha1_unnormalized + beta1_unnormalized)
  return(cbind(alpha0, alpha1, beta1))
}

dprior <- function(thetas, ...){
  if (is.null(dim(thetas))) thetas <- matrix(thetas, nrow = 1)
  logdensities <- dexp(thetas[,1], rate = 1, log = TRUE)
  logdensities <- logdensities + log(thetas[,2] + thetas[,3] < 1 & thetas[,2] > 0 & thetas[,3] > 0)
  return(logdensities)
}

generate_randomness <- function(nobservations){
  return(list())
}

robservation <- function(nobservations, theta, parameters, randomness){
  obs <- rnorm(nobservations, mean = 0, sd = 1)
  variance <- rep(1, nobservations)
  alpha0 <- theta[1]
  alpha1 <- theta[2]
  beta1 <- theta[3]
  for (t in 2:nobservations){
    variance[t] <- alpha0 + alpha1 * obs[t-1]^2 + beta1 * variance[t-1]
    obs[t] <- sqrt(variance[t]) * obs[t]
  }
  return(obs)
}

garch_model <- list(rprior = rprior,
                   dprior = dprior,
                   generate_randomness = generate_randomness,
                   robservation = robservation,
                   parameter_names = c("alpha_0", "alpha_1", "beta_1"),
                   thetadim = 3, ydim = 1)

true_theta <- c(0.1, 0.15, 0.6)
nobservations <- 250
obs <- robservation(nobservations, true_theta, list(), list())

lagvalue <- 3
lag_obs <- create_lagmatrix(matrix(obs, nrow = 1), lagvalue)
lag_obs <- lag_obs[,(lagvalue+1):ncol(lag_obs)]

eps <- 0.05
# compute regularized transport distance between delayed embeddings of time series
compute_d <- function(theta, transportiterations = 100){
  fake_rand <- garch_model$generate_randomness(nobservations)
  fake_obs <- garch_model$robservation(nobservations, theta, garch_model$parameters, fake_rand)
  lag_fake_obs <- create_lagmatrix(matrix(fake_obs, nrow = garch_model$ydim), lagvalue)
  C <- cost_matrix_L2(lag_obs, lag_fake_obs[,(lagvalue+1):nobservations,drop=FALSE])
  epsilon <- eps * median(C)
  equalw <- rep(1/(nobservations-lagvalue), (nobservations-lagvalue))
  wass <- wasserstein(equalw, equalw, C, epsilon, transportiterations)
  d <- as.numeric(wass$distances)
  if (is.na(d)){
    d <- Inf
  }
  return(d)
}

proposal <- mixture_proposal()

param_algo <- list(nthetas = 1024, nmoves = 1, proposal = proposal,
                   nsteps = 50, minimum_diversity = 0.5, R = 2, maxtrials = 1000)

filename <- paste0("~/Dropbox/ABCD/Results/arch/garch.n",
                   nobservations, "delaywasserstein.wsmc_rhit.RData")
# results <- wsmc_rhit(compute_d, garch_model, param_algo, savefile = filename)
# wsmc.df <- wsmc_to_dataframe(results, garch_model$parameter_names)
# nsteps <- max(wsmc.df$step)
# save(wsmc.df, results, nsteps, file = filename)

load(filename)
wsmc.df <- wsmc_to_dataframe(results, garch_model$parameter_names)
nsteps <- max(wsmc.df$step)
#
g <- ggplot(wsmc.df, aes(x = alpha_0, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g1 <- g + xlab(expression(alpha[0])) + geom_vline(xintercept = true_theta[1])
g1

g <- ggplot(wsmc.df, aes(x = alpha_1, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g2 <- g + xlab(expression(alpha[1])) + geom_vline(xintercept = true_theta[2])
g2

g <- ggplot(wsmc.df, aes(x = beta_1, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g3 <- g + xlab(expression(beta[1])) + geom_vline(xintercept = true_theta[3])
g3

g <- ggplot(wsmc.df, aes(x = alpha_0, y = alpha_1, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(alpha[0])) + ylab(expression(alpha[1]))
g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
g

g <- ggplot(wsmc.df, aes(x = alpha_1, y = beta_1, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(alpha[1])) + ylab(expression(beta[1]))
g <- g + geom_vline(xintercept = true_theta[2]) + geom_hline(yintercept = true_theta[3])
g

# g + geom_rug(alpha = 0.2)
g + geom_rug(data = wsmc.df %>% filter(step == nsteps), alpha = 0.2)


dist.df <- foreach(irep = 1:nrow(results$thetas_history[[nsteps]]), .combine = rbind) %dorng%{
  c(compute_d(results$thetas_history[[nsteps]][irep,]), compute_d(true_theta))
}

dist.df <- melt(data.frame(dist.df))
dist.df %>% head

ggplot(dist.df, aes(x = value, group = variable, fill = variable, colour = variable)) + geom_density(aes(y = ..density..), alpha = 0.5)
