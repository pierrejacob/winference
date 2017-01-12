library(winference)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(foreach)
library(doMC)
library(doRNG)
library(reshape2)
registerDoMC(cores = 2)
rm(list = ls())
setmytheme()

set.seed(11)

get_blowfly <- function(){
# prior https://github.com/wittawatj/k2abc/blob/master/code/blowfly/sample_from_prior_blowflydata.m
  rprior <- function(nparticles,...){
    log_P = 2 + 2.*rnorm(nparticles,0,1);
    log_delta = -1 + 0.4*rnorm(nparticles,0,1);
    log_N0 = 5 + 0.5*rnorm(nparticles,0,1);
    log_sigd = -0.5 + rnorm(nparticles,0,1);
    log_sigp = -0.5 + rnorm(nparticles, 0,1);
    log_tau = 2 + rnorm(nparticles,0,1);
    return(cbind(log_P, log_delta, log_N0, log_sigd, log_sigp, log_tau))
  }
  dprior <- function(thetas,...){
    if (is.null(dim(thetas))) thetas <- matrix(thetas, nrow = 1)
    densityevals <- dnorm(thetas[,1], mean = 2, sd = 2, log = TRUE)
    densityevals <- densityevals + dnorm(thetas[,2], mean = -1, sd = 0.4, log = TRUE)
    densityevals <- densityevals + dnorm(thetas[,3], mean = 5, sd = 0.5, log = TRUE)
    densityevals <- densityevals + dnorm(thetas[,4], mean = -0.5, sd = 1, log = TRUE)
    densityevals <- densityevals + dnorm(thetas[,5], mean = -0.5, sd = 1, log = TRUE)
    densityevals <- densityevals + dnorm(thetas[,6], mean = 2, sd = 1, log = TRUE)
    return(densityevals)
  }
  #
  generate_randomness <- function(nobservations){
    return(list())
  }
  # generate data https://github.com/wittawatj/k2abc/blob/master/code/blowfly/gendata_pop_dyn_eqn.m
  robservation <- function(nobservations, theta, parameters, randomness){
    P <- exp(theta[1])
    delta <- exp(theta[2])
    N0 <- exp(theta[3])
    sigma_d <- exp(theta[4])
    sigma_p <- exp(theta[5])
    tau <- floor(exp(theta[6]))
    if (tau == 0) tau <- 1
    burnin <- 50
    obs <- rep(0, nobservations+burnin+tau)
    obs[1:tau] <- 180
    eps_s <- rgamma(nobservations+burnin, shape = 1/(sigma_d^2), scale = sigma_d^2)
    e_s <- rgamma(nobservations+burnin, shape = 1/(sigma_p^2), scale = sigma_p^2)
    for (i in 1:(nobservations+burnin)){
      t <- i + tau
      eps <- eps_s[i]
      e <- e_s[i]
      tau_t <- t - tau
      obs[t] <- P * obs[tau_t] * exp(-obs[tau_t] / N0) * e + obs[t-1] * exp(-delta * eps);
    }
    return(obs[(burnin+tau+1):(nobservations+burnin+tau)])
  }
  #
  model <- list(rprior = rprior,
                dprior = dprior,
                generate_randomness = generate_randomness,
                robservation = robservation,
                parameter_names = c("logP", "logdelta", "logN0", "logsigma_d", "logsigma_p", "logtau"),
                thetadim = 6, ydim = 1)
  return(model)
}

# data at https://github.com/wittawatj/k2abc/blob/master/code/blowfly/blowflydata.m
flydata <- read.csv("~/Dropbox/ABCD/Results/blowfly/flydata.csv", header = FALSE, sep = ",")
obs <- flydata[,2]
plot(obs, type = "l")
target <- get_blowfly()
nobservations <- length(obs)
test_theta <- c(3.76529501, -1.03266828, 5.46587492, -0.40094812, -0.96334847,  log(7))
fake <- target$robservation(nobservations, test_theta, list(), list())
matplot(cbind(fake, obs), type = "l")

lagvalue <- 5

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

proposal <- mixture_proposal()

param_algo <- list(nthetas = 2048, nmoves = 1, proposal = proposal,
                   nsteps = 40, minimum_diversity = 0.5, R = 2, maxtrials = 1000)

filename <- paste0("~/Dropbox/ABCD/Results/blowfly/blowfly.n",
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

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = logP, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(log(P)))
g
#
# g <- ggplot(wsmc.df, aes(x = phi, group = step)) + geom_density(aes(y = ..density..), colour = "darkgrey")
# g <- g +  theme(legend.position = "none")
# g <- g + xlab(expression(phi))
# g

g <- ggplot(wsmc.df, aes(x = logP, y = logdelta, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
# g <- g + xlab(expression(log(r))) + ylab(expression(phi))
g <- g + geom_vline(xintercept = test_theta[1]) + geom_hline(yintercept = test_theta[2])
g

g <- ggplot(wsmc.df, aes(x = logN0, y = logtau, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
# g <- g + xlab(expression(log(r))) + ylab(expression(sigma[e]))
g <- g + geom_vline(xintercept = test_theta[3]) + geom_hline(yintercept = test_theta[6])
g

g <- ggplot(wsmc.df, aes(x = logsigma_d, y = logsigma_p, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + geom_vline(xintercept = test_theta[4]) + geom_hline(yintercept = test_theta[5])
g


sim_series <- foreach(i = 1:2048, .combine = rbind) %dorng% {
  z <- target$robservation(nobservations, results$thetas_history[[nsteps]][i,], list(), list())
  z
}
mean_series <- apply(sim_series, 2, mean)


plot(obs, type = "l")
lines(mean_series)
