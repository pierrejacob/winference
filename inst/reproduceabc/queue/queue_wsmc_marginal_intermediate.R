library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()

set.seed(11)

prefix = ""

load(paste0(prefix, "50.intermediateobs.neal.Rdata"))
obs <- matrix(obs, nrow = 1)
obs_sorted = sort(obs)
nobs <- ncol(obs)
#
target <- get_queue()
target$simulate <- function(theta){
  return(matrix(target$robservation(nobs, theta, target$parameters), nrow = 1))
}

#x_obs <- cumsum(obs[1,])
#
param_algo <- list(nthetas = 2048, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 100000)
#compute_d <- get_hilbert_to_y(obs)
compute_d = function(y){
  sort_y = sort(y)
  mean(abs(sort_y-obs_sorted))
}


filename <- paste0(prefix, "queue_intermediate_wsmc_marginal.RData")
#results <- wsmc(compute_d, target, param_algo, savefile = filename, maxsim = 10^7)
load(filename)
results <- wsmc_continue(results, savefile = filename, maxsim = 2*10^7)
# load(filename)


target$rprior <- function(ntheta, parameters){
  theta1 <- runif(n = ntheta, min = 0, max = min(obs[1,]))
  theta2minus1 <- runif(n = ntheta, min = 0, max = 10)
  theta3 <- runif(n = ntheta, min = 0, max = 1/3)
  return(cbind(theta1, theta2minus1, theta3))
}
#
target$dprior <- function(thetas, parameters){
  evals <- dunif(thetas[,1], min = 0, max = min(obs[1,]), log = TRUE)
  evals <- evals + dunif(thetas[,2], min = 0, max = 10, log = TRUE)
  evals <- evals + dunif(thetas[,3], min = 0, max = 1/3, log = TRUE)
  return(evals)
}

filename <- paste0(prefix, "queue_intermediate_wsmc_marginal_constraints.RData")
#results <- wsmc(compute_d, target, param_algo, savefile = filename, maxsim = 10^7)
load(filename)
results <- wsmc_continue(results, savefile = filename, maxsim = 2*10^7)
# load(filename)


#
# plot_marginal(results, 1)
# plot_marginal(results, 2)
# plot_marginal(results, 3)
# #
# # get posterior mean and variance
# thetas <- tail(results$thetas_history, 1)[[1]]
# thetas[,3] <- log(thetas[,3])
# theta_mean <- colMeans(thetas)
# theta_cov <- cov(thetas)
#
# theta_mean
# theta_cov
