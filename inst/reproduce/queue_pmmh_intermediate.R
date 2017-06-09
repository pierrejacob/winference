library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()

set.seed(11)
prefix <- ""

load(paste0(prefix, "50.intermediateobs.neal.Rdata"))
obs <- matrix(obs, nrow = 1)
nobs <- ncol(obs)
#
target <- get_queue()
x_obs <- cumsum(obs[1,])
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

#
target$rinit <- function(nparticles, theta){
  return(matrix(rexp(n = nparticles, rate = theta[3]), ncol = 1))
}

target$rtransition <- function(xparticles, time, theta){
  return(xparticles + matrix(rexp(n = nrow(xparticles), rate = theta[3]), ncol = 1))
}

target$dobs <- function(observations, time, xparticles, theta){
  if (time > 1){
    return(dunif(observations[1,time], min = theta[1] + pmax(0, xparticles[,1] - x_obs[time-1]),
                 max = theta[1] + theta[2] + pmax(0, xparticles[,1] - x_obs[time-1]), log = TRUE))
  } else {
    return(dunif(observations[1,time], min = theta[1] + xparticles[,1], max = theta[1] + theta[2] + xparticles[,1], log = TRUE))
  }
}

## load result from WSMC
filename <- paste0(prefix, "queue_intermediate_wsmc_marginal.RData")
load(filename)
# get posterior mean and variance
thetas <- tail(results$thetas_history, 1)[[1]]
theta_mean <- colMeans(thetas)
theta_cov <- cov(thetas)
#
theta_init <- thetas[sample(x = 1:nrow(thetas), 2, replace = T),]
while (any(is.infinite(target$dprior(theta_init, target$parameters)))){
  theta_init <- thetas[sample(x = 1:nrow(thetas), 2, replace = T),]
}

tuning_parameters <- list(niterations = 100000, nparticles = 2^12, nchains = nrow(theta_init),
                          cov_proposal = theta_cov / 3,
                          adaptation = 10000, init_chains = theta_init)
savefile <- paste0(prefix, "queue_intermediate_pmmh.RData")
res <- pmmh(obs, target, tuning_parameters, savefile = savefile)
save(res, file = savefile)
# load(file = savefile)
# res <- pmmh_results
# mcmc.df <- mhchainlist_to_dataframe(res$chains)
# mcmc.df <- mcmc.df %>% filter(iteration < res$iteration)
# mcmc.df %>% tail
# ggplot(mcmc.df %>% filter(iteration %% 10 == 1), aes(x = iteration, y = X.1, group = ichain, colour = factor(ichain))) + geom_line(alpha = 0.5)
# ggplot(mcmc.df %>% filter(iteration %% 10 == 1), aes(x = iteration, y = X.2, group = ichain, colour = factor(ichain))) + geom_line(alpha = 0.5)
# ggplot(mcmc.df %>% filter(iteration %% 10 == 1), aes(x = iteration, y = X.3, group = ichain, colour = factor(ichain))) + geom_line(alpha = 0.5)
# ggplot(mcmc.df %>% filter(iteration %% 10 == 1), aes(x = X.1, y = X.2, group = ichain, colour = factor(ichain))) + geom_point(alpha = 0.5)
# ggplot(mcmc.df %>% filter(iteration %% 10 == 1), aes(x = X.1, y = X.3, group = ichain, colour = factor(ichain))) + geom_point(alpha = 0.5)
# ggplot(mcmc.df %>% filter(iteration %% 10 == 1), aes(x = X.2, y = X.3, group = ichain, colour = factor(ichain))) + geom_point(alpha = 0.5)
#
# ggplot(mcmc.df, aes(x = X.1, group = ichain)) + geom_density(aes(y = ..density..))
# ggplot(mcmc.df, aes(x = X.2, group = ichain)) + geom_density(aes(y = ..density..))
# ggplot(mcmc.df, aes(x = X.3, group = ichain)) + geom_density(aes(y = ..density..))
#
#
# mcmc.df %>% group_by(ichain) %>% filter(iteration > 3000) %>% summarise(m1 = mean(X.1), m2 = mean(X.2), m3 = mean(log(X.3)),
#                                                                         sd1 = sd(X.1), sd2 = sd(X.2), sd3 = sd(log(X.3)))
