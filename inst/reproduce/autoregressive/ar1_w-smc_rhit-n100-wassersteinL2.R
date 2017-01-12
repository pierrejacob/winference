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
target <- get_autoregressive()
#
# number of observations
nobservations <- 100
load(file = "~/Dropbox/ABCD/Results/data/ar1data.RData")
obs <- obs[1:nobservations]

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

param_algo <- list(nthetas = 1024, nmoves = 5, proposal = proposal,
                   nsteps = 20, minimum_diversity = 0.5, R = 2, maxtrials = 1000)


filename <- paste0("~/Dropbox/ABCD/Results/autoregressive/ar1data.n",
                   nobservations, ".wsmc_rhit-wassersteinL2.RData")
# results <- wsmc_rhit(compute_d_wasserstein, target, param_algo, savefile = filename)
# wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
# nsteps <- max(wsmc.df$step)
# save(wsmc.df, results, nsteps, file = filename)
load(filename)
wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps <- max(wsmc.df$step)

wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps <- max(wsmc.df$step)

g <- ggplot(wsmc.df, aes(x = rho, group = step)) + geom_density(aes(y = ..density..), colour = "darkgrey")
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(rho))
g

g <- ggplot(wsmc.df, aes(x = logsigma, group = step)) + geom_density(aes(y = ..density..), colour = "darkgrey")
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(log(sigma)))
g

g <- ggplot(wsmc.df, aes(x = rho, y = logsigma, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(rho)) + ylab(expression(log(sigma)))
g


filename <- paste0("~/Dropbox/ABCD/Results/autoregressive/ar1data.n", nobservations, ".metropolis.RData")
load(filename)
chainlist_to_dataframe <- function(chains_list){
  nchains <- length(chains_list)
  niterations <- nrow(chains_list[[1]])
  chaindf <- foreach (i = 1:nchains, .combine = rbind) %do% {
    data.frame(ichain = rep(i, niterations), iteration = 1:niterations, X = chains_list[[i]])
  }
  return(chaindf)
}
chaindf <- chainlist_to_dataframe(mh$chains)
# plot
chaindf.melt <- melt(chaindf, id.vars = c("ichain", "iteration"))
chain.bycomponent.df <- chaindf.melt %>% filter(iteration > 2000) %>% spread(variable, value)

# g + geom_density2d(data=chain.bycomponent.df, aes(x = X.1, y = X.2, colour = NULL, group = NULL),
                   # colour = "black")
g + geom_point(data=chain.bycomponent.df %>% filter(iteration %% 10 == 1), aes(x = X.1, y = X.2, colour = NULL, group = NULL),
               colour = "black")

g <- ggplot(wsmc.df %>% filter(step > 10), aes(x = rho, y = logsigma, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + theme(legend.position = "none")
g <- g + xlab(expression(rho)) + ylab(expression(log(sigma)))
g + geom_point(data=chain.bycomponent.df %>% filter(iteration %% 10 == 1), aes(x = X.1, y = X.2, colour = NULL, group = NULL),
               colour = "orange", alpha = 0.1)

g <- ggplot(wsmc.df %>% filter(step == 15), aes(x = rho, y = logsigma, colour = step, group = step))
g <- g + geom_density2d()
g <- g + theme(legend.position = "none")
g <- g + xlab(expression(rho)) + ylab(expression(log(sigma)))
g + geom_density2d(data=chain.bycomponent.df, aes(x = X.1, y = X.2, colour = NULL, group = NULL),
               colour = "black")

g <- ggplot(wsmc.df %>% filter(step == 15), aes(x = rho))
g <- g + geom_histogram(aes(y = ..density..), alpha = 0.5)
g <- g + theme(legend.position = "none")
g <- g + xlab(expression(rho))
g <- g + geom_histogram(data=chain.bycomponent.df, aes(x = X.1, y = ..density..), fill = "red", alpha = 0.2)
g

g <- ggplot(wsmc.df %>% filter(step == 15), aes(x = logsigma))
g <- g + geom_histogram(aes(y = ..density..), alpha = 0.5)
g <- g + theme(legend.position = "none")
g <- g + xlab(expression(log(sigma)))
g <- g + geom_histogram(data=chain.bycomponent.df, aes(x = X.2, y = ..density..), fill = "red", alpha = 0.2)
g
