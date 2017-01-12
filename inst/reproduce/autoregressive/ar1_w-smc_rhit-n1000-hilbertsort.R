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
nobservations <- 1000
load(file = "~/Dropbox/ABCD/Results/data/ar1data.RData")
obs <- obs[1:nobservations]

# now using the embeddings, and Hilbert sort
lagvalue <- 1
lag_obs <- create_lagmatrix(matrix(obs, nrow = 1), lagvalue)
lag_obs <- lag_obs[,(lagvalue+1):ncol(lag_obs)]
order_obs <- hilbert_order(lag_obs)
orderded_obs <- lag_obs[,order_obs]
# compute regularized transport distance between delayed embeddings of time series
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

param_algo <- list(nthetas = 1024, nmoves = 5, proposal = proposal,
                   nsteps = 20, minimum_diversity = 0.5, R = 2, maxtrials = 1000)


filename <- paste0("~/Dropbox/ABCD/Results/autoregressive/ar1data.n",
                   nobservations, ".wsmc_rhit-hilbert.RData")
# results <- wsmc_rhit(compute_d_hilbert, target, param_algo, savefile = filename)
# wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
# nsteps <- max(wsmc.df$step)
# save(wsmc.df, results, nsteps, file = filename)
load(filename)
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


##### load posterior samples obtained with MH
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

g <- ggplot(wsmc.df %>% filter(step == 17), aes(x = rho, y = logsigma, colour = factor(step), group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + theme(legend.position = "none")
g <- g + xlab(expression(rho)) + ylab(expression(log(sigma)))
g + geom_point(data=chain.bycomponent.df %>% filter(iteration %% 10 == 1), aes(x = X.1, y = X.2, colour = NULL, group = NULL),
               colour = "black")

g <- ggplot(wsmc.df %>% filter(step == 17), aes(x = rho, y = logsigma, colour = factor(step), group = step))
g <- g + geom_density2d()
g <- g + theme(legend.position = "none")
g <- g + xlab(expression(rho)) + ylab(expression(log(sigma)))
g + geom_density2d(data=chain.bycomponent.df, aes(x = X.1, y = X.2, colour = NULL, group = NULL),
                   colour = "black")



