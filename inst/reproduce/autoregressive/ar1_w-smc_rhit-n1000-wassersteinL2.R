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
                   nsteps = 20, minimum_diversity = 0.5, R = 2, maxtrials = 1000)


filename <- paste0("~/Dropbox/ABCD/Results/autoregressive/ar1data.n",
                   nobservations, ".wsmc_rhit.Lag1.wassersteinL2.RData")
# results <- wsmc_rhit(compute_d_wasserstein, target, param_algo, savefile = filename)
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


#
# compute_d_wasserstein(true_theta)
#
# param_algo$compute_d <- compute_d_wasserstein
# filename <- paste0("~/Dropbox/ABCD/Results/ar1.smc1lagwasserstein.n", nobservations, ".RData")
# # results_wasserstein <- smcsampler(target, param_algo, filename)
# load(filename)
#
# results_wasserstein <- results
# th <- results_wasserstein$thetas_history
# nsteps <- length(th) - 1
# df <- data.frame()
# for (i in 1:(nsteps+1)){
#   df_ <- data.frame(cbind(th[[i]][,1:target$thetadim], rep(i, nrow(th[[i]]))))
#   names(df_) <- c("phi", "sigma", "step")
#   df <- rbind(df, df_)
# }
# names(df) <- c("phi", "sigma", "step")
#
# g <- ggplot(df, aes(x = phi, y = sigma, colour = step, group = step))
# g <- g + geom_point(alpha = 0.5)
# g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2), low = "orange", mid = "black", high = "white") + theme(legend.position = "none")
# g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
# g
#
# gthreshold <- qplot(x = 1:(length(results_wasserstein$threshold_history)), y = results_wasserstein$threshold_history, geom = "line")
# gthreshold <- gthreshold + xlab("step") + ylab("threshold")
# gthreshold
#
#
# thetas <- th[[nsteps+1]]
# wassersteindistances  <- results_wasserstein$distances_history[[nsteps+1]]
# summary(wassersteindistances)
# hilbertdistances <- as.numeric(foreach(i = 1:nrow(thetas), .combine = c) %dorng% {
#   theta <- thetas[i,]
#   compute_d_hilbert(theta)
# })
#
#
# plot(wassersteindistances, hilbertdistances)
