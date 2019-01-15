library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()
set.seed(11)
# model
target <- get_autoregressive()
#
# number of observations
nobservations <- 1000
nparticles = 2048
prefix <- ""
load(file = paste0(prefix, "ar1data.RData"))
obs <- obs[1:nobservations]

#
lagvalue <- 1
lag_obs <- create_lagmatrix(matrix(obs, nrow = 1), lagvalue)
lag_obs <- lag_obs[,seq(from=1, to=ncol(lag_obs), by=2)]


compute_d <- get_transport_to_y(lag_obs)
compute_distance <- function(y_sim){
  lag_y_sim <- create_lagmatrix(matrix(y_sim, nrow = 1), lagvalue)
  lag_y_sim <- lag_y_sim[,seq(from=1, to=ncol(lag_y_sim), by=2)]
  compute_d(lag_y_sim)
}

target$simulate <- function(theta){
  return(matrix(target$robservation(nobservations, theta, target$parameters, target$generate_randomness(nobservations)), nrow = 1))
}

thetas <- target$rprior(1, target$parameters)
y_sim <- target$simulate(thetas[1,])
compute_distance(y_sim)

param_algo <- list(nthetas = nparticles, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 1000)

filename <- paste0(prefix, "ar1.n", nobservations, ".wsmc_delay1.RData")
results <- wsmc(compute_distance, target, param_algo, savefile = filename, maxsim = 10^6)
load(filename)
# results <- wsmc_continue(results, savefile = filename, maxtime = 10*60*60)
#
# load(filename)
# wsmc.df <- wsmc_to_dataframe(results)
# nsteps <- max(wsmc.df$step)
# plot_bivariate_polygon(results, 1, 2)
# plot_bivariate(results, 1, 2)
# plot_threshold_time(results)


