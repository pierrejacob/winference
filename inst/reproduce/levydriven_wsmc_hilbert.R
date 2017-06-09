library(winference)
registerDoParallel(cores = 6)
rm(list = ls())
setmytheme()

set.seed(13)

target <- get_levydriven()

# number of observations
prefix <- ""
load(file = paste0(prefix, "levydrivendata.RData"))

nobservations <- 10000
obs <- obs[1:nobservations]
plot(obs, type = "l")

lagvalue <- 1
# #
lag_obs <- create_lagmatrix(matrix(obs, nrow = 1), lagvalue)
lag_obs <- lag_obs[,seq(from=1, to=ncol(lag_obs), by=2)]

compute_hilbert <- get_hilbert_to_y(lag_obs)

compute_d <- function(z){
  fake_obs <- create_lagmatrix(matrix(z, nrow = 1), lagvalue)
  fake_obs <- fake_obs[,seq(from=1, to=ncol(fake_obs), by=2)]
  return(compute_hilbert(fake_obs))
}

target$simulate <- function(theta){
  r <- target$generate_randomness(nobservations)
  return(target$robservation(nobservations, theta, list(), r))
}

param_algo <- list(nthetas = 1024, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 1000)

filename <- paste0(prefix, "levydriven.n", nobservations, ".lag", lagvalue, ".wsmc.hilbert.RData")
results <- wsmc(compute_d, target, param_algo, savefile = filename, maxtime = 60*60)
load(file = filename)
# results <- wsmc_continue(results, savefile = filename, maxtime = 30*60)

grid.arrange(plot_threshold_time(results) + scale_y_log10(), plot_ncomputed(results))
