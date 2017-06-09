library(winference)
registerDoParallel(cores = 6)
rm(list = ls())
setmytheme()
set.seed(11)
# model
target <- get_autoregressive()
#
# number of observations
nobservations <- 1000
prefix <- ""
load(file = paste0(prefix, "ar1data.RData"))
obs <- obs[1:nobservations]
compute_d <- get_hilbert_to_y(matrix(obs, nrow = 1))

target$simulate <- function(theta){
  return(matrix(target$robservation(nobservations, theta, target$parameters, target$generate_randomness(nobservations)), nrow = 1))
}

thetas <- target$rprior(1, target$parameters)
y_sim <- target$simulate(thetas[1,])
compute_d(y_sim)

param_algo <- list(nthetas = 1024, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 1000)

filename <- paste0(prefix, "ar1.n", nobservations, ".wsmc_marginal.RData")
results <- wsmc(compute_d, target, param_algo, savefile = filename, maxtime = 60*60)
# load(filename)
# results <- wsmc_continue(results, savefile = filename, maxtime = 2*60*60)
#
# load(filename)
# wsmc.df <- wsmc_to_dataframe(results)
# nsteps <- max(wsmc.df$step)
# plot_bivariate_polygon(results, 1, 2)
# plot_bivariate(results, 1, 2)


