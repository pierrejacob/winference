library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()
set.seed(11)
prefix <- ""
target <- get_cosine()

nobservations <- 500

load(paste0(prefix, "cosinedata.RData"))
obs <- matrix(obs[1:nobservations], nrow = 1)
target$simulate <- function(theta){
  return(matrix(target$robservation(nobservations, theta, target$parameters, target$generate_randomness(nobservations)), nrow = 1))
}

param_algo <- list(nthetas = 1024, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 1e5)

compute_d <- function(z){
  return(sqrt(sum((z[1,] - obs[1,])^2)))
}


filename <- paste0(prefix, "cosine_wsmc_euclidean.n", nobservations, ".RData")
results <- wsmc(compute_d, target, param_algo, savefile = filename, maxsimulation = 5e6)
load(filename)
# results <- wsmc_continue(results, savefile = filename, maxsimulation = 1e6)

# plot_marginal(results, 1)
# plot_marginal(results, 2)
# plot_marginal(results, 3)
# plot_marginal(results, 4)
