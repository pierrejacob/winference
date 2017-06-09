library(winference)
registerDoParallel(cores = 4)
rm(list = ls())
setmytheme()
set.seed(11)
target <- get_toggleswitch()
# number of observations
nobservations <- 2000
load(file = "toggleswitchdata.RData")
obs <- obs[1:nobservations]
obs_sorted <- sort(obs)

# function to compute distance between observed data and data generated given theta
compute_d <- get_hilbert_to_y(matrix(obs, nrow = 1))

target$simulate <- function(theta) matrix(target$robservation(nobservations, theta, target$parameters, target$generate_randomness(nobservations)), nrow = 1)

y_sim <- target$simulate(target$rprior(1, target$parameters))
compute_d(y_sim)

param_algo <- list(nthetas = 2048, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 1000)

filename <- paste0("toggleswitchwsmc.n", nobservations, ".RData")
load(filename)
results <- wsmc(compute_d, target, param_algo, savefile = filename, maxtime = 60*60)
# results <- wsmc_continue(results, savefile = filename, maxsim = 1e6)

