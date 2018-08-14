library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()
set.seed(11)

prefix = ""

target <- get_toggleswitch()

# number of observations
nobservations <- 2000
load(file = paste0(prefix,"toggleswitchdata.RData"))
obs <- obs[1:nobservations]
obs_sorted <- sort(obs)

# function to compute distance between observed data and data generated given theta
compute_d = function(y){
  sort_y = sort(y)
  mean(abs(sort_y-obs_sorted))
}

target$simulate <- function(theta) matrix(target$robservation(nobservations, theta, target$parameters, target$generate_randomness(nobservations)), nrow = 1)

#test
y_sim <- target$simulate(target$rprior(1, target$parameters))
compute_d(y_sim)

param_algo <- list(nthetas = 2048, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 1000)

filename <- paste0(prefix, "toggleswitchwsmc.n", nobservations, ".RData")
#results <- wsmc(compute_d, target, param_algo, savefile = filename, maxsimulation = 1e6)
load(filename)
#results <- wsmc_continue(results, savefile = filename, maxsim = 1e6)
# plot_bivariate(results, 1, 2)
