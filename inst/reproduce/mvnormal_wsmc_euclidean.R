library(winference)
registerDoParallel(cores = 8)
rm(list = ls())
setmytheme()

set.seed(11)

doRun <- FALSE
max_time <- 30*60
d <- 2
target <- get_multivariate_normal(d)
target$parameters$tau <- 5
nobservations <- 100
p <- 1
prefix <- ""

obsfile <- paste0(prefix, "mvnormaldata.d", d, ".n", nobservations, ".RData")
load(obsfile)

# function to simulate data
target$simulate <- function(theta){
  return(target$robservation(nobservations, theta, target$parameters))
}
# euclidean distance
# and also Euclidean distance
ground_p <- 2
deuclidean <- function(z){
  return(mean(apply(abs(z - obs), 2, function(v) (sum(v^ground_p))^(1/ground_p))^p)^(1/p))
}

# common algorithmic parameters
param_algo <- list(nthetas = 1024, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 1e5)

filename <- paste0(prefix, "mvnormalwsmc.d", d, ".n", nobservations, ".euclidean.RData")
results <- wsmc(deuclidean, target, param_algo, maxsimulation = 10^6, savefile = filename)
load(filename)
# results <- wsmc_continue(results, savefile = filename, maxtime = 10*60)

#
# load(filename)
# plot_threshold_time(results) + geom_point()
# mle <- rowMeans(obs)
# plot_bivariate(results, 1, 2, from = 10) + geom_vline(xintercept = mle[1]) + geom_hline(yintercept = mle[2])
# plot_marginal(results, 1, from = 10)
