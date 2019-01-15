library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()

set.seed(11)

doRun <- FALSE
max_time <- 30*60
d <- 2
target <- get_multivariate_normal(d)
target$parameters$tau <- 5
nobservations <- 100
nparticles <- 2048
p <- 1
prefix <- ""

obsfile <- paste0(prefix, "mvnormaldata.d", d, ".n", nobservations, ".RData")
load(obsfile)

# function to simulate data
target$simulate <- function(theta){
  return(target$robservation(nobservations, theta, target$parameters))
}
# wasserstein distance
wdistance <- get_transport_to_y(obs, p = p)
#
param_algo <- list(nthetas = nparticles, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 1e5)

filename <- paste0(prefix, "mvnormalwsmc.d", d, ".n", nobservations, ".wasserstein.RData")
results <- wsmc(wdistance, target, param_algo, maxsimulation = 10^6, savefile = filename)
load(filename)
# results <- wsmc_continue(results, savefile = filename, maxsimulation = 800000)
#
# load(filename)
# plot_threshold_time(results) + geom_point()
# mle <- rowMeans(obs)
# plot_bivariate(results, 1, 2, from = 10) + geom_vline(xintercept = mle[1]) + geom_hline(yintercept = mle[2])
# plot_marginal(results, 1, from = 10)

# library(microbenchmark)
# microbenchmark(wdistance(target$simulate(true_theta)), times = 1000)
