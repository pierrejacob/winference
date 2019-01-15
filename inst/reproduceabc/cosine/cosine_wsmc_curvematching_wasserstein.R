library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()
set.seed(11)

prefix = ""

target <- get_cosine()

nobservations <- 100

load(paste0(prefix, "cosinedata.RData"))
obs <- matrix(obs[1:nobservations], nrow = 1)
target$simulate <- function(theta){
  return(matrix(target$robservation(nobservations, theta, target$parameters, target$generate_randomness(nobservations)), nrow = 1))
}

param_algo <- list(nthetas = 2048, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 1e5)

lambda <- 1
multiplier <- lambda*(max(obs[1,]) - min(obs[1,]))
augment <- function(series) rbind(multiplier * (1:length(series))/length(series), series)
augmented_obs <- augment(obs)
plot(augmented_obs[1,],augmented_obs[2,])
sorted_augmented_obs <- augmented_obs[,hilbert_order(augmented_obs)]

compute_d <- function(y_fake){
  augmented_y_fake <- augment(y_fake)
  sink("/dev/null")
  z <- exact_transport_distance(augmented_obs, augmented_y_fake, p = 1, ground_p = 2)
  sink(NULL)
  return(z)
}

# library(microbenchmark)
# microbenchmark(y_sim <- target$simulate(true_theta),compute_d(y_sim))
# y_sim <- target$simulate(true_theta)
# compute_d(y_sim)

filename <- paste0(prefix, "cosine_wsmc_curvematching.wasserstein.lambda", lambda, ".n", nobservations, ".RData")
# results <- wsmc(compute_d, target, param_algo, savefile = filename, maxsimulation = 1e6)
load(filename)
# results$compute_d = compute_d
# ncomp = sum(results$ncomputed)
# results <- wsmc_continue(results, savefile = filename, maxsimulation = 1e6)
