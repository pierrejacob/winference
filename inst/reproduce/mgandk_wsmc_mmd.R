library(winference)
registerDoParallel(cores = 6)
rm(list = ls())
setmytheme()
set.seed(11)
target <- get_mgandk()
prefix <- ""

nobservations <- 500
load(paste0(prefix, "mgandkdata.RData"))
obs <- obs[,1:nobservations]
target$simulate <- function(theta) target$robservation(nobservations, theta)

compute_d <- get_mmd_to_y(obs)
y_sim <- target$simulate(target$rprior(1, target$parameters)[1,])
compute_d(y_sim)

param_algo <- list(nthetas = 1024, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 100000)

filename <- paste0(prefix, "mgandk.wsmc.n", nobservations, ".mmd.RData")
results <- wsmc(compute_d, target, param_algo, savefile = filename, maxsimulation = 1e6)
load(filename)
# results <- wsmc_continue(results, savefile = filename, maxsimulation = 1e6)
# plot_marginal_time(results, 9, from_step = 50)
# qplot(x = cumsum(results$ncomputed), y = results$threshold_history, geom = "line") + scale_x_log10() + scale_y_log10()
# qplot(x = tail(results$distances_history, 1)[[1]], geom = "histogram")

