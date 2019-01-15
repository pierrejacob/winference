#+ presets, echo = FALSE, warning = FALSE, message = FALSE
library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()
set.seed(11)
target <- get_gandk()
# number of observations
nobservations <- 250
prefix = ""
load(paste0(prefix, "gandkdata.RData"))
obs <- obs[1:nobservations]
sort_obs = sort(obs)


#compute_d <- get_hilbert_to_y(matrix(obs, nrow = 1))

compute_d = function(y){
  sort_y = sort(y)
  mean(abs(sort_y-sort_obs))
}

target$simulate <- function(theta){
  return(matrix(target$robservation(nobservations, theta, target$parameters, target$generate_randomness(nobservations)), nrow = 1))
}

# M=10000
# ts = rep(0,M)
# for(i in 1:M){
#   y_sim = target$simulate(target$rprior(1, target$parameters))
#   t=proc.time()
#   compute_d_alt(y_sim)
#   t=proc.time() -t
#   ts[i] = t[3]
# }

param_algo <- list(nthetas = 2048, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 100000)
#t = proc.time()
filename <- paste0(prefix, "gandkwsmc.n", nobservations, ".RData")
results <- wsmc(compute_d, target, param_algo, savefile = filename, maxsim = 2.4e6)
#t = proc.time() - t
load(filename)
#results <- wsmc_continue(results, savefile = filename, maxtime = 14*60*60)


# load(filename)
# wsmc.df <- wsmc_to_dataframe(results)
# nsteps <- max(wsmc.df$step)
#
# # plot_bivariate_polygon(results, 1, 2)
# # plot_bivariate_polygon(results, 3, 4)
#
# library(gridExtra)
# grid.arrange(plot_marginal_time(results, 1),
#   plot_marginal_time(results, 2),
#   plot_marginal_time(results, 3),
#   plot_marginal_time(results, 4), nrow = 2)


