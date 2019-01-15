library(microbenchmark)
library(winference)
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

# wasserstein distance
wdistance <- get_transport_to_y(obs, p = p)

# euclidean distance
ground_p <- 2
deuclidean <- function(z){
  return(mean(apply(abs(z - obs), 2, function(v) (sum(v^ground_p))^(1/ground_p))^p)^(1/p))
}

# distance between summary
summary_obs <- rowMeans(obs)
dsummary <- function(z){
  summary_z <- rowMeans(z)
  return(mean(abs(summary_z - summary_obs)))
}

{sink("/dev/null"); wtime = invisible(microbenchmark(wdistance(target$simulate(true_theta)), times = 10000)); sink();}
wtime

etime = microbenchmark(deuclidean(target$simulate(true_theta)), times = 10000)
etime

stime = microbenchmark(dsummary(target$simulate(true_theta)), times = 10000)
stime

simtime = microbenchmark(target$simulate(true_theta), times = 10000)
simtime

# # #result from timing (average values in seconds)
# wtime = 0.0123
# etime = 0.000380
# stime = 0.000064
# simtime = 0.000040
#
# wtime/etime
# wtime/stime
#
# #load simulation results
# filename <- paste0(prefix, "mvnormalwsmc.d", d, ".n", nobservations, ".wasserstein.RData")
# load(filename)
# results_wasserstein <- results
# tail(results_wasserstein$compute_times,n=1)
#
#
# filename <- paste0(prefix, "mvnormalwsmc.d", d, ".n", nobservations, ".euclidean.RData")
# load(filename)
# results_euclidean <- results
# tail(results_euclidean$compute_times,n=1)
#
#
# filename <- paste0(prefix, "mvnormalwsmc.d", d, ".n", nobservations, ".summary.RData")
# load(filename)
# results_summary <- results
# tail(results_summary$compute_times,n=1)
#
#
# filename <- paste0(prefix, "mvnormalrejection.d", d, ".n", nobservations, ".summary.RData")
# load(filename)
# results_rejection_summary <- results
# t

