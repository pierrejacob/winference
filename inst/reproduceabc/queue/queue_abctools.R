library(abctools)
library(winference)
rm(list = ls())
registerDoParallel(cores = detectCores())
setmytheme()
set.seed(11)
target <- get_queue()

prefix = ""

nobservations <- 50
load(paste0(prefix, "50.intermediateobs.neal.RData"))
obs = matrix(obs, nrow = 1)


target$simulate <- function(theta){
  return(matrix(target$robservation(nobservations, theta, target$parameters), nrow = 1))
}


load(paste0(prefix, "queue_intermediate_wsmc_marginal.RData"))
results_was = results
#ndists = cumsum(results$ncomputed)
#nsim = ndists[34]
# nsim_approx = 10^6
# nsim = ndists[which.min(abs(nsim_approx - ndists))]
nsim = 10^7
m = 20 #determines number of order stats in initial summary stat
l = 1 #determines how many powers of the order stats are included in the initial summary stat

obs_subset = sort(obs)[c(1,(1:nobservations)[(nobservations/m)*(1:(m-2))+0.5],nobservations)]


############### Without constraint ###############

# t = proc.time()
# thetas = target$rprior(nsim, target$parameters)
#
# sumstats = apply(thetas, MARGIN = 1, function(theta) {
#   y_fake = target$simulate(theta)
#   sort_y_fake = sort(y_fake)
#   subset_y_fake = sort_y_fake[c(1,(1:nobservations)[(nobservations/m)*(1:(m-2))+0.5],nobservations)]
#   })
# sumstats = t(sumstats)
#
# # sumstats = foreach(i = 1:nsim, .combine = rbind) %dorng%{
# #     y_fake = target$simulate(theta)
# #     sort_y_fake = sort(y_fake)
# #     subset_y_fake = sort_y_fake[c(1,(1:nobservations)[(nobservations/m)*(1:(m-1))],nobservations)]
# #     return(subset_y_fake)
# # }

filename_sumstats = paste0("~/Documents/Harvard/ABCD/code/queueing/sumstats.queueing.abctools.n", nobservations,".m",m,".l",l,".nsim",nsim,".RData")
#filename_sumstats = paste0(prefix, "sumstats.queueing.abctools.n", nobservations,".m",m,".l",l,".nsim",nsim,".RData")
# savethis = list(sumstats,thetas)
# save(savethis, file = filename_sumstats)

load(filename_sumstats)
sumstats = savethis[[1]]
thetas = savethis[[2]]
#thetas[,2] = thetas[,2] + thetas[,1] #alternative approach
N = results_was$param_algo$nthetas
#tol = N/nsim
#tol = 1000/nsim
#tol = 500/nsim
tol = 100/nsim


#
# tfs <- list(function(x){cbind(x)})
# saabc <- semiauto.abc(obs = obs_subset, param = thetas, sumstats = sumstats,
#                       satr = tfs, overlap = TRUE, saprop = 1,
#                       abcprop = 1, tol = tol, method = "rejection",
#                       final.dens = TRUE)

#t = proc.time() - t
#results = list(thetas = saabc$post.sample, compute_time = t, threshold = tol, norder = m+1, npower = l, nsim = nsim)
filename = paste0(prefix, "queueing.abctools.n", nobservations,".m",m,".l",l,".nsim",nsim,".RData")
#save(results, file = filename)
load(filename)

# hist(saabc$post.sample[,1],prob=T)
# hist(saabc$post.sample[,2],prob=T)
# hist(saabc$post.sample[,3],prob=T)

#### Timing for experiment with 10^7 simulations ####
# > results$compute_time
# user   system  elapsed
# 3321.374  429.954 5616.018





# ############### With constraint ###############
#
# target$rprior <- function(ntheta, parameters){
#   theta1 <- runif(n = ntheta, min = 0, max = min(obs[1,]))
#   theta2minus1 <- runif(n = ntheta, min = 0, max = 10)
#   theta3 <- runif(n = ntheta, min = 0, max = 1/3)
#   return(cbind(theta1, theta2minus1, theta3))
# }
# #
# target$dprior <- function(thetas, parameters){
#   evals <- dunif(thetas[,1], min = 0, max = min(obs[1,]), log = TRUE)
#   evals <- evals + dunif(thetas[,2], min = 0, max = 10, log = TRUE)
#   evals <- evals + dunif(thetas[,3], min = 0, max = 1/3, log = TRUE)
#   return(evals)
# }
#
#
# # t = proc.time()
# # thetas = target$rprior(nsim, target$parameters)
# #
# # sumstats = apply(thetas, MARGIN = 1, function(theta) {
# #   y_fake = target$simulate(theta)
# #   sort_y_fake = sort(y_fake)
# #   subset_y_fake = sort_y_fake[c(1,(1:nobservations)[(nobservations/m)*(1:(m-2))+0.5],nobservations)]
# #   })
# # sumstats = t(sumstats)
# #
# # # sumstats = foreach(i = 1:nsim, .combine = rbind) %dorng%{
# # #     y_fake = target$simulate(theta)
# # #     sort_y_fake = sort(y_fake)
# # #     subset_y_fake = sort_y_fake[c(1,(1:nobservations)[(nobservations/m)*(1:(m-1))],nobservations)]
# # #     return(subset_y_fake)
# # # }
#
# filename_sumstats = paste0("~/Documents/Harvard/ABCD/Code/queueing/sumstats.constrained.queueing.abctools.n", nobservations,".m",m,".l",l,".nsim",nsim,".RData")
# #savethis = list(sumstats,thetas)
# #save(savethis, file = filename_sumstats)
#
# load(filename_sumstats)
# sumstats = savethis[[1]]
# thetas = savethis[[2]]
# #N = results_was$param_algo$nthetas
# N = 2048
# tol = N/nsim
# #tol = 1000/nsim
# #tol = 500/nsim
# #tol = 100/nsim
#
#
#
# tfs <- list(function(x){cbind(x)})
# saabc <- semiauto.abc(obs = obs_subset, param = thetas, sumstats = sumstats,
#                       satr = tfs, overlap = TRUE, saprop = 1,
#                       abcprop = 1, tol = tol, method = "rejection",
#                       final.dens = TRUE)
#
# #t = proc.time() - t
# results = list(thetas = saabc$post.sample, compute_time = t, threshold = tol, norder = m+1, npower = l, nsim = nsim)
# filename = paste0(prefix, "constrained.queueing.abctools.n", nobservations,".m",m,".l",l,".nsim",nsim,".RData")
# save(results, file = filename)
# #load(filename)
#
