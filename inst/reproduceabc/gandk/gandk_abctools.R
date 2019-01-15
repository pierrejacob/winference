library(abctools)
library(winference)
rm(list = ls())
registerDoParallel(cores = detectCores())
setmytheme()
set.seed(11)
target <- get_gandk()

prefix = ""

nobservations <- 250
load(paste0(prefix, "gandkdata.RData"))
obs <- obs[1:nobservations]

target$simulate <- function(theta){
  return(matrix(target$robservation(nobservations, theta, target$parameters, target$generate_randomness(nobservations)), nrow = 1))
}


# load(paste0(prefix, "gandkwsmc.n", nobservations, ".RData"))
# results_was = results
# ndists = cumsum(results$ncomputed)
# nsim = ndists[34]

# nsim_approx = 10^6
# nsim = ndists[which.min(abs(nsim_approx - ndists))]

nsim = 2.4*10^6

m = 10 #determines number of order stats in initial summary stat
l = 4 #determines how many powers of the order stats are included in the initial summary stat

t = proc.time()
thetas = target$rprior(nsim, target$parameters)

obs_subset = sort(obs)[c(1,(1:nobservations)[(nobservations/m)*(1:(m-1))],nobservations)]

sumstats = apply(thetas, MARGIN = 1, function(theta) {
  y_fake = target$simulate(theta)
  sort_y_fake = sort(y_fake)
  subset_y_fake = sort_y_fake[c(1,(1:nobservations)[(nobservations/m)*(1:(m-1))],nobservations)]
  })
sumstats = t(sumstats)

# sumstats = foreach(i = 1:nsim, .combine = rbind) %dorng%{
#     y_fake = target$simulate(theta)
#     sort_y_fake = sort(y_fake)
#     subset_y_fake = sort_y_fake[c(1,(1:nobservations)[(nobservations/m)*(1:(m-1))],nobservations)]
#     return(subset_y_fake)
# }

filename_sumstats = paste0(prefix, "sumstats.gandk.abctools.n", nobservations,".m",m,".l",l,".nsim",nsim,".RData")
save(list(sumstats,thetas), file = filename_sumstats)

load(filename_sumstats)
N = results_was$param_algo$nthetas
tol = N/nsim
#tol = 1000/nsim
#tol = 500/nsim
#tol = 100/nsim


tfs <- list(function(x){cbind(x, x^2, x^3, x^4)})
saabc <- semiauto.abc(obs = obs_subset, param = thetas, sumstats = sumstats,
                      satr = tfs, overlap = TRUE, saprop = 1,
                      abcprop = 1, tol = tol, method = "rejection",
                      final.dens = TRUE)

t = proc.time() - t
t
results = list(thetas = saabc$post.sample, compute_time = t, threshold = tol, norder = m+1, npower = l, nsim = nsim)
filename = paste0(prefix, "gandk.abctools.n", nobservations,".m",m,".l",l,".nsim",nsim,".RData")
save(results, file = filename)
load(filename)

# hist(saabc$post.sample[,1],prob=T)
# hist(saabc$post.sample[,2],prob=T)
# hist(saabc$post.sample[,3],prob=T)
# hist(saabc$post.sample[,4],prob=T)



