library(winference)
rm(list = ls())
registerDoParallel(cores = detectCores())
setmytheme()

set.seed(11)

doRun <- FALSE
max_time <- 30*60
d <- 2
target <- get_multivariate_normal(d)
target$parameters$tau <- 5
nobservations <- 100
nparticles <- 2048
maxsim = 10^6
p <- 1
prefix <- ""

obsfile <- paste0(prefix, "mvnormaldata.d", d, ".n", nobservations, ".RData")
load(obsfile)

# function to simulate data
target$simulate <- function(theta){
  return(target$robservation(nobservations, theta, target$parameters))
}
# distance function
wdistance <- get_transport_to_y(obs, p = p)

filename <- paste0(prefix, "mvnormalrejection.d", d, ".n", nobservations, ".wasserstein.RData")

naccept = 2048
t = proc.time()
prior_theta = matrix(rnorm(d*maxsim,0,5), nrow=maxsim, ncol=d)
distances = rep(0,maxsim)
for(i in 1:maxsim){
  ysim = target$simulate(prior_theta[i,])
  {sink("/dev/null"); distances[i] = wdistance(ysim); sink();}
}
# distances = foreach(i = 1:maxsim, .combine = c) %dorng% {
#   ysim = target$simulate(prior_theta[i,])
#   return(wdistance(ysim))
# }
sort_distances = sort(distances)
results = prior_theta[distances <= sort_distances[naccept],]
t = proc.time() - t
t

save(results, t, file = filename)

