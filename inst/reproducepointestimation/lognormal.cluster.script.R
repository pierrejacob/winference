## computes MEWE on B=1000 bootstrap version of a data set
## File to be run on a cluster preferably!

rm(list = ls())
library(parallel)
# parallel RNG using L'Ecuyer et al (2002)
RNGkind("L'Ecuyer-CMRG") # L'Ecuyer CMRG required for multiple streams
igrid <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(1) # initial seed
for (i in 1:igrid){
  .Random.seed <- nextRNGStream(.Random.seed) # compute appropriate stream
}

source("lognormal_functions.R")

#Pick m to be larger than or equal to max(n) and a multiple of each entry in n 
#M = 50    #Number of data sets from dgp
n = 10^3  #Number of samples in each data set from dgp
N = 20    #Number of synthetic data sets for finding MEWE
m = 10^4  #Number of samples in the synthetic data sets
B = 1000  #Number of bootstrap samples per data set from dgp
alpha = 0.05 #Desired significance level, compute the alpha/2 and 1-alpha/2 quantiles

filename = paste0("output/lognormal_bootstrap_run",igrid,"_N",N,"_m",m,"_n",n,"_B",B,".Rdata")

# Bootstrap confidence intervals for the A and B parameters
time_elapsed = proc.time()
#Temporary storage for the bootstrapped estimators and other information
bs_store = matrix(0, B, thetadim)#target$thetadim)
bs_count_store = rep(0, B)

#Generate the observations, sort, and replicate observations to match length of the synthetic data
obs_rand = target$generate_randomness(n)
obs = target$robservation(true_theta,obs_rand)
sort_obs = sort(obs)
sort_obs_mult = rep(sort_obs, each = m/n)

#Generate the synthetic randomness, sort. 
randomness = t(sapply(1:N, function(x) target$generate_randomness(m)))
# sort_randomness = t(apply(randomness, 1, sort))

#Define the objective to be minimized to find the MEWE
obj1 = function(theta){
  if(theta[2] < 0){
    out = 100
  } else{
    wass_dists = apply(randomness, 1, function(x) metricL1(sort_obs_mult, sort(target$robservation(theta, x))))
    out = mean(wass_dists)
  }
  return(out)
}

#Find the MEWE
mewe = optim(true_theta,obj1)
mewe$par

# #### Bootstrap
#Generate multinomials
multinom = matrix(sample(1:n, size = n*B, replace = T), nrow = B)
sort_multinom = t(apply(multinom, 1, sort))

#Resample B times from the empirical distribution, sorted
emp_sort_obs = t(apply(sort_multinom, 1, function(x) sort_obs[x]))

#For each of the resampled data sets, find the MEWE.
for(i in 1:B){
  #Replicate observations to match length of the synthetic data
  emp_obs_mult = rep(emp_sort_obs[i,], each = m/n)
  
  #Define the objective function
  obj_bs = function(theta){
    if(theta[2] < 0){
      out = 100
    } else{
      wass_dists = apply(randomness, 1, function(x) metricL1(emp_obs_mult, sort(target$robservation(theta, x))))
      out = mean(wass_dists)
    }
    return(out)
  }
  
  #Find minimum, initializing at the MEWE
  mewe_bs = try(optim(mewe$par,obj_bs))
  if (inherits(mewe_bs, "try-error")){
    mewe_bs <- list(par = rep(NA, thetadim), count = Inf)
  }
  
  #Store the stuff
  bs_store[i,] = mewe_bs$par
  bs_count_store[i] = mewe_bs$count[1]
}

time_elapsed = proc.time() - time_elapsed

save(bs_store, bs_count_store, time_elapsed, file = filename)
