library(doParallel)
registerDoParallel(cores = detectCores())
rm(list = ls())
set.seed(11)

source("gamma_normal_functions.R")
gen_obs_data = function(n){rgamma(n,10,5)}

#Pick m to be larger than or equal to max(n) and a multiple of each entry in n
M = 400   #Number of data sets from dgp
n = 1000   #Number of samples in each data set from dgp
N = 20    #Number of synthetic data sets for finding MEWE
m = 10^4  #Number of samples in the synthetic data sets
B = 1000  #Number of bootstrap samples per data set from dgp
alpha = 0.05 #Desired significance level, compute the alpha/2 and 1-alpha/2 quantiles

prefix = ""
filename = paste0(prefix,"gamma_normal_bootstrap_M",M,"_N",N,"_m",m,"_n",n,".Rdata")

# #Calculate the "true" parameter once.
# nt = 10^8
# mt = 10^8
#
# obs = gen_obs_data(nt)
# randomness = target$generate_randomness(mt)
#
# sort_obs = sort(obs)
# sort_randomness = sort(randomness)
#
# obj1 = function(theta){
#   wass_dist = metricL1(sort_obs,(theta[2]*sort_randomness+theta[1]))
#   return(wass_dist)
# }
#
# init = c(2,0.6)
# true_param = optim(init,obj1)  #Equal to c(1.9641217,0.6209369)
true_param =  c(1.9641217,0.6209369)

#Bootstrap confidence intervals

mewe_bootstrap = foreach(rep = 1:M) %dorng% {

  #Temporary storage for the bootstrapped estimators and other information
  bs_store = matrix(0,B,target$thetadim)
  bs_count_store = rep(0,B)


  #Generate the observations, sort, and replicate observations to match length of the synthetic data
  obs = gen_obs_data(n)
  sort_obs = sort(obs)
  sort_obs_mult = rep(sort_obs, each = m/n)

  #Generate the synthetic randomness, sort
  randomness = t(sapply(1:N, function(x) target$generate_randomness(m)))
  sort_randomness = t(apply(randomness, 1, sort))

  #Define the objective to be minimized to find the MEWE
  obj1 = function(theta){
    wass_dists = apply(sort_randomness, 1, function(x) metricL1(sort_obs_mult,(theta[2]*x+theta[1])))
    out = mean(wass_dists)
    return(out)
  }

  #Find the MEWE
  init = c(runif(1,-3,5),runif(1,0.1,2))
  mewe = optim(init,obj1)

  #### Bootstrap

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
      wass_dists = apply(sort_randomness, 1, function(x) metricL1(emp_obs_mult,(theta[2]*x+theta[1])))
      out = mean(wass_dists)
      return(out)
    }

    #Find minimum, initializing at the MEWE
    mewe_bs = optim(mewe$par,obj_bs)

    #Store the stuff
    bs_store[i,] = mewe_bs$par
    bs_count_store[i] = mewe_bs$count[1]

  }

  #Compute the empirical quantiles of the bootstrap samples, both separately for each parameter and with Bonferroni correction.
  q_mu = quantile(bs_store[,1],probs = c(alpha/4,alpha/2,1-alpha/2,1-alpha/4))
  q_sigma = quantile(bs_store[,2],probs = c(alpha/4,alpha/2,1-alpha/2,1-alpha/4))

  ci_mu = q_mu[2:3]
  ci_sigma = q_sigma[2:3]
  ci_both = rbind(q_mu[c(1,4)],q_sigma[c(1,4)])

  #Summary of the count distribution
  counts = summary(bs_count_store)

  #Summarize output in a list
  out = list(mewe$par,ci_mu,ci_sigma,ci_both)
  return(out)
}

save(mewe_bootstrap, file = filename)


load(filename)

ci_mu = t(sapply(mewe_bootstrap, function(x) x[[2]]))
mu_star = true_param[1]
cover_mu = (ci_mu[,1] <= mu_star) & (mu_star <= ci_mu[,2])
coverage_mu = sum(cover_mu)/M
coverage_mu

ci_sigma = t(sapply(mewe_bootstrap, function(x) x[[3]]))
sigma_star = true_param[2]
cover_sigma = (ci_sigma[,1] <= sigma_star) & (sigma_star <= ci_sigma[,2])
coverage_sigma = sum(cover_sigma)/M
coverage_sigma

ci_both_mu = t(sapply(mewe_bootstrap, function(x) x[[4]][1,]))
ci_both_sigma = t(sapply(mewe_bootstrap, function(x) x[[4]][2,]))
cover_both = (ci_both_mu[,1] <= mu_star) & (mu_star <= ci_both_mu[,2]) &
  (ci_both_sigma[,1] <= sigma_star) & (sigma_star <= ci_both_sigma[,2])
coverage_both = sum(cover_both)/M
coverage_both
