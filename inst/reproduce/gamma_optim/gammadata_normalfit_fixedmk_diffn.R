#+ presets, echo = FALSE, warning = FALSE, message = FALSE
library(winference)
library(ggplot2)
library(gridExtra)
library(ggthemes)
library(dplyr)
library(foreach)
library(doMC)
library(doRNG)
library(reshape2)
registerDoMC(cores = 4)
rm(list = ls())
setmytheme()
set.seed(11)


target <- get_normal()
gen_obs_data = function(n){rgamma(n,10,5)}
metricL1 <- function(xvec,yvec)  mean(abs(xvec - yvec))
metricL2 <- function(xvec,yvec)  sqrt(mean((xvec - yvec)^2))

#Pick m to be larger than or equal to max(n) and a multiple of each entry in n 
M = 1000
N = 20
m = 10^4
n = c(50,100,250,500,1000,5000,10000)


######### Fits both location and scale. See bottom of file for location only.

filename <- paste0("~/Dropbox/ABCD/Results/gamma_optim/gamma_optim_data_m",m,"_k",N,"_bound",bound,".RData")

t = proc.time()
mewe_gamma = foreach(rep = 1:M) %dorng% {

  #Allocate space for output
  mewe1_store = matrix(0,length(n),target$thetadim)
  mewe1_runtimes = rep(0,length(n))
  mewe1_evals = rep(0,length(n))

  mewe2_store = matrix(0,length(n),target$thetadim)
  mewe2_runtimes = rep(0,length(n))
  mewe2_evals = rep(0,length(n))

  #generate all observations and sets of randomness to be used
  obs_all = gen_obs_data(max(n))
  sort_randomness = t(sapply(1:N, function(x) sort(target$generate_randomness(m))))

  for(i in 1:(length(n)) ){
    #Subset observations and sort
    obs = obs_all[1:n[i]]
    sort_obs = sort(obs)
    sort_obs_mult = rep(sort_obs, each = m/n[i])

    #Define optimization objectives
    mewe1_obj = function(theta){
      wass_dists = apply(sort_randomness, MARGIN = 1 , function(x) metricL1(sort_obs_mult,(theta[2]*x+theta[1])))
      out = mean(wass_dists)
      return(out)
    }
    mewe2_obj = function(theta){
      wass_dists = apply(sort_randomness, MARGIN = 1 , function(x) metricL2(sort_obs_mult,(theta[2]*x+theta[1])))
      out = mean(wass_dists)
      return(out)
    }

    #Initial point for optimization
    init = c(runif(1,-3,5),runif(1,0.1,2))

    #Optimization
    t_mewe1 = proc.time()
    mewe1 = optim(init,mewe1_obj)
    t_mewe1 = proc.time() - t_mewe1

    t_mewe2 = proc.time()
    mewe2 = optim(init,mewe2_obj)
    t_mewe2 = proc.time() - t_mewe2

    #Save the results
    mewe1_store[i,] = mewe1$par
    mewe1_runtimes[i] = t_mewe1[3]
    mewe1_evals[i] = mewe1$count

    mewe2_store[i,] = mewe2$par
    mewe2_runtimes[i] = t_mewe2[3]
    mewe2_evals[i] = mewe2$count
  }
  #Format into list
  mewe1_list = list(mewe1_store,mewe1_runtimes,mewe1_evals)
  mewe2_list = list(mewe2_store,mewe2_runtimes,mewe2_evals)

  #Output
  list(mewe1_list,mewe2_list)
}
t = proc.time() - t

save(mewe_gamma,file = filename)


# mewe1 = t(sapply(mewe_cauchy, function(x) x[[1]][[1]]))
# mewe2bdd = t(sapply(mewe_cauchy, function(x) x[[3]][[1]]))
# 
# mewe1_dummy = lapply(1:length(n), function(k) t(cbind(mewe1[,k],rep(n[k],M))))
# mewe1_dummy = matrix(unlist(mewe1_dummy), ncol = 2, byrow = TRUE)
# df_mewe1 = data.frame(mewe1_dummy)
# names(df_mewe1) = c("mu","n")
# df_mewe1$n = as.factor(df_mewe1$n)

