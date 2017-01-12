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
# number of observations
n = 100
load(file = "~/Dropbox/ABCD/Results/data/gammadata.RData")
obs = obs[1:n]
sort_obs = sort(obs)
metricL1 <- function(xvec,yvec)  mean(abs(xvec - yvec))


#Choose m to be a multiple of n (or the other way around)
M = 500
N = c(1,5,10,20,50,100,500,1000)   #This is what we refer to as k in the paper
m = c(10,20,50,100,300,1000,5000,10000)

filename <- paste0("~/Dropbox/ABCD/Results/gamma_optim/gamma_optim_data_n",n,"_interplay_k_m.RData")

# t = proc.time()
# mewe_k_m = foreach(rep = 1:M) %dorng% {
#
#   #Store the data
#   mewe_mu_store = matrix(0,length(N),length(m))
#   mewe_sigma_store = matrix(0,length(N),length(m))
#   mewe_time_store = matrix(0,length(N),length(m))
#   count_evaluations_store = matrix(0,length(N),length(m))
#
#   #Generate all the randomness needed.
#   randomness = t(sapply(1:max(N), function(k) target$generate_randomness(max(m))))
#
#   #For each data set size m, find the MEWE using different numbers of sets of randomness N to approximate the expectation.
#   for(i in 1:length(m)){
#
#     #Subset all the sets of randomness to be of the required size (m), and sort the subsetted randomness.
#     use_randomness = randomness[,1:(m[i])]
#     sort_randomness = t(apply(use_randomness, MARGIN = 1, FUN = function(x) sort(x)))
#
#     #Make sure the observed and synthetic data are of the same length
#     if(m[i] > n){
#       sort_obs_mult = rep(sort_obs, each = m[i]/n)
#     } else{
#       sort_obs_mult = sort_obs
#     }
#     if(m[i] <= n){
#       sort_randomness_mult = t(apply(sort_randomness, MARGIN = 1, FUN = function(x) rep(x, each = n/m[i])))
#     } else{
#       sort_randomness_mult = sort_randomness
#     }
#
#     for(j in 1:length(N)){
#
#       #Subset data again to the required number of sets of randomness
#       if(N[j]==1){
#         sort_randomness_sub = t(as.matrix(sort_randomness_mult[1:N[j],]))
#       } else{
#         sort_randomness_sub = sort_randomness_mult[1:N[j],]
#       }
#
#       #Define the objective function defining the MEWE (choose number of sets of randomness).
#       mewe_objective = function(theta){
#         wass_dists = apply(sort_randomness_sub, MARGIN = 1 , function(x) metricL1(sort_obs_mult,(theta[2]*x+theta[1])))
#         out = mean(wass_dists)
#         return(out)
#       }
#
#       #Optimize the objective to find the MEWE.
#       init = c(runif(1,-4,4),runif(1,0.1,6))
#       optim_time = proc.time()
#       out = optim(init,mewe_objective)#,lower=c(-Inf,0.01),upper=c(Inf,Inf)) #Takes longer with specified domain.
#       optim_time = proc.time() - optim_time
#
#       mewe = out$par
#       objective_evals = out$count
#
#       #Store the results
#       mewe_mu_store[j,i] = mewe[1]
#       mewe_sigma_store[j,i] = mewe[2]
#       mewe_time_store[j,i] = optim_time[3]
#       count_evaluations_store[j,i] = objective_evals[1]
#
#     }
#   }
#
#   output = list(mewe_mu_store,mewe_sigma_store,mewe_time_store,count_evaluations_store)
#   return(output)
#
# }
# t = proc.time() - t
#
# #Use huge m to find MWE
# mm = 10^8
# randomness= target$generate_randomness(mm)
# sort_randomness = sort(randomness)
# sort_obs_mult = rep(sort_obs, each = mm/n)
# mewe_objective = function(theta){
#   out = metricL1(sort_obs_mult,(theta[2]*sort_randomness+theta[1]))
#   return(out)
# }
# init = c(runif(1,-4,4),runif(1,0.1,6))
# optim_time = proc.time()
# mwe = optim(init,mewe_objective)$par
# optim_time = proc.time() - optim_time
#
# save(mwe,mewe_k_m,file = filename)

load(filename)

df_mwek = data.frame(mu = mwe[1], sigma = mwe[2], k=1)
df_mwem = data.frame(mu = mwe[1], sigma = mwe[2], m=1)

mewe = lapply(1:length(m),
              function(i) lapply(1:length(N),
                                 function(j) t(sapply(1:M,
                                                      function(k) c(mewe_k_m[[k]][[1]][j,i],mewe_k_m[[k]][[2]][j,i])))))

df_mewe_list = list()
for(i in 1:length(m)){
  dummy = lapply(1:length(N), function(k) t(cbind(mewe[[i]][[k]],rep(k,M))))
  dummy = matrix(unlist(dummy), ncol = 3, byrow = TRUE)
  dummy = data.frame(dummy)
  names(dummy) = c("mu","sigma","k")
  df_mewe_list[[i]] = assign(paste("df_mewe_m",m[i],sep=""),dummy)
}


## Plot all levels of k for different m
#On the same axes
g = list()
for(i in 1:length(m)){
  g[[i]] <- ggplot(data = df_mewe_list[[i]], aes(x = mu, y = sigma, colour = k, group = k)) + ylim(0, 1.5) + xlim(1,3)
  g[[i]] <- g[[i]] + geom_point(alpha = 0.5)
  g[[i]] <- g[[i]] + scale_colour_gradient2(midpoint = 5) + theme(legend.position = "none", plot.title = element_text(size=10))
  g[[i]] <- g[[i]] + theme(axis.title.x=element_text(size=14), axis.text.x=element_text(size=12),axis.title.y=element_text(size=14), axis.text.y=element_text(size=12))
  g[[i]] <- g[[i]] + xlab(expression(mu)) + ylab(expression(sigma)) + ggtitle(paste("m = ",m[i],sep=""))
  g[[i]] <- g[[i]] + geom_point(data = df_mwek, aes(x = mu, y = sigma), color="black")
}
do.call(grid.arrange, c(g, ncol=4))

#On scaled axes
g = list()
for(i in 1:length(m)){
  g[[i]] <- ggplot(data = df_mewe_list[[i]] %>% filter(sigma> 0.2), aes(x = mu, y = sigma, colour = k, group = k)) #+ ylim(0, 1.5) + xlim(1,3)
  g[[i]] <- g[[i]] + geom_point(alpha = 0.5)
  g[[i]] <- g[[i]] + scale_colour_gradient2(midpoint = 5) + theme(legend.position = "none", plot.title = element_text(size=10))
  g[[i]] <- g[[i]] + theme(axis.title.x=element_text(size=14), axis.text.x=element_text(size=12),axis.title.y=element_text(size=14), axis.text.y=element_text(size=12))
  g[[i]] <- g[[i]] + xlab(expression(mu)) + ylab(expression(sigma)) + ggtitle(paste("m = ",m[i],sep=""))
  g[[i]] <- g[[i]] + geom_point(data = df_mwek, aes(x = mu, y = sigma), color="black")
}
do.call(grid.arrange, c(g, ncol=4))



### Stratify on k instead
mewe_k = lapply(1:length(N),
                function(j) lapply(1:length(m),
                                   function(i) t(sapply(1:M,
                                                        function(k) c(mewe_k_m[[k]][[1]][j,i],mewe_k_m[[k]][[2]][j,i])))))

df_mewe_k_list = list()
for(i in 1:length(N)){
  dummy = lapply(1:length(m), function(k) t(cbind(mewe_k[[i]][[k]],rep(k,M))))
  dummy = matrix(unlist(dummy), ncol = 3, byrow = TRUE)
  dummy = data.frame(dummy)
  names(dummy) = c("mu","sigma","m")
  df_mewe_k_list[[i]] = assign(paste("df_mewe_k",N[i],sep=""),dummy)
}


### Plot all levels of m for different k
#On the same axes
g = list()
for(i in 1:length(N)){
  g[[i]] <- ggplot(data = df_mewe_k_list[[i]] %>% filter(sigma> 0.2), aes(x = mu, y = sigma, colour = m, group = m)) + ylim(0.2, 1.4) + xlim(1.2,2.6)
  g[[i]] <- g[[i]] + geom_point(alpha = 0.5)
  g[[i]] <- g[[i]] + scale_colour_gradient2(midpoint = 5) + theme(legend.position = "none", plot.title = element_text(size=10))
  g[[i]] <- g[[i]] + theme(axis.title.x=element_text(size=14), axis.text.x=element_text(size=12),axis.title.y=element_text(size=14), axis.text.y=element_text(size=12))
  g[[i]] <- g[[i]] + xlab(expression(mu)) + ylab(expression(sigma)) + ggtitle(paste("k = ",N[i],sep=""))
  g[[i]] <- g[[i]] + geom_point(data = df_mwem, aes(x = mu, y = sigma), color="black")
}
do.call(grid.arrange, c(g, ncol=4))

#On scaled axes
g = list()
for(i in 1:length(N)){
  g[[i]] <- ggplot(data = df_mewe_k_list[[i]] %>% filter(sigma > 0.2), aes(x = mu, y = sigma, colour = m, group = m)) #+ ylim(0.2, 1.4) + xlim(1.2,2.6)
  g[[i]] <- g[[i]] + geom_point(alpha = 0.5)
  g[[i]] <- g[[i]] + scale_colour_gradient2(midpoint = 5) + theme(legend.position = "none", plot.title = element_text(size=10))
  g[[i]] <- g[[i]] + theme(axis.title.x=element_text(size=14), axis.text.x=element_text(size=12),axis.title.y=element_text(size=14), axis.text.y=element_text(size=12))
  g[[i]] <- g[[i]] + xlab(expression(mu)) + ylab(expression(sigma)) + ggtitle(paste("k = ",N[i],sep=""))
  g[[i]] <- g[[i]] + geom_point(data = df_mwem, aes(x = mu, y = sigma), color="black")
}
do.call(grid.arrange, c(g, ncol=4))


#Runtimes
average_runtimes = matrix(0,length(N),length(m))
for(i in 1:M){
  average_runtimes = mewe_k_m[[i]][[3]] + average_runtimes
}
average_runtimes = 1/M*average_runtimes

#Fixed k, runtime as function of m
plot(log(m),log(average_runtimes[4,]))
plot(log(m),log(average_runtimes[8,]))
plot(m,average_runtimes[8,])


#Fixed m, runtime as function of k
plot(log(N),log(average_runtimes[,4]))
plot(log(N),log(average_runtimes[,8]))
