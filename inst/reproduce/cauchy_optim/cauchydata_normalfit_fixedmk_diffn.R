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
gen_obs_data = function(n){rcauchy(n)}
metricL1 <- function(xvec,yvec)  mean(abs(xvec - yvec))
metricL2 <- function(xvec,yvec)  sqrt(mean((xvec - yvec)^2))
metricBDD = function(xvec,yvec,metric,bound) min(bound,metric(xvec,yvec))

#Pick m to be larger than or equal to max(n) and a multiple of each entry in n 
M = 1000
N = 20
m = 10^4
n = c(50,100,250,500,1000,5000,10000)
bound = 10


######### Fits both location and scale. See bottom of file for location only.

filename <- paste0("~/Dropbox/ABCD/Results/cauchy_optim/cauchy_ls_optim_data_m",m,"_k",N,"_bound",bound,".RData")

# t = proc.time()
# mewe_cauchy = foreach(rep = 1:M) %dorng% {
# 
#   #Allocate space for output
#   mewe1_store = matrix(0,length(n),target$thetadim)
#   mewe1_runtimes = rep(0,length(n))
#   mewe1_evals = rep(0,length(n))
# 
#   mewe1bdd_store = matrix(0,length(n),target$thetadim)
#   mewe1bdd_runtimes = rep(0,length(n))
#   mewe1bdd_evals = rep(0,length(n))
# 
#   mewe2bdd_store = matrix(0,length(n),target$thetadim)
#   mewe2bdd_runtimes = rep(0,length(n))
#   mewe2bdd_evals = rep(0,length(n))
# 
#   #generate all observations and sets of randomness to be used
#   obs_all = gen_obs_data(max(n))
#   sort_randomness = t(sapply(1:N, function(x) sort(target$generate_randomness(m))))
# 
#   for(i in 1:(length(n)) ){
#     #Subset observations and sort
#     obs = obs_all[1:n[i]]
#     sort_obs = sort(obs)
#     sort_obs_mult = rep(sort_obs, each = m/n[i])
# 
#     #Define optimization objectives
#     mewe1_obj = function(theta){
#       wass_dists = apply(sort_randomness, MARGIN = 1 , function(x) metricL1(sort_obs_mult,(theta[2]*x+theta[1])))
#       out = mean(wass_dists)
#       return(out)
#     }
#     mewe1bdd_obj = function(theta){
#       wass_dists = apply(sort_randomness, MARGIN = 1 , function(x) metricBDD(sort_obs_mult,(theta[2]*x+theta[1]),metricL1,bound))
#       out = mean(wass_dists)
#       return(out)
#     }
#     mewe2bdd_obj = function(theta){
#       wass_dists = apply(sort_randomness, MARGIN = 1 , function(x) metricBDD(sort_obs_mult,(theta[2]*x+theta[1]),metricL2,bound))
#       out = mean(wass_dists)
#       return(out)
#     }
# 
#     #Initial point for optimization
#     init = c(runif(1,-5,5),runif(1,0.2,3))
# 
#     #Optimization
#     t_mewe1 = proc.time()
#     mewe1 = optim(init,mewe1_obj)
#     t_mewe1 = proc.time() - t_mewe1
# 
#     t_mewe1bdd = proc.time()
#     mewe1bdd = optim(init,mewe1bdd_obj)
#     t_mewe1bdd = proc.time() - t_mewe1bdd
# 
#     t_mewe2bdd = proc.time()
#     mewe2bdd = optim(init,mewe2bdd_obj)
#     t_mewe2bdd = proc.time() - t_mewe2bdd
# 
#     #Save the results
#     mewe1_store[i,] = mewe1$par
#     mewe1_runtimes[i] = t_mewe1[3]
#     mewe1_evals[i] = mewe1$count
# 
#     mewe1bdd_store[i,] = mewe1bdd$par
#     mewe1bdd_runtimes[i] = t_mewe1bdd[3]
#     mewe1bdd_evals[i] = mewe1bdd$count
# 
#     mewe2bdd_store[i,] = mewe2bdd$par
#     mewe2bdd_runtimes[i] = t_mewe2bdd[3]
#     mewe2bdd_evals[i] = mewe2bdd$count
#   }
#   #Format into list
#   mewe1_list = list(mewe1_store,mewe1_runtimes,mewe1_evals)
#   mewe1bdd_list = list(mewe1bdd_store,mewe1bdd_runtimes,mewe1bdd_evals)
#   mewe2bdd_list = list(mewe2bdd_store,mewe2bdd_runtimes,mewe2bdd_evals)
# 
#   #Output
#   list(mewe1_list,mewe1bdd_list,mewe2bdd_list)
# }
# t = proc.time() - t
# 
# save(mewe_cauchy,file = filename)
# 
# 
# mewe1 = t(sapply(mewe_cauchy, function(x) x[[1]][[1]]))
# mewe1bdd = t(sapply(mewe_cauchy, function(x) x[[2]][[1]]))
# mewe2bdd = t(sapply(mewe_cauchy, function(x) x[[3]][[1]]))
# 
# mewe1_dummy = lapply(1:length(n), function(k) t(cbind(mewe1[,k],rep(n[k],M))))
# mewe1_dummy = matrix(unlist(mewe1_dummy), ncol = 2, byrow = TRUE)
# df_mewe1 = data.frame(mewe1_dummy)
# names(df_mewe1) = c("mu","n")
# df_mewe1$n = as.factor(df_mewe1$n)









######### Fits location only

#filename <- paste0("~/Dropbox/ABCD/Results/cauchy_optim/cauchy_optim_data_m",m,"_k",N,"_bound",bound,".RData")
#
# t = proc.time()
# mewe_cauchy = foreach(rep = 1:M) %dorng% {
#   
#   #Allocate space for output
#   mewe1_store = rep(0,length(n))
#   mewe1_runtimes = rep(0,length(n))
#   mewe1_evals = rep(0,length(n))
#   
#   mewe1bdd_store = rep(0,length(n))
#   mewe1bdd_runtimes = rep(0,length(n))
#   mewe1bdd_evals = rep(0,length(n))
#   
#   mewe2bdd_store = rep(0,length(n))
#   mewe2bdd_runtimes = rep(0,length(n))
#   mewe2bdd_evals = rep(0,length(n))
#   
#   #generate all observations and sets of randomness to be used
#   obs_all = gen_obs_data(max(n))
#   sort_randomness = t(sapply(1:N, function(x) sort(target$generate_randomness(m))))
#   
#   for(i in 1:(length(n)) ){
#     #Subset observations and sort
#     obs = obs_all[1:n[i]]
#     sort_obs = sort(obs)
#     sort_obs_mult = rep(sort_obs, each = m/n[i])
#     
#     #Define optimization objectives
#     mewe1_obj = function(theta){
#       wass_dists = apply(sort_randomness, MARGIN = 1 , function(x) metricL1(sort_obs_mult,(x+theta)))
#       out = mean(wass_dists)
#       return(out)
#     }
#     mewe1bdd_obj = function(theta){
#       wass_dists = apply(sort_randomness, MARGIN = 1 , function(x) metricBDD(sort_obs_mult,(x+theta),metricL1,bound))
#       out = mean(wass_dists)
#       return(out)
#     }
#     mewe2bdd_obj = function(theta){
#       wass_dists = apply(sort_randomness, MARGIN = 1 , function(x) metricBDD(sort_obs_mult,(x+theta),metricL2,bound))
#       out = mean(wass_dists)
#       return(out)
#     }
#     
#     #Initial point for optimization
#     init = runif(1,-5,5)
#     
#     #Optimization
#     t_mewe1 = proc.time()
#     mewe1 = optim(init,mewe1_obj)
#     t_mewe1 = proc.time() - t_mewe1
#     
#     t_mewe1bdd = proc.time()
#     mewe1bdd = optim(init,mewe1bdd_obj)
#     t_mewe1bdd = proc.time() - t_mewe1bdd
#     
#     t_mewe2bdd = proc.time()
#     mewe2bdd = optim(init,mewe2bdd_obj)
#     t_mewe2bdd = proc.time() - t_mewe2bdd
#     
#     #Save the results
#     mewe1_store[i] = mewe1$par
#     mewe1_runtimes[i] = t_mewe1[3]
#     mewe1_evals[i] = mewe1$count
#     
#     mewe1bdd_store[i] = mewe1bdd$par
#     mewe1bdd_runtimes[i] = t_mewe1bdd[3]
#     mewe1bdd_evals[i] = mewe1bdd$count
#     
#     mewe2bdd_store[i] = mewe2bdd$par
#     mewe2bdd_runtimes[i] = t_mewe2bdd[3]
#     mewe2bdd_evals[i] = mewe2bdd$count
#   }
#   #Format into list
#   mewe1_list = list(mewe1_store,mewe1_runtimes,mewe1_evals)
#   mewe1bdd_list = list(mewe1bdd_store,mewe1bdd_runtimes,mewe1bdd_evals)
#   mewe2bdd_list = list(mewe2bdd_store,mewe2bdd_runtimes,mewe2bdd_evals)
#   
#   #Output
#   list(mewe1_list,mewe1bdd_list,mewe2bdd_list)
# }
# t = proc.time() - t
# 
# save(mewe_cauchy,file = filename)
# 
# 
# mewe1 = t(sapply(mewe_cauchy, function(x) x[[1]][[1]]))
# mewe1bdd = t(sapply(mewe_cauchy, function(x) x[[2]][[1]]))
# mewe2bdd = t(sapply(mewe_cauchy, function(x) x[[3]][[1]]))
# 
# mewe1_dummy = lapply(1:length(n), function(k) t(cbind(mewe1[,k],rep(n[k],M))))
# mewe1_dummy = matrix(unlist(mewe1_dummy), ncol = 2, byrow = TRUE)
# df_mewe1 = data.frame(mewe1_dummy)
# names(df_mewe1) = c("mu","n")
# df_mewe1$n = as.factor(df_mewe1$n)
# 
# p = ggplot(df_mewe1, aes(x=n, y=mu))
# p = p + geom_point(alpha=0.05)
# p = p + theme(plot.title = element_text(size = rel(0.6)))
# p = p + xlab("n") + ylab(expression(mu))
# p = p + theme(axis.title.x=element_text(size=a.title.x), axis.text.x=element_text(size=a.text.x),
#               axis.title.y=element_text(size=a.title.y), axis.text.y=element_text(size=a.text.y))
# p = p + theme(axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
#               axis.title.y = element_text(margin = unit(c(0, 2, 0, 0), "mm")),
#               plot.title = element_text(margin = unit(c(0, 2, 0, 0), "mm")))
# p
# 
# 
# ### This is data from a model that fits both location and scale.
# mewe1 = readRDS("mewe1_normal_ls_cauchy_M1000_N20_m10000.RDS")
# df_mewe1 = lapply(1:(length(n)), function(k) t(cbind(mewe1[[k]],rep(k,M))))
# df_mewe1 = matrix(unlist(df_mewe1), ncol = 3, byrow = TRUE)
# df_mewe1 = data.frame(df_mewe1)
# names(df_mewe1) = c("mu","sigma","n")
# 
# df_mewe1_scaled = lapply(1:(length(n)), function(k) t(cbind(sqrt(n[k])*mewe1[[k]],rep(k,M))))
# df_mewe1_scaled = matrix(unlist(df_mewe1_scaled), ncol = 3, byrow = TRUE)
# df_mewe1_scaled = data.frame(df_mewe1_scaled)
# names(df_mewe1_scaled) = c("mu","sigma","n")
# 
# #width of plots
# w = 10
# #height of plots
# h = 5
# 
# #limits for axes
# yliml = 0
# ylimu = 1.5
# xliml = 1
# xlimu = 3
# 
# #legend position
# leg.pos = "none"
# 
# #fontsizes
# title.size = 10
# a.title.x = 10
# a.text.x = 8
# a.title.y = 10
# a.text.y = 8
# 
# g <- ggplot(data = df_mewe1, aes(x = mu, y = sigma, colour = n, group = n)) #+ ylim(0, 1.5) + xlim(1,3)
# g <- g + geom_point(alpha = 0.5)
# g <- g + scale_colour_gradient2(midpoint = 4) + theme(legend.position = leg.pos, plot.title = element_text(size=title.size))
# g <- g + theme(axis.title.x=element_text(size=a.title.x), axis.text.x=element_text(size=a.text.x),
#                          axis.title.y=element_text(size=a.title.y), axis.text.y=element_text(size=a.text.y))
# g <- g + theme(axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
#                          axis.title.y = element_text(margin = unit(c(0, 2, 0, 0), "mm")),
#                          plot.title = element_text(margin = unit(c(0, 2, 0, 0), "mm")))
# g <- g + xlab(expression(mu)) + ylab(expression(sigma)) + ggtitle(paste("m = ",m,sep=""))
# g
# 
# 
# g <- ggplot(data = df_mewe1_scaled, aes(x = mu, group = n)) #+ ylim(0, 1.5) + xlim(1,3)
# g = g + geom_density(aes(y=..density..))
# g <- g + scale_colour_gradient2(midpoint = 4) + theme(legend.position = leg.pos, plot.title = element_text(size=title.size))
# g <- g + theme(axis.title.x=element_text(size=a.title.x), axis.text.x=element_text(size=a.text.x),
#                axis.title.y=element_text(size=a.title.y), axis.text.y=element_text(size=a.text.y))
# g <- g + theme(axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
#                axis.title.y = element_text(margin = unit(c(0, 2, 0, 0), "mm")),
#                plot.title = element_text(margin = unit(c(0, 2, 0, 0), "mm")))
# g <- g + xlab(expression(mu)) + ylab(expression(sigma)) + ggtitle(paste("m = ",m,sep=""))
# g
# 
# g <- ggplot(data = df_mewe1_scaled, aes(x = sigma, group = n)) #+ ylim(0, 1.5) + xlim(1,3)
# g = g + geom_density(aes(y=..density..))
# g



