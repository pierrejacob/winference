library(doParallel)
library(ggplot2)
registerDoParallel(cores = detectCores())
rm(list = ls())
theme_set(theme_bw())
theme_update(axis.text.x = element_text(size = 20),
             axis.text.y = element_text(size = 20),
             axis.title.x = element_text(size = 25, margin=margin(20,0,0,0)),
             axis.title.y = element_text(size = 25, angle = 90, margin = margin(0,20,0,0)),
             legend.text = element_text(size = 20),
             legend.title = element_text(size = 20),
             title = element_text(size = 30),
             strip.text = element_text(size = 25),
             strip.background = element_rect(fill="white"),
             panel.spacing = unit(2, "lines"),
             legend.position = "bottom")
set.seed(11)

prefix = ""

source("gamma_normal_functions.R")

gen_obs_data = function(n){rgamma(n,10,5)}


#Pick m to be larger than or equal to max(n) and a multiple of each entry in n
M = 1000
N = 20
m = 10^4
n = c(50,100,250,500,1000,5000,10000)


######### Fits both location and scale.

filename <- paste0(prefix,"gamma_normal_m",m,"_k",N,"_M",M,".RData")

t = proc.time()
mewe_gamma = foreach(rep = 1:M, .combine = rbind) %dorng% {

  #Allocate space for output
  mewe1_store = matrix(0,length(n),target$thetadim)
  mewe1_runtimes = rep(0,length(n))
  mewe1_evals = rep(0,length(n))

  mle_store = matrix(0,length(n),target$thetadim)

  #generate all observations and sets of randomness to be used
  obs_all = gen_obs_data(max(n))
  sort_randomness = t(sapply(1:N, function(x) sort(target$generate_randomness(m))))

  for(i in 1:(length(n)) ){
    #Subset observations and sort
    obs = obs_all[1:n[i]]
    sort_obs = sort(obs)
    sort_obs_mult = rep(sort_obs, each = m/n[i])

    #Compute the MLE
    mu.mle = mean(obs)
    sigma.mle = sqrt((n[i]-1)/n[i])*sd(obs)
    mle_store[i,] = c(mu.mle,sigma.mle)

    #Define optimization objectives
    mewe1_obj = function(theta){
      wass_dists = apply(sort_randomness, MARGIN = 1 , function(x) metricL1(sort_obs_mult,(theta[2]*x+theta[1])))
      out = mean(wass_dists)
      return(out)
    }

    #Initial point for optimization
    init = c(runif(1,-3,5),runif(1,0.1,2))

    #Optimization
    t_mewe1 = proc.time()
    mewe1 = optim(init,mewe1_obj)
    t_mewe1 = proc.time() - t_mewe1

    #Save the results
    mewe1_store[i,] = mewe1$par
    mewe1_runtimes[i] = t_mewe1[3]
    mewe1_evals[i] = mewe1$count
  }
  #Output
  output = cbind(mewe1_store,mewe1_runtimes,mewe1_evals,mle_store,n,(1:length(n)))
  output
}
t = proc.time() - t
df.mewe = data.frame(mewe_gamma)
names(df.mewe) = c("mu","sigma","runtime","fn.evals","mu.mle","sigma.mle","n","gp")
save(df.mewe,file = filename)


load(filename)

gp.order = order(df.mewe$n)
df.mewe = df.mewe[gp.order,]

g <- ggplot(data = df.mewe, aes(x = mu, y = sigma, colour = gp, group = gp)) + ylim(0.35, 0.9) + xlim(1.6,2.4)
g <- g + geom_point(alpha = 0.5)
g <- g +scale_colour_gradient2(midpoint = 4)
g <- g + xlab(expression(gamma)) + ylab(expression(sigma)) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "gamma_mewe_mu_vs_sigma.png"), plot = g, width = 5, height = 5, dpi = 300)

g <- ggplot(data = df.mewe, aes(x = mu.mle, y = sigma.mle, colour = gp, group = gp)) + ylim(0.35, 0.9) + xlim(1.6,2.4)
g <- g + geom_point(alpha = 0.5)
g <- g +scale_colour_gradient2(midpoint = 4)
g <- g + xlab(expression(gamma)) + ylab(expression(sigma)) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "gamma_mle_mu_vs_sigma.png"), plot = g, width = 5, height = 5, dpi = 300)


g <- ggplot(data = df.mewe %>% filter(n == 10000), aes(x = mu)) + geom_density(aes(y=..density..))
g = g + geom_density(aes(x = mu.mle, y=..density..),linetype="dashed") + xlab(expression(gamma))
g
ggsave(filename = paste0(prefix, "gamma_mu_mle_vs_mewe_n10000.pdf"), plot = g, width = 5, height = 5)

g <- ggplot(data = df.mewe %>% filter(n == 10000), aes(x = sigma)) + geom_density(aes(y=..density..)) #+ ylim(0, 1.5) + xlim(1,3)
g = g + geom_density(aes(x = sigma.mle, y=..density..),linetype="dashed") + xlab(expression(sigma))
g
ggsave(filename = paste0(prefix, "gamma_sigma_mle_vs_mewe_n10000.pdf"), plot = g, width = 5, height = 5)
