library(winference)
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
target <- list()
# generate random variables used to compute a synthetic dataset
target$generate_randomness <- function(nobservations){
  return(rnorm(nobservations))
}
# function to compute a dataset for each theta value, given fixed randomness
target$robservation <- function(theta, randomness){
  observations <- theta[1] + randomness * theta[2]
  return(observations)
}

gen_obs_data = function(n){rcauchy(n)}

metricL1 <- function(xvec,yvec)  mean(abs(xvec - yvec))

#Pick m to be larger than or equal to max(n) and a multiple of each entry in n
M = 1000
N = 20
m = 10^4
n = c(50,100,250,500,1000,5000,10000)

filename <- paste0(prefix,"cauchy_optim_m",m,"_k",N,"_M",M,".RData")

t = proc.time()
mewe_cauchy = foreach(rep = 1:M, .combine = rbind) %dorng% {

  #Allocate space for output
  mewe1_store = matrix(0,length(n),target$thetadim)
  mewe1_runtimes = rep(0,length(n))
  mewe1_evals = rep(0,length(n))

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

    #Initial point for optimization
    init = c(runif(1,-5,5),runif(1,0.2,3))

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
  output = cbind(mewe1_store,mewe1_runtimes,mewe1_evals,n)
  output
}
t = proc.time() - t
df.mewe = data.frame(mewe_cauchy)
names(df.mewe) = c("mu","sigma","runtime","fn.evals","n")
save(df.mewe,file = filename)

#calculate limiting MWE
mm = 10^8
nn = 10^8
N = 1

filename2 = paste0(prefix,"mwe_limit_m",mm,"_n",nn,".Rdata")

obs = gen_obs_data(nn)
sort_obs_mult = sort(obs)
sort_randomness = t(sapply(1:N, function(x) sort(target$generate_randomness(mm))))

#Define optimization objectives
mewe1_obj = function(theta){
  wass_dists = apply(sort_randomness, MARGIN = 1 , function(x) metricL1(sort_obs_mult,(theta[2]*x+theta[1])))
  out = mean(wass_dists)
  return(out)
}

init = c(0,2.2)

#Optimization
t_mewe1 = proc.time()
mewe1 = optim(init,mewe1_obj)
t_mewe1 = proc.time() - t_mewe1

#Save the results
mewe1_store = mewe1$par
mewe1_runtimes = t_mewe1[3]
mewe1_evals = mewe1$count

save(mewe1,t_mewe1,file = filename2)


# ################# Plots #################
load(filename)
load(filename2)

df.mewe$mu.scaled = df.mewe$mu*sqrt(df.mewe$n)
df.mewe$sigma.scaled = (df.mewe$sigma-mewe1$par[2])*sqrt(df.mewe$n)

df.mewe$gp = rep(1:7,M)
gp.order = order(df.mewe$n)
df.mewe = df.mewe[gp.order,]


g <- ggplot(data = df.mewe, aes(x = mu, y = sigma, colour = gp, group = gp)) #+ ylim(0, 1.5) + xlim(1,3)
g <- g + geom_point(alpha = 0.5)
g <- g +scale_colour_gradient2(midpoint = 4)
g <- g + xlab(expression(gamma)) + ylab(expression(sigma)) + theme(legend.position = "none")
g <- g + geom_vline(xintercept=mewe1$par[1]) + geom_hline(yintercept=mewe1$par[2])
g
ggsave(filename = paste0(prefix, "cauchy_mu_vs_sigma.png"), plot = g, width = 5, height = 5, dpi = 300)


g <- ggplot(data = df.mewe, aes(x = mu.scaled, group = gp, color =gp)) #+ ylim(0, 1.5) + xlim(1,3)
g = g + geom_density(aes(y=..density..), size = 1) + xlab(expression(gamma))
g = g + scale_colour_gradient2(midpoint = 4)
g = g + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "cauchy_rescaled_mu.pdf"), plot = g, width = 5, height = 5)

g <- ggplot(data = df.mewe, aes(x = sigma.scaled, group = gp, color =gp)) #+ ylim(0, 1.5) + xlim(1,3)
g = g + geom_density(aes(y=..density..), size = 1) + xlab(expression(sigma))
g = g + scale_colour_gradient2(midpoint = 4)
g = g + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "cauchy_rescaled_sigma.pdf"), plot = g, width = 5, height = 5)



