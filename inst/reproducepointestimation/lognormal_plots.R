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

prefix = "plots/"

L <- 10

target <- list()

target$generate_randomness <- function(nobservations){
  return(rnorm(nobservations*L))
}

target$robservation <- function(theta, randomness){
  normals_ <- theta[1] + theta[2] * randomness
  lognormals_ <- exp(normals_)
  return(rowSums(matrix(lognormals_, ncol = L)))
}

metricL1 = function(xvec, yvec)  mean(abs(xvec - yvec))

true_theta <- c(0,1)
target$thetadim = 2 #Inference on all parameters


#Pick m to be larger than or equal to max(n) and a multiple of each entry in n
M = 1000
N = 20
m = 10^4
n = c(50,100,250,500,1000,5000,10000)

filename <- paste0(prefix,"lognormal_pierre_optim_m",m,"_k",N,"_M",M,".RData")


# t = proc.time()
mewe_lognormal = foreach(rep = 1:M, .combine = rbind) %dorng% {
  # rep <- 1
  #Allocate space for output
  mewe_store = matrix(0,length(n),target$thetadim)
  mewe_runtimes = rep(0,length(n))
  mewe_evals = rep(0,length(n))

  #generate all observations and sets of randomness to be used
  obs_rand = target$generate_randomness(max(n))
  obs_all = target$robservation(true_theta,obs_rand)

  #Generate the synthetic randomness, sort.
  randomness = t(sapply(1:N, function(x) target$generate_randomness(m)))

  for(i in 1:(length(n))){
    print(i)
    #Subset observations and sort
    obs = obs_all[1:n[i]]
    sort_obs = sort(obs)
    sort_obs_mult = rep(sort_obs, each = m/n[i])

    #Define the objective to be minimized to find the MEWE
    obj1 = function(theta){
      if(theta[2] < 0){
        out = 10^6
      } else{
        wass_dists = apply(randomness, 1, function(x) metricL1(sort_obs_mult, sort(target$robservation(theta, x))))
        out = mean(wass_dists)
      }
      return(out)
    }

    #Optimization
    t_mewe = proc.time()
    mewe = optim(true_theta,obj1)
    t_mewe = proc.time() - t_mewe

    #Save the results
    mewe_store[i,] = mewe$par
    mewe_runtimes[i] = t_mewe[3]
    mewe_evals[i] = mewe$counts[1]
  }

  #Output
  output = cbind(mewe_store,mewe_runtimes,mewe_evals,n,(1:length(n)))
  output
}


  # t = proc.time() - t
df.mewe = data.frame(mewe_lognormal)
names(df.mewe) = c("mu","sigma","runtime","fn.evals","n","gp")
save(df.mewe,t,file = filename)

# ################# Plots #################
load(filename)

df.mewe$mu.scaled = (df.mewe$mu - true_theta[1])*sqrt(df.mewe$n)
df.mewe$sigma.scaled = (df.mewe$sigma - true_theta[2])*sqrt(df.mewe$n)

gp.order = order(df.mewe$n)
df.mewe = df.mewe[gp.order,]


g <- ggplot(data = df.mewe, aes(x = mu, y = sigma, colour = gp, group = gp)) #+ ylim(0, 1.5) + xlim(1,3)
g <- g + geom_point(alpha = 0.5)
g <- g +scale_colour_gradient2(midpoint = 4)
g <- g + xlab(expression(mu)) + ylab(expression(sigma)) + theme(legend.position = "none")
g <- g + geom_vline(xintercept=true_theta[1]) + geom_hline(yintercept=true_theta[2])
g
ggsave(filename = paste0(prefix, "lognormal_mu_vs_sigma.png"), plot = g, width = 7, height = 5, dpi = 300)


g <- ggplot(data = df.mewe, aes(x = mu.scaled, group = gp, color = gp)) #+ ylim(0, 1.5) + xlim(1,3)
g = g + scale_colour_gradient2(midpoint = 4)
g = g + geom_density(aes(y=..density..), size = 1) + xlab(expression(mu)) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "lognormal_rescaled_mu.pdf"), plot = g, width = 7, height = 5)


