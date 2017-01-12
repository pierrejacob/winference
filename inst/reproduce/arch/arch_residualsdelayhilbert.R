library(winference)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(foreach)
library(doMC)
library(doRNG)
library(reshape2)
registerDoMC(cores = 6)
rm(list = ls())
setmytheme()

set.seed(11)

target <- get_arch()

nobservations <- 1000
load(file = "~/Dropbox/ABCD/Results/data/archdata.RData")
obs <- obs[1:nobservations]
plot(obs, type= "l")


reconstitute <- function(observations, theta){
  alpha0 <- theta[1]
  alpha1 <- theta[2]
  nobservations <- length(observations)
  variance <- rep(1, nobservations)
  residuals <- rep(0, nobservations)
  residuals[1] <- observations[1] / sqrt(variance[1])
  for (t in 2:nobservations){
    variance[t] <- alpha0 + alpha1 * observations[t-1]^2
    residuals[t] <-  observations[t] / sqrt(variance[t])
  }
  return(residuals)
}

eps <- reconstitute(obs, true_theta)
hist(eps, nclass = 100, prob = TRUE)
curve(dnorm(x), add = TRUE)


compute_d <- function(theta, metric = metricL2){
  eps <- reconstitute(obs, theta)
  eps_sorted <- sort(eps)
  fake_eps <- rnorm(length(obs))
  fake_eps_sorted <- sort(fake_eps)
  return(metric(eps_sorted, fake_eps_sorted))
}

proposal <- mixture_proposal()

param_algo <- list(nthetas = 1024, nmoves = 1, proposal = proposal,
                   nsteps = 20, minimum_diversity = 0.5, R = 2, maxtrials = 1000)

filename <- paste0("~/Dropbox/ABCD/Results/arch/arch.n",
                   nobservations, "residuals.wsmc_rhit.RData")
results <- wsmc_rhit(compute_d, target, param_algo, savefile = filename)
wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps <- max(wsmc.df$step)
save(wsmc.df, results, nsteps, file = filename)

load(filename)
wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps <- max(wsmc.df$step)
#
g <- ggplot(wsmc.df, aes(x = alpha_0, group = step)) + geom_density(aes(y = ..density..), colour = "darkgrey")
g <- g +  theme(legend.position = "none")
g1 <- g + xlab(expression(alpha[0])) + geom_vline(xintercept = true_theta[1])
g1

g <- ggplot(wsmc.df, aes(x = alpha_1, group = step)) + geom_density(aes(y = ..density..), colour = "darkgrey")
g <- g +  theme(legend.position = "none")
g2 <- g + xlab(expression(alpha[1])) + geom_vline(xintercept = true_theta[2])
g2


g <- ggplot(wsmc.df, aes(x = alpha_0, y = alpha_1, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(alpha[0])) + ylab(expression(alpha[1]))
g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
g + geom_rug(alpha = 0.2)

# library(gridExtra)
# empty <- ggplot()+geom_point(aes(1,1), colour="white")+
#   theme(axis.ticks=element_blank(),
#         panel.background=element_blank(),
#         axis.text.x=element_blank(), axis.text.y=element_blank(),
#         axis.title.x=element_blank(), axis.title.y=element_blank())
# grid.arrange(g1, empty, g, g2+coord_flip(), ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

# g + geom_rug(alpha = 0.2)
g + geom_rug(data = wsmc.df %>% filter(step == nsteps), alpha = 0.2)
true_theta
theta <- results$thetas_history[[nsteps]][150,]
eps <- reconstitute(obs, theta)
plot(eps, type= "l")
hist(eps, prob = T)
curve(dnorm(x), add = T)

dist.df <- foreach(irep = 1:nrow(results$thetas_history[[nsteps]]), .combine = rbind) %dorng%{
  c(compute_d(results$thetas_history[[nsteps]][irep,]), compute_d(true_theta))
}

dist.df <- melt(data.frame(dist.df))
dist.df %>% head

ggplot(dist.df, aes(x = value, group = variable, fill = variable, colour = variable)) + geom_density(aes(y = ..density..), alpha = 0.5)
