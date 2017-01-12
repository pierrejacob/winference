#+ presets, echo = FALSE, warning = FALSE, message = FALSE
library(winference)
library(doMC)
library(doRNG)
library(dplyr)
library(ggthemes)

registerDoMC(cores = 6)

rm(list = ls())
set.seed(13)
setmytheme()
#
# model
target <- get_cosine()
#
# number of observations
nobservations <- 500
load(file = "~/Dropbox/ABCD/Results/data/cosinedata.RData")
obs <- obs[1:nobservations]
lagvalue <- 1

# plot(obs, type = "l")
reconstitute <- function(observations, theta){
  nobservations <- length(observations)
  eps <- (observations - exp(theta[4]) * cos(2*pi*theta[1]*(1:nobservations) + theta[2])) / exp(theta[3])
  return(eps)
}

#
regularizer <- 0.05
compute_d <- function(theta, transportiterations = 100){
  eps <- reconstitute(obs, theta)
  lag_obs <- create_lagmatrix(matrix(eps, nrow = 1), lagvalue)
  lag_obs <- lag_obs[,(lagvalue+1):ncol(lag_obs)]
  # order_obs <- hilbert_order(lag_obs)
  # orderded_obs <- lag_obs[,order_obs]
  fake_eps <- rnorm(length(obs))
  lag_fake_obs <- create_lagmatrix(matrix(fake_eps, nrow = 1), lagvalue)
  lag_fake_obs <- lag_fake_obs[,(lagvalue+1):ncol(lag_fake_obs)]
  # order_fake_obs <- hilbert_order(lag_fake_obs)
  # orderded_fake_obs <- lag_fake_obs[,order_fake_obs]
  C <- cost_matrix_L2(lag_obs, lag_fake_obs)
  epsilon <- regularizer * median(C)
  equalw <- rep(1/(nobservations-lagvalue), (nobservations-lagvalue))
  wass <- wasserstein(equalw, equalw, C, epsilon, transportiterations)
  return(as.numeric(wass$distances))
}

proposal <- mixture_proposal()

param_algo <- list(nthetas = 1024, nmoves = 1, proposal = proposal,
                   nsteps = 40, minimum_diversity = 0.5, R = 2, maxtrials = 1000)

filename <- paste0("~/Dropbox/ABCD/Results/cosine/cosinedata.n",
                   nobservations, "residualsdelay.wsmc_rhit-wasserstein.RData")
# results <- wsmc_rhit(compute_d, target, param_algo, savefile = filename)
# wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
# nsteps <- max(wsmc.df$step)
# save(wsmc.df, results, nsteps, file = filename)
target$parameter_names
load(filename)
wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps <- max(wsmc.df$step)
#
g <- ggplot(wsmc.df, aes(x = omega, y = logsigma, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(omega)) + ylab(expression(log(sigma)))
g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[3])
g

g <- ggplot(wsmc.df, aes(x = phi, y = logA, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(phi)) + ylab(expression(log(A)))
g <- g + geom_vline(xintercept = true_theta[2]) + geom_hline(yintercept = true_theta[4])
g

g <- ggplot(wsmc.df, aes(x = phi, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(phi))
g <- g + geom_vline(xintercept = true_theta[2])
g

g <- qplot(x = 1:length(results$threshold_history), y = results$threshold_history, geom = "line") +
  scale_y_log10()
g <- g + xlab("step") + ylab("threshold")
g



# example of a theta very different from true_theta but
# yielding nearly same distance
# theta <-  c(0.04241519, 3.00471029, 0.57630632, 0.11245227)
# eps <- reconstitute(obs, theta)
# plot(eps, type = "l")
# hist(eps, prob = T, nclass = 30)
# curve(dnorm(x), add =T)
#

###
dist.df <- foreach(irep = 1:nrow(results$thetas_history[[nsteps]]), .combine = rbind) %dorng%{
  c(compute_d(results$thetas_history[[nsteps]][irep,]), compute_d(true_theta))
}

dist.df <- melt(data.frame(dist.df))
dist.df %>% head

ggplot(dist.df, aes(x = value, group = variable, fill = variable, colour = variable)) + geom_density(aes(y = ..density..), alpha = 0.5)

g <- ggplot(wsmc.df, aes(x = omega, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(omega))
g <- g + geom_vline(xintercept = true_theta[1])
g
# ggsave(filename = "~/Dropbox/ABCD/draft4/cosine_withlag3_omega.pdf", plot = g, width = 5, height = 5 )

g <- ggplot(wsmc.df, aes(x = phi, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(phi))
g <- g + geom_vline(xintercept = true_theta[2])
g
# ggsave(filename = "~/Dropbox/ABCD/draft4/cosine_withlag3_phi.pdf", plot = g, width = 5, height = 5 )

g <- ggplot(wsmc.df, aes(x = logsigma, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(log(sigma)))
g <- g + geom_vline(xintercept = true_theta[3])
g
# ggsave(filename = "~/Dropbox/ABCD/draft4/cosine_withlag3_logsigma.pdf", plot = g, width = 5, height = 5 )

g <- ggplot(wsmc.df, aes(x = logA, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(log(A)))
g <- g + geom_vline(xintercept = true_theta[4])
g
# ggsave(filename = "~/Dropbox/ABCD/draft4/cosine_withlag3_logA.pdf", plot = g, width = 5, height = 5 )

