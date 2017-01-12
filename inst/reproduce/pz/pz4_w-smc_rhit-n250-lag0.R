#+ presets, echo = FALSE, warning = FALSE, message = FALSE
library(winference)
library(doMC)
library(doRNG)
library(dplyr)
library(ggthemes)

registerDoMC(cores = 6)

rm(list = ls())
set.seed(11)
setmytheme()
#
# model
target <- get_pz_4param()
#
# number of observations
nobservations <- 250
load(file = "~/Dropbox/ABCD/Results/data/pzdata.RData")
obs <- obs[1:nobservations]
# plot(obs, type = "l")
# now using the embeddings, and Hilbert sort
lagvalue <- 0


obs_sorted <- sort(obs)
# function to compute distance between observed data and data generated given theta
compute_d <- function(theta, metric = metricL2){
  fake_rand <- target$generate_randomness(nobservations)
  fake_obs <- target$robservation(nobservations, theta, target$parameters, fake_rand)
  fake_obs_sorted <- sort(fake_obs)
  return(metric(obs_sorted, fake_obs_sorted))
}

proposal <- mixture_proposal()

param_algo <- list(nthetas = 1024, nmoves = 1, proposal = proposal,
                   nsteps = 25, minimum_diversity = 0.5, R = 2, maxtrials = 1000)


filename <- paste0("~/Dropbox/ABCD/Results/pz/pzdata.n",
                   nobservations, ".L", lagvalue, ".wsmc_rhit-hilbert.RData")
# results <- wsmc_rhit(compute_d, target, param_algo, savefile = filename)
# wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
# nsteps <- max(wsmc.df$step)
# save(wsmc.df, results, nsteps, file = filename)
load(filename)
wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps <- max(wsmc.df$step)

target$parameter_names

# g <- ggplot(wsmc.df, aes(x = omega, group = step)) + geom_density(aes(y = ..density..), colour = "darkgrey")
# g <- g +  theme(legend.position = "none")
# g <- g + xlab(expression(omega))
# g
#
# g <- ggplot(wsmc.df, aes(x = phi, group = step)) + geom_density(aes(y = ..density..), colour = "darkgrey")
# g <- g +  theme(legend.position = "none")
# g <- g + xlab(expression(phi))
# g

g <- ggplot(wsmc.df, aes(x = mu_alpha, y = sigma_alpha, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(mu[alpha])) + ylab(expression(sigma[alpha]))
g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
g

g <- ggplot(wsmc.df, aes(x = c, y = e, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(c)) + ylab(expression(e))
g <- g + geom_vline(xintercept = true_theta[3]) + geom_hline(yintercept = true_theta[4])
g

g <- qplot(x = 1:length(results$threshold_history), y = results$threshold_history, geom = "line") +
  scale_y_log10()
g <- g + xlab("step") + ylab("threshold")
g

dist.df <- foreach(irep = 1:nrow(results$thetas_history[[nsteps]]), .combine = rbind) %dorng%{
  c(compute_d(results$thetas_history[[nsteps]][irep,]), compute_d(true_theta))
}

dist.df <- melt(data.frame(dist.df))
dist.df %>% head

ggplot(dist.df, aes(x = value, group = variable, fill = variable, colour = variable)) + geom_density(aes(y = ..density..), alpha = 0.5)


