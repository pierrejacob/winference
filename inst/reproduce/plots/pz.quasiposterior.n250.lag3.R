#+ presets, echo = FALSE, warning = FALSE, message = FALSE
library(winference)
library(doMC)
library(doRNG)
library(dplyr)
library(ggthemes)

registerDoMC(cores = 4)

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
lagvalue <- 3

filename <- paste0("~/Dropbox/ABCD/Results/pz/pzdata.n250.L3.wsmc_rhit-wassersteinL2.RData")
load(filename)
wsmc.wasserstein.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps.wasserstein <- max(wsmc.wasserstein.df$step)


filename <- paste0("~/Dropbox/ABCD/Results/pz/pzdata.n250.L3.wsmc_rhit-hilbert.RData")
load(filename)
wsmc.hilbert.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps.hilbert <- max(wsmc.hilbert.df$step)

# joint
g <- ggplot(wsmc.wasserstein.df %>% filter(step < 16), aes(x = mu_alpha, y = sigma_alpha, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(16/2)) + theme(legend.position = "none")
g <- g + xlab(expression(mu[alpha])) + ylab(expression(sigma[alpha]))
g <- g + geom_hline(yintercept = true_theta[2]) + geom_vline(xintercept = true_theta[1])
g

ggsave(filename = "~/Dropbox/ABCD/draft4/pz_musigma_quasiposterior.n250.lag3.wasserstein.png", plot = g, height = 5, width = 7, dpi = 150)

g <- ggplot(wsmc.hilbert.df %>% filter(step < 16), aes(x = mu_alpha, y = sigma_alpha, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(16/2)) + theme(legend.position = "none")
g <- g + xlab(expression(mu[alpha])) + ylab(expression(sigma[alpha]))
g <- g + geom_hline(yintercept = true_theta[2]) + geom_vline(xintercept = true_theta[1])
g

g <- ggplot(wsmc.wasserstein.df %>% filter(step < 16), aes(x = c, y = e, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(16/2)) + theme(legend.position = "none")
g <- g + xlab(expression(c)) + ylab(expression(e))
g <- g + geom_hline(yintercept = true_theta[4]) + geom_vline(xintercept = true_theta[3])
g

ggsave(filename = "~/Dropbox/ABCD/draft4/pz_ce_quasiposterior.n250.lag3.wasserstein.png", plot = g, height = 5, width = 7, dpi = 150)

g <- ggplot(wsmc.hilbert.df %>% filter(step < 16), aes(x = c, y = e, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(16/2)) + theme(legend.position = "none")
g <- g + xlab(expression(c)) + ylab(expression(e))
g <- g + geom_hline(yintercept = true_theta[4]) + geom_vline(xintercept = true_theta[3])
g


# marginals
g <- ggplot(wsmc.wasserstein.df %>% filter(step == 15), aes(x = mu_alpha, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(mu[alpha]))
g <- g + geom_density(data=wsmc.hilbert.df %>% filter(step == 15), linetype = 2, colour = "red")
g <- g + geom_vline(xintercept = true_theta[1])
g

g <- ggplot(wsmc.wasserstein.df %>% filter(step == 15), aes(x = sigma_alpha, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(sigma[alpha]))
g <- g + geom_density(data=wsmc.hilbert.df %>% filter(step  == 15), linetype = 2, colour = "red")
g <- g + geom_vline(xintercept = true_theta[2])
g

g <- ggplot(wsmc.wasserstein.df %>% filter(step == 15), aes(x = c, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(c))
g <- g + geom_density(data=wsmc.hilbert.df %>% filter(step == 15), linetype = 2, colour = "red")
g <- g + geom_vline(xintercept = true_theta[3])
g

g <- ggplot(wsmc.wasserstein.df %>% filter(step == 15), aes(x = e, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(e))
g <- g + geom_density(data=wsmc.hilbert.df %>% filter(step == 15), linetype = 2, colour = "red")
g <- g + geom_vline(xintercept = true_theta[4])
g
