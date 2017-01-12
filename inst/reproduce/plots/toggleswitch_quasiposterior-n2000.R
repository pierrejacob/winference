#+ presets, echo = FALSE, warning = FALSE, message = FALSE
library(winference)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(foreach)
library(doMC)
library(doRNG)
library(reshape2)
library(tidyr)
registerDoMC(cores = 6)
rm(list = ls())
setmytheme()
set.seed(11)
target <- get_toggleswitch()
# number of observations
nobservations <- 2000
load(file = "~/Dropbox/ABCD/Results/data/toggleswitchdata.RData")

observations <- obs[1:nobservations]
g <- qplot(x = observations, geom = "blank") + geom_histogram(aes(y = ..density..)) + xlab(expression(y))
# ggsave(filename = "~/Dropbox/ABCD/draft5/toggleswitch_data.pdf", plot = g, height = 5, width = 5)

filename <- paste0("~/Dropbox/ABCD/Results/toggleswitch/toggleswitchdata.n",
                   nobservations, ".wsmc_rhit.RData")
#
load(filename)
wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps <- max(wsmc.df$step)
target$parameter_names

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = alpha_1, group = step, colour = step))
g <- g + theme(legend.position = "none")
g <- g + geom_density(aes(y = ..density..)) + xlab(expression(alpha[1])) + geom_vline(xintercept = true_theta[1])
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
g
ggsave(filename = "~/Dropbox/ABCD/draft5/toggleswitch_marginal1_quasiposterior.n2000.pdf", plot = g, height = 5, width = 5)

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = alpha_2, group = step, colour = step))
g <- g + theme(legend.position = "none")
g <- g + geom_density(aes(y = ..density..)) + xlab(expression(alpha[2])) + geom_vline(xintercept = true_theta[2])
g
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
ggsave(filename = "~/Dropbox/ABCD/draft5/toggleswitch_marginal2_quasiposterior.n2000.pdf", plot = g, height = 5, width = 5)

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = beta_1, group = step, colour = step))
g <- g + theme(legend.position = "none")
g <- g + geom_density(aes(y = ..density..)) + xlab(expression(beta[1])) + geom_vline(xintercept = true_theta[3])
g
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
ggsave(filename = "~/Dropbox/ABCD/draft5/toggleswitch_marginal3_quasiposterior.n2000.pdf", plot = g, height = 5, width = 5)

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = beta_2, group = step, colour = step))
g <- g + theme(legend.position = "none")
g <- g + geom_density(aes(y = ..density..)) + xlab(expression(beta[2])) + geom_vline(xintercept = true_theta[4])
g
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
ggsave(filename = "~/Dropbox/ABCD/draft5/toggleswitch_marginal4_quasiposterior.n2000.pdf", plot = g, height = 5, width = 5)

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = mu, group = step, colour = step))
g <- g + theme(legend.position = "none")
g <- g + geom_density(aes(y = ..density..)) + xlab(expression(mu)) + geom_vline(xintercept = true_theta[5])
g
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
ggsave(filename = "~/Dropbox/ABCD/draft5/toggleswitch_marginal5_quasiposterior.n2000.pdf", plot = g, height = 5, width = 5)

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = sigma, group = step, colour = step))
g <- g + theme(legend.position = "none")
g <- g + geom_density(aes(y = ..density..)) + xlab(expression(sigma)) + geom_vline(xintercept = true_theta[6])
g
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
ggsave(filename = "~/Dropbox/ABCD/draft5/toggleswitch_marginal6_quasiposterior.n2000.pdf", plot = g, height = 5, width = 5)

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = gamma, group = step, colour = step))
g <- g + theme(legend.position = "none")
g <- g + geom_density(aes(y = ..density..)) + xlab(expression(gamma)) + geom_vline(xintercept = true_theta[7])
g
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
ggsave(filename = "~/Dropbox/ABCD/draft5/toggleswitch_marginal7_quasiposterior.n2000.pdf", plot = g, height = 5, width = 5)

#
# g <- ggplot(wsmc.df, aes(x = alpha_1, y = alpha_2, colour = step, group = step))
# g <- g + geom_point(alpha = 0.5)
# g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
# g <- g + xlab(expression(alpha[1])) + ylab(expression(alpha[2]))
# # g <- g + geom_point(data=NULL, aes(x = true_theta[1], y = true_theta[2], colour = NULL, group = NULL),
# # size = 5)
# g <- g + geom_hline(yintercept = true_theta[2]) + geom_vline(xintercept = true_theta[1])
# g
#
#
# g <- ggplot(wsmc.df, aes(x = beta_1, y = beta_2, colour = step, group = step))
# g <- g + geom_point(alpha = 0.5)
# g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
# g <- g + xlab(expression(beta[1])) + ylab(expression(beta[2]))
# g + geom_point(data=NULL, aes(x = true_theta[3], y = true_theta[4], colour = NULL, group = NULL),
#                size = 5)
#
# g <- ggplot(wsmc.df, aes(x = mu, y = sigma, colour = step, group = step))
# g <- g + geom_point(alpha = 0.5)
# g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
# g <- g + xlab(expression(mu)) + ylab(expression(sigma))
# g + geom_point(data=NULL, aes(x = true_theta[5], y = true_theta[6], colour = NULL, group = NULL),
#                size = 5)
#
# g <- ggplot(wsmc.df, aes(x = sigma, y = gamma, colour = step, group = step))
# g <- g + geom_point(alpha = 0.5)
# g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
# g <- g + xlab(expression(sigma)) + ylab(expression(gamma))
# g + geom_point(data=NULL, aes(x = true_theta[6], y = true_theta[7], colour = NULL, group = NULL),
#                size = 5)
#
# g <- ggplot(wsmc.df %>% filter(step > 20), aes(x = gamma, group = step)) + geom_density(aes(y = ..density..))
# g <- g + theme(legend.position = "none")
# g <- g + xlab(expression(gamma))
# g
