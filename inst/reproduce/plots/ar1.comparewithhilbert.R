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
target <- get_autoregressive()
#
# number of observations
nobservations <- 1000
load(file = "~/Dropbox/ABCD/Results/data/ar1data.RData")

### load the results obtained with delay reconstruction, with a lag of 1
### with hilbert sort
filename <- paste0("~/Dropbox/ABCD/Results/autoregressive/ar1data.n1000.wsmc-hilbert.RData")
load(filename)
wsmc.hilbert.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps.hilbert <- max(wsmc.hilbert.df$step)

### load the results obtained with delay reconstruction, with a lag of 1
### with Wasserstein distance
filename <- paste0("~/Dropbox/ABCD/Results/autoregressive/ar1data.n1000.wsmc_rhit.Lag1.wassersteinL2.RData")
load(filename)
wsmc.wasserstein.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps.wasserstein <- max(wsmc.wasserstein.df$step)


### now plot figures
# g <- ggplot(wsmc.wasserstein.df, aes(x = rho, y = logsigma, colour = step, group = step))
# g <- g + geom_point(alpha = 0.5)
# g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
# g <- g + xlab(expression(phi)) + ylab(expression(log(sigma))) +  xlim(-1,1) + ylim(-4,4)
# g
# ggsave(filename = "~/Dropbox/ABCD/draft4/ar1_identification_withoutlag.pdf", plot = g, width = 7, height = 5 )
# g <- ggplot(wsmc.wasserstein.df, aes(x = rho, y = logsigma, colour = step, group = step))
# g <- g + geom_point(alpha = 0.5)
# g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
# g <- g + xlab(expression(phi)) + ylab(expression(log(sigma))) +  xlim(-1,1) + ylim(-4,4)
# g
# ggsave(filename = "~/Dropbox/ABCD/draft4/ar1_identification_withlag1.pdf", plot = g, width = 7, height = 5 )

g <- ggplot(wsmc.wasserstein.df %>% filter(step > 2, step < 11), aes(x = rho, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(phi))
g <- g + geom_density(data=wsmc.hilbert.df %>% filter(step > 2, step < 11), colour = "red", linetype = 2)
g
ggsave(filename = "~/Dropbox/ABCD/draft4/ar1_withhilbert_logsigma.pdf", plot = g, width = 7, height = 5 )

g <- ggplot(wsmc.wasserstein.df %>% filter(step > 2, step < 11), aes(x = logsigma, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(phi))
g <- g + geom_density(data=wsmc.hilbert.df %>% filter(step > 2, step < 11), colour = "red", linetype = 2)
g

ggsave(filename = "~/Dropbox/ABCD/draft4/ar1_withhilbert_phi.pdf", plot = g, width = 7, height = 5 )

