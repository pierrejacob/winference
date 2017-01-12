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
target <- get_cosine()
#
# number of observations
nobservations <- 500
load(file = "~/Dropbox/ABCD/Results/data/cosinedata.RData")
obs <- obs[1:nobservations]

filename <- paste0("~/Dropbox/ABCD/Results/cosine/cosinedata.n",
                   nobservations, "residualsdelay.wsmc_rhit-wasserstein.RData")

load(filename)
wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps <- max(wsmc.df$step)

g <- ggplot(wsmc.df %>% filter(step >= 0), aes(x = omega, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(omega))
g <- g + geom_vline(xintercept = true_theta[1])
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
g
ggsave(filename = "~/Dropbox/ABCD/draft5/cosine_residanddelay_omega.pdf", plot = g, width = 5, height = 5 )

g <- ggplot(wsmc.df %>% filter(step >= 0), aes(x = phi, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(phi))
g <- g + geom_vline(xintercept = true_theta[2])
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
g
ggsave(filename = "~/Dropbox/ABCD/draft5/cosine_residanddelay_phi.pdf", plot = g, width = 5, height = 5 )

g <- ggplot(wsmc.df %>% filter(step >= 0), aes(x = logsigma, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(log(sigma)))
g <- g + geom_vline(xintercept = true_theta[3])
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
g
ggsave(filename = "~/Dropbox/ABCD/draft5/cosine_residanddelay_logsigma.pdf", plot = g, width = 5, height = 5 )

g <- ggplot(wsmc.df %>% filter(step >= 0), aes(x = logA, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(log(A)))
g <- g + geom_vline(xintercept = true_theta[4])
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
g
ggsave(filename = "~/Dropbox/ABCD/draft5/cosine_residanddelay_logA.pdf", plot = g, width = 5, height = 5 )

#
#
# g <- ggplot(wsmc.df, aes(x = omega, y = logsigma, colour = step, group = step))
# g <- g + geom_point(alpha = 0.5)
# g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
# g <- g + xlab(expression(omega)) + ylab(expression(log(sigma)))
# g
#
# g <- ggplot(wsmc.df, aes(x = phi, y = logA, colour = step, group = step))
# g <- g + geom_point(alpha = 0.5)
# g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
# g <- g + xlab(expression(phi)) + ylab(expression(log(A)))
# g
#
