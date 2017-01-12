library(winference)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(foreach)
library(doMC)
library(doRNG)
library(reshape2)
library(tidyr)
registerDoMC(cores = 4)
rm(list = ls())
setmytheme()
set.seed(11)
target <- get_gandk()
# number of observations
nobservations <- 250
load(file = "~/Dropbox/ABCD/Results/data/gandkdata.RData")
obs <- obs[1:nobservations]
obs_sorted <- sort(obs)

filename <- paste0("~/Dropbox/ABCD/Results/gandk/gandkdata.n",
                   nobservations, ".wsmc_rhit.RData")
#
load(filename)
wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps <- max(wsmc.df$step)

g <- qplot(x = 1:length(results$threshold_history), y = results$threshold_history, geom = "line") +
  scale_y_log10()
g <- g + xlab("step") + ylab("threshold")
g
ggsave(filename = "~/Dropbox/ABCD/draft5/gandk_thresholds.pdf", plot = g, height = 5, width = 7)


g <- ggplot(wsmc.df, aes(x = A, y = B, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(A)) + ylab(expression(B))
g <- g + geom_hline(yintercept = true_theta[2]) + geom_vline(xintercept = true_theta[1])
# g
ggsave(filename = "~/Dropbox/ABCD/draft5/gandk_AB_quasiposterior.n250.png", plot = g, height = 5, width = 7, dpi = 150)

g <- ggplot(wsmc.df, aes(x = g, y = k, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(g)) + ylab(expression(k))
g <- g + geom_hline(yintercept = true_theta[4]) + geom_vline(xintercept = true_theta[3])
# g
ggsave(filename = "~/Dropbox/ABCD/draft5/gandk_gk_quasiposterior.n250.png", plot = g, height = 5, width = 7, dpi = 150)


g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = g, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none")
g <- g + geom_vline(xintercept = true_theta[3])
g <- g + xlab(expression(g))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
g

ggsave(filename = "~/Dropbox/ABCD/draft5/gandk_gmarginal_quasiposterior.n250.pdf", plot = g, height = 5, width = 7)

