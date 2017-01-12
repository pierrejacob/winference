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
target <- get_levydriven()
load(file = "~/Dropbox/ABCD/Results/data/levydrivendata.RData")
load(file = "~/Dropbox/ABCD/Results/levydriven/levydriven.n10000.L1.wsmc_rhit-hilbert.RData")
wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps <- max(wsmc.df$step)
#
results$threshold_history
target$parameter_names

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = mu, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(mu)) + geom_vline(xintercept = true_theta[1])
g

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = beta, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(beta)) + geom_vline(xintercept = true_theta[2])
g

g <- ggplot(wsmc.df %>% filter(step > 10), aes(x = xi, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(xi)) + geom_vline(xintercept = true_theta[3])
g

g <- ggplot(wsmc.df %>% filter(step > 30), aes(x = omega2, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(omega^2)) + geom_vline(xintercept = true_theta[4])
g

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = lambda, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + xlab(expression(lambda)) + geom_vline(xintercept = true_theta[5])
# g + scale_colour_gradient2(midpoint = floor(20), mid = "white")
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
g
ggsave(filename = "~/Dropbox/ABCD/draft5/levydriven.n10000.L1.hilbert.lambda.pdf", plot = g, height = 5, width = 7)

g <- qplot(x = 1:length(results$threshold_history), y = results$threshold_history, geom = "line") +
  scale_y_log10()
g <- g + xlab("step") + ylab("threshold")
g

##

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = mu, y = beta, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(20)) + theme(legend.position = "none")
g <- g + xlab(expression(mu)) + ylab(expression(beta))
g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
g
ggsave(filename = "~/Dropbox/ABCD/draft5/levydriven.n10000.L1.hilbert.mubeta.png", plot = g, height = 5, width = 7, dpi = 150)

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = xi, y = omega2, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(20)) + theme(legend.position = "none")
g <- g + xlab(expression(xi)) + ylab(expression(omega^2))
g <- g + geom_vline(xintercept = true_theta[3]) + geom_hline(yintercept = true_theta[4])
g
ggsave(filename = "~/Dropbox/ABCD/draft5/levydriven.n10000.L1.hilbert.xiomega2.png", plot = g, height = 5, width = 7, dpi = 150)

g <- ggplot(wsmc.df, aes(x = omega2, y = lambda, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + ylab(expression(lambda)) + xlab(expression(omega^2))
g <- g + geom_vline(xintercept = true_theta[4]) + geom_hline(yintercept = true_theta[5])
g

# ggsave(filename = "~/Dropbox/ABCD/draft5/levydriven.n10000.L1.hilbert.png", plot = g, height = 5, width = 7, dpi = 150)
