library(winference)
registerDoParallel(cores = 6)
rm(list = ls())
setmytheme()
set.seed(11)
target <- get_toggleswitch()
# number of observations
nobservations <- 2000
fig.height <- 5
fig.width <- 5

prefix <- ""
load(file = paste0(prefix, "toggleswitchdata.RData"))
obs <- obs[1:nobservations]
obs_sorted <- sort(obs)

g <- qplot(x = obs, geom = "blank") + geom_histogram(aes(y = ..density..))
g <- g + xlab("y")
g
ggsave(filename = paste0(prefix, "toggleswitch_data.pdf"), plot = g, width = fig.width, height = fig.height)

filename <- paste0(prefix, "toggleswitchwsmc.n", nobservations, ".RData")
load(filename)

plot_threshold_time(results)
plot_ncomputed(results)

nsteps <- length(results$thetas_history)

wsmc.df <- wsmc_to_dataframe(results)

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = alpha_1, colour = step, group = step))
g <- g + geom_density(aes(y = ..density..))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue") + theme(legend.position = "none")
g <- g + xlab(expression(alpha[1]))
g <- g + geom_vline(xintercept = true_theta[1])
g
ggsave(filename = paste0(prefix, "toggleswitch_marginal1.pdf"), plot = g, width = fig.width, height = fig.height)

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = alpha_2, colour = step, group = step))
g <- g + geom_density(aes(y = ..density..))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue") + theme(legend.position = "none")
g <- g + xlab(expression(alpha[2]))
g <- g + geom_vline(xintercept = true_theta[2])
g
ggsave(filename = paste0(prefix, "toggleswitch_marginal2.pdf"), plot = g, width = fig.width, height = fig.height)


g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = beta_1, colour = step, group = step))
g <- g + geom_density(aes(y = ..density..))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue") + theme(legend.position = "none")
g <- g + xlab(expression(beta[1]))
g <- g + geom_vline(xintercept = true_theta[3])
g
ggsave(filename = paste0(prefix, "toggleswitch_marginal3.pdf"), plot = g, width = fig.width, height = fig.height)

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = beta_2, colour = step, group = step))
g <- g + geom_density(aes(y = ..density..))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue") + theme(legend.position = "none")
g <- g + xlab(expression(beta[2]))
g <- g + geom_vline(xintercept = true_theta[4])
g
ggsave(filename = paste0(prefix, "toggleswitch_marginal4.pdf"), plot = g, width = fig.width, height = fig.height)

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = mu, colour = step, group = step))
g <- g + geom_density(aes(y = ..density..))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue") + theme(legend.position = "none")
g <- g + xlab(expression(mu))
g <- g + geom_vline(xintercept = true_theta[5])
g
ggsave(filename = paste0(prefix, "toggleswitch_marginal5.pdf"), plot = g, width = fig.width, height = fig.height)

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = sigma, colour = step, group = step))
g <- g + geom_density(aes(y = ..density..))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue") + theme(legend.position = "none")
g <- g + xlab(expression(sigma))
g <- g + geom_vline(xintercept = true_theta[6])
g
ggsave(filename = paste0(prefix, "toggleswitch_marginal6.pdf"), plot = g, width = fig.width, height = fig.height)

g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = gamma, colour = step, group = step))
g <- g + geom_density(aes(y = ..density..))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue") + theme(legend.position = "none")
g <- g + xlab(expression(gamma))
g <- g + geom_vline(xintercept = true_theta[7])
g
ggsave(filename = paste0(prefix, "toggleswitch_marginal7.pdf"), plot = g, width = fig.width, height = fig.height)





g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = alpha_1, y = alpha_2, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(alpha[1])) + ylab(expression(alpha[2]))
g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
# g
ggsave(paste0(prefix, "toggleswitch_bivar12.png"), plot = g, width = 5, height = 5, dpi = 150)

g <- ggplot(wsmc.df, aes(x = beta_1, y = beta_2, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(beta[1])) + ylab(expression(beta[2]))
g <- g + geom_vline(xintercept = true_theta[3]) + geom_hline(yintercept = true_theta[4])
ggsave(paste0(prefix, "toggleswitch_bivar34.png"), plot = g, width = 5, height = 5, dpi = 150)

g <- ggplot(wsmc.df, aes(x = mu, y = sigma, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(mu)) + ylab(expression(sigma))
g <- g + geom_vline(xintercept = true_theta[5]) + geom_hline(yintercept = true_theta[6])
ggsave(paste0(prefix, "toggleswitch_bivar56.png"), plot = g, width = 5, height = 5, dpi = 150)

g <- ggplot(wsmc.df, aes(x = sigma, y = gamma, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(sigma)) + ylab(expression(gamma))
g <- g + geom_vline(xintercept = true_theta[6]) + geom_hline(yintercept = true_theta[7])
ggsave(paste0(prefix, "toggleswitch_bivar67.png"), plot = g, width = 5, height = 5, dpi = 150)


# g
