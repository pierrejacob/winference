library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()
set.seed(11)

fig.height <- 5
fig.width <- 5

my_colors <- get_my_colors()

prefix = ""

target <- get_toggleswitch()
# number of observations
nobservations <- 2000
load(file = paste0(prefix,"toggleswitchdata.RData"))
obs <- obs[1:nobservations]
obs_sorted <- sort(obs)

g <- qplot(x = obs, geom = "blank") + geom_histogram(aes(y = ..density..))
g <- g + xlab("y")
g
ggsave(filename = paste0(prefix, "toggleswitch_data.pdf"), plot = g, width = fig.width, height = fig.height)
#ggsave(filename = paste0(prefix, "toggleswitch_data.png"), plot = g, width = fig.width, height = fig.height, dpi = 150)

# Wasserstein SMC
filename = paste0(prefix, "toggleswitchwsmc.n", nobservations, ".RData")
load(filename)
results_was = results
wsmc.df = wsmc_to_dataframe(results_was)
nsteps = tail(wsmc.df$step,n=1)
step = nsteps
#step = 33
wsmc.df = wsmc.df[wsmc.df$step == step,]


# Wasserstein SMC with constraint
filename = paste0(prefix,"toggleswitchw.summary.smc.n", nobservations, ".RData")
load(filename)
results_summary = results
summary.df = wsmc_to_dataframe(results_summary)
snsteps = tail(summary.df$step,n=1)
sstep = snsteps
#sstep = 27
summary.df = summary.df[summary.df$step == sstep,]

g <- ggplot(data = wsmc.df, aes(x = alpha_1)) + geom_density(aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + geom_density(data = summary.df, aes(y = ..density.., fill = "Summary", colour = "Summary"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab(expression(alpha[1]))
g <- g + geom_label(data = data.frame(x = c(22.95, 21.1), y = c(1.1, 0.5), method = c("Wasserstein","Summary")),
                    aes(x = x, y = y, colour = method, label = method), size = 7) + theme(legend.position = "none")
g <- g + geom_vline(xintercept = true_theta[1])
g
ggsave(filename = paste0(prefix, "compare.toggleswitch_marginal1.pdf"), plot = g, width = fig.width, height = fig.height)


g <- ggplot(data = wsmc.df, aes(x = alpha_2)) + geom_density(aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + geom_density(data = summary.df, aes(y = ..density.., fill = "Summary", colour = "Summary"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab(expression(alpha[2]))
g <- g + geom_vline(xintercept = true_theta[2])
g <- g + geom_label(data = data.frame(x = c(35, 18), y = c(0.022, 0.037), method = c("Wasserstein","Summary")),
                    aes(x = x, y = y, colour = method, label = method), size = 7) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "compare.toggleswitch_marginal2.pdf"), plot = g, width = fig.width, height = fig.height)


g <- ggplot(data = wsmc.df, aes(x = beta_1)) + geom_density(aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + geom_density(data = summary.df, aes(y = ..density.., fill = "Summary", colour = "Summary"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab(expression(beta[1]))
g <- g + geom_vline(xintercept = true_theta[3])
g <- g + geom_label(data = data.frame(x = c(4.2, 4.45), y = c(0.75, 0.55), method = c("Wasserstein","Summary")),
                    aes(x = x, y = y, colour = method, label = method), size = 7) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "compare.toggleswitch_marginal3.pdf"), plot = g, width = fig.width, height = fig.height)


g <- ggplot(data = wsmc.df, aes(x = beta_2)) + geom_density(aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + geom_density(data = summary.df, aes(y = ..density.., fill = "Summary", colour = "Summary"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab(expression(beta[2]))
g <- g + geom_label(data = data.frame(x = c(2.5, 3.3), y = c(0.35, 0.5), method = c("Wasserstein","Summary")),
                    aes(x = x, y = y, colour = method, label = method), size = 7) + theme(legend.position = "none")
g <- g + geom_vline(xintercept = true_theta[4])
g
ggsave(filename = paste0(prefix, "compare.toggleswitch_marginal4.pdf"), plot = g, width = fig.width, height = fig.height)


g <- ggplot(data = wsmc.df, aes(x = mu)) + geom_density(aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + geom_density(data = summary.df, aes(y = ..density.., fill = "Summary", colour = "Summary"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab(expression(mu))
g <- g + geom_label(data = data.frame(x = c(293, 290), y = c(0.045, 0.019), method = c("Wasserstein","Summary")),
                    aes(x = x, y = y, colour = method, label = method), size = 7) + theme(legend.position = "none")
g <- g + geom_vline(xintercept = true_theta[5])
g
ggsave(filename = paste0(prefix, "compare.toggleswitch_marginal5.pdf"), plot = g, width = fig.width, height = fig.height)


g <- ggplot(data = wsmc.df, aes(x = sigma)) + geom_density(aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + geom_density(data = summary.df, aes(y = ..density.., fill = "Summary", colour = "Summary"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab(expression(sigma))
g <- g + geom_label(data = data.frame(x = c(0.125, 0.13), y = c(15, 4.8), method = c("Wasserstein","Summary")),
                    aes(x = x, y = y, colour = method, label = method), size = 7) + theme(legend.position = "none")
g <- g + geom_vline(xintercept = true_theta[6])
g
ggsave(filename = paste0(prefix, "compare.toggleswitch_marginal6.pdf"), plot = g, width = fig.width, height = fig.height)


g <- ggplot(data = wsmc.df, aes(x = gamma)) + geom_density(aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + geom_density(data = summary.df, aes(y = ..density.., fill = "Summary", colour = "Summary"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab(expression(gamma))
g <- g + geom_label(data = data.frame(x = c(0.28, 0.3), y = c(15, 4.7), method = c("Wasserstein","Summary")),
                    aes(x = x, y = y, colour = method, label = method), size = 7) + theme(legend.position = "none")
g <- g + geom_vline(xintercept = true_theta[7])
g
ggsave(filename = paste0(prefix, "compare.toggleswitch_marginal7.pdf"), plot = g, width = fig.width, height = fig.height)
