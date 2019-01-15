library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()
set.seed(11)
target <- get_gandk()

fig.height <- 5
fig.width <- 5

my_colors <- get_my_colors()

prefix <- ""

nobservations = 250
load(paste0(prefix, "gandkdata.RData"))
obs <- obs[1:nobservations]

# Wasserstein SMC
filename = paste0(prefix, "gandkwsmc.n", nobservations, ".RData")
load(filename)
results_was = results
wsmc.df = wsmc_to_dataframe(results_was)
nsteps = tail(wsmc.df$step,n=1)
#step = nsteps
step = 20
wsmc.df = wsmc.df[wsmc.df$step>=step,]


# MCMC
mhfile <- paste0(prefix, "gandkmcmc.n", nobservations, "mh.RData")
load(mhfile)
mcmc.df <- mhchainlist_to_dataframe(mh$chains)
names(mcmc.df) <- c("ichain", "iteration", target$parameter_names)
burnin = 50000
mcmc.df = mcmc.df[mcmc.df$iteration>burnin,3:6]


g <- ggplot(mcmc.df, aes(x = A)) + geom_density(aes(y = ..density.., fill = "Posterior"), alpha = 0.5)
g <- g + geom_density(data = wsmc.df, aes(x = A, colour = step, group = step), alpha = 0.5)
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue", name = "")
g <- g + scale_fill_manual(name = "", values = my_colors) + xlab("a") + theme(legend.position = "none")
g <- g + geom_vline(xintercept = true_theta[1])
g
ggsave(filename = paste0(prefix, "conv.gandk_marginal1.pdf"), plot = g, width = fig.width, height = fig.height)

g <- ggplot(mcmc.df, aes(x = B)) + geom_density(aes(y = ..density.., fill = "Posterior"), alpha = 0.5)
g <- g + geom_density(data = wsmc.df, aes(x = B, colour = step, group = step), alpha = 0.5)
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue", name = "")
g <- g + scale_fill_manual(name = "", values = my_colors) + xlab("b") + theme(legend.position = "none")
g <- g + geom_vline(xintercept = true_theta[2])
g
ggsave(filename = paste0(prefix, "conv.gandk_marginal2.pdf"), plot = g, width = fig.width, height = fig.height)

g <- ggplot(mcmc.df, aes(x = g)) + geom_density(aes(y = ..density.., fill = "Posterior"), alpha = 0.5)
g <- g + geom_density(data = wsmc.df, aes(x = g, colour = step, group = step), alpha = 0.5)
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue", name = "")
g <- g + scale_fill_manual(name = "", values = my_colors) + xlab("g") + theme(legend.position = "none")
g <- g + geom_vline(xintercept = true_theta[3])
g <- g + coord_cartesian(xlim = c(1, 6))
g
ggsave(filename = paste0(prefix, "conv.gandk_marginal3.pdf"), plot = g, width = fig.width, height = fig.height)

g <- ggplot(mcmc.df, aes(x = k)) + geom_density(aes(y = ..density.., fill = "Posterior"), alpha = 0.5)
g <- g + geom_density(data = wsmc.df, aes(x = k, colour = step, group = step), alpha = 0.5)
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue", name = "")
g <- g + scale_fill_manual(name = "", values = my_colors) + xlab("k") + theme(legend.position = "none")
g <- g + geom_vline(xintercept = true_theta[4])
g
ggsave(filename = paste0(prefix, "conv.gandk_marginal4.pdf"), plot = g, width = fig.width, height = fig.height)
