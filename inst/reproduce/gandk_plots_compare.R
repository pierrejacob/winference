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
step = 34
wsmc.df = wsmc.df[wsmc.df$step == step,]

# abctools package
nsim = cumsum(results_was$ncomputed)[step]
m = 10; l = 4
filename = paste0(prefix, "gandk.abctools.n", nobservations,".m",m,".l",l,".nsim",nsim,".RData")
load(filename)
results_abctools <- results
abctools.df = data.frame(results_abctools$thetas)
col.names = colnames(wsmc.df)[1:4]
colnames(abctools.df) = col.names


# MCMC
mhfile <- paste0(prefix, "gandkmcmc.n", nobservations, "mh.RData")
load(mhfile)
mcmc.df <- mhchainlist_to_dataframe(mh$chains)
names(mcmc.df) <- c("ichain", "iteration", target$parameter_names)
burnin = 50000
mcmc.df = mcmc.df[mcmc.df$iteration>burnin,3:6]


prefix <- paste0(prefix,"compare.")

g <- ggplot(mcmc.df, aes(x = A)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
g <- g + geom_density(data = wsmc.df, aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + geom_density(data = abctools.df, aes(y = ..density.., fill = "Semi-auto", colour = "Semi-auto"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab("a")
g <- g + geom_label(data = data.frame(x = c(2.63, 2.73, 3.35), y = c(3.9, 5.5, 1.6), method = c("Wasserstein", "Posterior","Semi-auto")),
                    aes(x = x, y = y, colour = method, label = method), size = 7) + theme(legend.position = "none")
g <- g + geom_vline(xintercept = true_theta[1])
g
#ggsave(filename = paste0(prefix, "gandk_marginal1.pdf"), plot = g, width = fig.width, height = fig.height)


g <- ggplot(mcmc.df, aes(x = B)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
g <- g + geom_density(data = wsmc.df %>% filter(step == step), aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + geom_density(data = abctools.df %>% filter(step == step), aes(y = ..density.., fill = "Semi-auto", colour = "Semi-auto"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab("b")
g <- g + geom_label(data = data.frame(x = c(1.85, 1.53, 1.9), y = c(1.8, 2.7, 0.6), method = c("Wasserstein", "Posterior","Semi-auto")),
                    aes(x = x, y = y, colour = method, label = method), size = 7) + theme(legend.position = "none")
g <- g + geom_vline(xintercept = true_theta[2])
g
#ggsave(filename = paste0(prefix, "gandk_marginal2.pdf"), plot = g, width = fig.width, height = fig.height)



g <- ggplot(mcmc.df, aes(x = g)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
g <- g + geom_density(data = wsmc.df %>% filter(step == step), aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + geom_density(data = abctools.df %>% filter(step == step), aes(y = ..density.., fill = "Semi-auto", colour = "Semi-auto"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab("g")
g <- g + geom_label(data = data.frame(x = c(5.05, 4.4, 5.5), y = c(0.6, 1.5, 0.3), method = c("Wasserstein", "Posterior","Semi-auto")),
                    aes(x = x, y = y, colour = method, label = method), size = 7) + theme(legend.position = "none")
g <- g + geom_vline(xintercept = true_theta[3])
g
#ggsave(filename = paste0(prefix, "gandk_marginal3.pdf"), plot = g, width = fig.width, height = fig.height)


g <- ggplot(mcmc.df, aes(x = k)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
g <- g + geom_density(data = wsmc.df %>% filter(step == step), aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + geom_density(data = abctools.df %>% filter(step == step), aes(y = ..density.., fill = "Semi-auto", colour = "Semi-auto"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab("k")
g <- g + geom_label(data = data.frame(x = c(1.6, 1.3, 2), y = c(3, 4.6, 0.65), method = c("Wasserstein", "Posterior","Semi-auto")),
                    aes(x = x, y = y, colour = method, label = method), size = 7) + theme(legend.position = "none")
g <- g + geom_vline(xintercept = true_theta[4])
g <- g + coord_cartesian(xlim = c(0, 4)) 
g
#ggsave(filename = paste0(prefix, "gandk_marginal4.pdf"), plot = g, width = fig.width, height = fig.height)


