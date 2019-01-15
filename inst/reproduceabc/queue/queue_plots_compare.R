library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()
set.seed(11)

prefix = ""

target <- get_queue()

nobservations = 50
load(paste0(prefix, "50.intermediateobs.neal.Rdata"))
true_theta = c(4,7,0.15)
obs <- matrix(obs, nrow = 1)
obs_sorted = sort(obs)
nobs <- ncol(obs)

fig.height <- 5
fig.width <- 5

my_colors <- get_my_colors()

# Wasserstein SMC
filename = paste0(prefix, "queue_intermediate_wsmc_marginal.RData")
load(filename)
results_was = results
wsmc.df = wsmc_to_dataframe(results_was)
wsmc.df$theta2 = wsmc.df$theta2minus1 + wsmc.df$theta1
nsteps = tail(wsmc.df$step,n=1)
#step = nsteps
step = 26


# Wasserstein SMC with constraint
filename = paste0(prefix, "queue_intermediate_wsmc_marginal_constraints.RData")
load(filename)
results_cwas = results
cwsmc.df = wsmc_to_dataframe(results_cwas)
cwsmc.df$theta2 = cwsmc.df$theta2minus1 + cwsmc.df$theta1
cnsteps = tail(cwsmc.df$step,n=1)
#cstep = cnsteps
cstep = 27


# abctools package
#nsim = cumsum(results_was$ncomputed)[step]
nsim = 10^7
m = 20; l = 1
filename = paste0(prefix, "queueing.abctools.n", nobservations,".m",m,".l",l,".nsim",nsim,".RData")
load(filename)
results_abctools <- results
abctools.df = data.frame(results_abctools$thetas)
#abctools.df[,2] = abctools.df[,2] - abctools.df[,1] #alternative
col.names = colnames(wsmc.df)[1:3]
colnames(abctools.df) = col.names
#abctools.df$theta2 = abctools.df$theta2minus1 + abctools.df$theta1


# # abctools package with constraint
# nsim = 10^7
# m = 20; l = 1
# filename = paste0(prefix, "constrained.queueing.abctools.n", nobservations,".m",m,".l",l,".nsim",nsim,".RData")
# load(filename)
# results_abctools.const <- results
# abctools.const.df = data.frame(results_abctools.const$thetas)
# col.names = colnames(wsmc.df)[1:3]
# colnames(abctools.const.df) = col.names
# abctools.const.df$theta2 = abctools.const.df$theta2minus1 + abctools.const.df$theta1


# PMMH
filename = paste0(prefix, "queue_intermediate_pmmh.RData")
load(file = filename)
results_pmmh = res
mcmc.df = mhchainlist_to_dataframe(results_pmmh$chains)
mcmc.df = mcmc.df %>% filter(iteration %% 10 == 1)
names(mcmc.df) = c("ichain","iteration","theta1","theta2minus1","theta3")
mcmc.df$theta2 = mcmc.df$theta2minus1 + mcmc.df$theta1



prefix = "compare."

g <- ggplot(mcmc.df, aes(x = theta1)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
g <- g + geom_density(data = wsmc.df %>% filter(step == step), aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + geom_density(data = cwsmc.df %>% filter(step == cstep), aes(y = ..density.., fill = "W + constraint", colour = "W + constraint"), alpha = 0.5)
g <- g + geom_density(data = abctools.df, aes(y = ..density.., fill = "Semi-auto", colour = "Semi-auto"), alpha = 0.5)
#g <- g + geom_density(data = abctools.df, aes(y = ..density.., fill = "SA + constraint", colour = "SA + constraint"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab(expression(theta[1]))
g <- g + geom_label(data = data.frame(x = c(4.37, 3.57, 3.75, 3.4), y = c(1.8, 3.8, 8, 1.75), method = c("Wasserstein","W + constraint", "Posterior","Semi-auto")),
                    aes(x = x, y = y, colour = method, label = method), size = 6) + theme(legend.position = "none")
g <- g + coord_cartesian(xlim = c(3.2, 4.6))
g <- g + geom_vline(xintercept = true_theta[1])
g
ggsave(filename = paste0(prefix, "queueing_marginal1.pdf"), plot = g, width = fig.width, height = fig.height)

# g <- ggplot(mcmc.df, aes(x = theta2)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
# g <- g + geom_density(data = wsmc.df %>% filter(step == step), aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
# g <- g + geom_density(data = cwsmc.df %>% filter(step == cstep), aes(y = ..density.., fill = "W + constraint", colour = "W + constraint"), alpha = 0.5)
# g <- g + geom_density(data = abctools.df, aes(y = ..density.., fill = "Semi-auto", colour = "Semi-auto"), alpha = 0.5)
# g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab(expression(theta[2]))
# #g <- g + geom_label(data = data.frame(x = c(5.3, 5.5, 4.7, 8.3), y = c(0.7, 1.5, 4, 0.3), method = c("Wasserstein","W + constraint", "Posterior","Semi-auto")),
# #                    aes(x = x, y = y, colour = method, label = method), size = 6) + theme(legend.position = "none")
# g <- g + geom_vline(xintercept = true_theta[2])
# g
# ggsave(filename = paste0(prefix, "queueing_marginal2.pdf"), plot = g, width = fig.width, height = fig.height)


g <- ggplot(mcmc.df, aes(x = theta2minus1)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
g <- g + geom_density(data = wsmc.df %>% filter(step == step), aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + geom_density(data = cwsmc.df %>% filter(step == cstep), aes(y = ..density.., fill = "W + constraint", colour = "W + constraint"), alpha = 0.5)
g <- g + geom_density(data = abctools.df, aes(y = ..density.., fill = "Semi-auto", colour = "Semi-auto"), alpha = 0.5)
#g <- g + geom_density(data = abctools.df, aes(y = ..density.., fill = "SA + constraint", colour = "SA + constraint"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab(expression(theta[2]-theta[1]))
g <- g + geom_label(data = data.frame(x = c(1.85, 4.3, 3.7, 5.3), y = c(0.85, 1.5, 4, 0.5), method = c("Wasserstein","W + constraint", "Posterior","Semi-auto")),
                    aes(x = x, y = y, colour = method, label = method), size = 6) + theme(legend.position = "none")
g <- g + geom_vline(xintercept = true_theta[2]-true_theta[1])
g <- g + coord_cartesian(xlim = c(1.2, 6))
g
ggsave(filename = paste0(prefix, "queueing_marginal2minus1.pdf"), plot = g, width = fig.width, height = fig.height)


g <- ggplot(mcmc.df, aes(x = theta3)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
g <- g + geom_density(data = wsmc.df %>% filter(step == step), aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + geom_density(data = cwsmc.df %>% filter(step == cstep), aes(y = ..density.., fill = "W + constraint", colour = "W + constraint"), alpha = 0.5)
g <- g + geom_density(data = abctools.df, aes(y = ..density.., fill = "Semi-auto", colour = "Semi-auto"), alpha = 0.5)
#g <- g + geom_density(data = abctools.df, aes(y = ..density.., fill = "SA + constraint", colour = "SA + constraint"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab(expression(theta[3]))
g <- g + geom_label(data = data.frame(x = c(0.072, 0.073, 0.25, 0.29), y = c(8, 14, 15, 5.3), method = c("Wasserstein","W + constraint", "Posterior","Semi-auto")),
                    aes(x = x, y = y, colour = method, label = method), size = 6) + theme(legend.position = "none")
g <- g + geom_vline(xintercept = true_theta[3])
g
ggsave(filename = paste0(prefix, "queueing_marginal3.pdf"), plot = g, width = fig.width, height = fig.height)
