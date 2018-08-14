library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()
set.seed(11)
prefix = ""
my_colors <- get_my_colors()

target <- get_cosine()
fig.width <- 7
fig.height <- 5

nobservations <- 100

load(paste0(prefix, "cosinedata.RData"))
obs <- matrix(obs[1:nobservations], nrow = 1)
filename <- paste0(prefix, "cosine_wsmc_euclidean.n", nobservations, ".RData")
load(filename)
results_euclidean <- results
wsmc.euclidean.df <- wsmc_to_dataframe(results_euclidean)
sum(results_euclidean$ncomputed)
step_euclidean = tail(wsmc.euclidean.df$step, n = 1)



lambda <- 2
# filename <- paste0(prefix, "cosine_wsmc_curvematching.hilbert.lambda", lambda, ".n", nobservations, ".RData")
filename <- paste0(prefix, "cosine_wsmc_curvematching.swap.lambda", lambda, ".n", nobservations, ".RData")
load(filename)
results_curvematching1 <- results
wsmc.curvematching1.df <- wsmc_to_dataframe(results_curvematching1)
sum(results_curvematching1$ncomputed)
step_cm = tail(wsmc.curvematching1.df$step, n = 1)


mhfile <- paste0(prefix, "cosine.mcmc.n", nobservations, ".RData")
load(mhfile)
burnin <- 50000
chain.df <- mhchainlist_to_dataframe(mh$chains)
chain.df <- chain.df %>% filter(iteration > burnin)
names(chain.df) <- c("ichain", "iteration", target$parameter_names)

wsmc.euclidean.df %>% head
g <- ggplot(wsmc.euclidean.df %>% filter(step == step_euclidean), aes(x = omega)) + geom_density(aes(y = ..density.., fill = "Euclidean", colour = "Euclidean"), alpha = 0.5)
g <- g + geom_density(data=wsmc.curvematching1.df %>% filter(step == step_cm), aes(y = ..density.., fill = "Curve matching", colour = "Curve matching"), alpha = 0.5)
g <- g + geom_density(data=chain.df, aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
# g <- g + geom_density(data=wsmc.marginal.df %>% filter(step == step_marginal), aes(y = ..density.., fill = "marginal", colour = "marginal"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors)
g <- g + xlab(expression(omega)) #+ geom_vline(xintercept = true_theta[1])
g <- g + geom_label(data = data.frame(x = c(0.013, 0.0131,0.0117), y = c(750, 1200, 1700), method = c("Posterior", "Curve matching", "Euclidean")),
                    aes(x = x, y = y, colour = method, label = method), size = 8) + theme(legend.position = "none")
g
#ggsave(filename = paste0(prefix, "cosine_omega.pdf"), plot = g, width = fig.width, height = fig.height)
#ggsave(filename = paste0(prefix, "cosine_omega.png"), plot = g, width = fig.width, height = fig.height, dpi = 150)

g <- ggplot(wsmc.euclidean.df %>% filter(step == step_euclidean), aes(x = phi)) + geom_density(aes(y = ..density.., fill = "Euclidean", colour = "Euclidean"), alpha = 0.5)
g <- g + geom_density(data=wsmc.curvematching1.df %>% filter(step == step_cm), aes(y = ..density.., fill = "Curve matching", colour = "Curve matching"), alpha = 0.5)
g <- g + geom_density(data=chain.df, aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
# g <- g + geom_density(data=wsmc.curvematching2.df, aes(y = ..density.., fill = "CM 2", colour = "CM 2"), alpha = 0.5)
# g <- g + geom_density(data=wsmc.marginal.df %>% filter(step == step_marginal), aes(y = ..density.., fill = "marginal", colour = "marginal"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors)
g <- g + xlab(expression(phi)) # + geom_vline(xintercept = true_theta[2])
g <- g + geom_label(data = data.frame(x = c(0.65, 0.6,1), y = c(2, 3, 5), method = c("Posterior", "Curve matching", "Euclidean")),
                    aes(x = x, y = y, colour = method, label = method), size = 8) + theme(legend.position = "none")
g
#ggsave(filename = paste0(prefix, "cosine_phi.pdf"), plot = g, width = fig.width, height = fig.height)
#ggsave(filename = paste0(prefix, "cosine_phi.png"), plot = g, width = fig.width, height = fig.height, dpi = 150)

g <- ggplot(wsmc.euclidean.df%>% filter(step == step_euclidean), aes(x = logsigma)) + geom_density(aes(y = ..density.., fill = "Euclidean", colour = "Euclidean"), alpha = 0.5)
g <- g + geom_density(data=wsmc.curvematching1.df %>% filter(step == step_cm), aes(y = ..density.., fill = "Curve matching", colour = "Curve matching"), alpha = 0.5)
g <- g + geom_density(data=chain.df, aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
# g <- g + geom_density(data=wsmc.curvematching2.df, aes(y = ..density.., fill = "CM 2", colour = "CM 2"), alpha = 0.5)
# g <- g + geom_density(data=wsmc.marginal.df %>% filter(step == step_marginal), aes(y = ..density.., fill = "marginal", colour = "marginal"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors)
g <- g + xlab(expression(log(sigma))) # + geom_vline(xintercept = true_theta[3])
g <- g + geom_label(data = data.frame(x = c(-0.5, -1,-2), y = c(4,  2, 1), method = c("Posterior", "Curve matching", "Euclidean")),
                    aes(x = x, y = y, colour = method, label = method), size = 8) + theme(legend.position = "none")

g
#ggsave(filename = paste0(prefix, "cosine_logsigma.pdf"), plot = g, width = fig.width, height = fig.height)
#ggsave(filename = paste0(prefix, "cosine_logsigma.png"), plot = g, width = fig.width, height = fig.height, dpi = 150)

g <- ggplot(wsmc.euclidean.df%>% filter(step == step_euclidean), aes(x = logA)) + geom_density(aes(y = ..density.., fill = "Euclidean", colour = "Euclidean"), alpha = 0.5)
g <- g + geom_density(data=wsmc.curvematching1.df %>% filter(step == step_cm), aes(y = ..density.., fill = "Curve matching", colour = "Curve matching"), alpha = 0.5)
g <- g + geom_density(data=chain.df, aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
# g <- g + geom_density(data=wsmc.curvematching2.df, aes(y = ..density.., fill = "CM 2", colour = "CM 2"), alpha = 0.5)
# g <- g + geom_density(data=wsmc.marginal.df %>% filter(step == step_marginal), aes(y = ..density.., fill = "marginal", colour = "marginal"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors)
g <- g + xlab(expression(log(A))) # + geom_vline(xintercept = true_theta[4])
g <- g + geom_label(data = data.frame(x = c(0.6, 0.96 ,0.64), y = c(3, 6.5, 6.5), method = c("Posterior", "Curve matching", "Euclidean")),
               aes(x = x, y = y, colour = method, label = method), size = 8) + theme(legend.position = "none")
g <- g + xlim(0.4, 1.1)
g
#ggsave(filename = paste0(prefix, "cosine_logA.pdf"), plot = g, width = fig.width, height = fig.height)
#ggsave(filename = paste0(prefix, "cosine_logA.png"), plot = g, width = fig.width, height = fig.height, dpi = 150)


