library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
my_colors <- get_my_colors()
setmytheme()

set.seed(11)
prefix <- ""
regime <- "intermediate"
fig.width <- 6
fig.height <- 6.5

load(paste0(prefix, "50.", regime, "obs.neal.Rdata"))
obs <- matrix(obs, nrow = 1)
nobs <- ncol(obs)
#
target <- get_queue()
x_obs <- cumsum(obs[1,])
#
filename <- paste0(prefix, "queue_", regime, "_wsmc_marginal.RData")
load(filename)
results_wsmc <- results

filename <- paste0(prefix, "queue_", regime, "_wsmc_marginal_constraints.RData")
load(filename)
results_wsmc_constraints <- results

#
load(paste0(prefix, "queue_", regime, "_alivepmmh.RData"))
if (exists("pmmh_results")){
  pmmh.alive.df <- mhchainlist_to_dataframe(pmmh_results$chains)
  pmmh.alive.df <- pmmh.alive.df %>% filter(iteration < pmmh_results$iteration)
  rm(pmmh_results)
} else {
  pmmh.alive.df <- mhchainlist_to_dataframe(res$chains)
}
load(paste0(prefix, "queue_", regime, "_pmmh.RData"))
if (exists("pmmh_results")){
  pmmh.df <- mhchainlist_to_dataframe(pmmh_results$chains)
  pmmh.df <- pmmh.df %>% filter(iteration < pmmh_results$iteration)
} else {
  pmmh.df <- mhchainlist_to_dataframe(res$chains)
}

names(pmmh.df) <- c("ichain", "iteration", target$parameter_names)
names(pmmh.alive.df) <- c("ichain", "iteration", target$parameter_names)

# ggplot(pmmh.df %>% filter(iteration %% 100 == 1), aes(x = iteration, y = theta1, group = ichain, colour = factor(ichain))) + geom_line(alpha = 0.5)
# ggplot(pmmh.alive.df %>% filter(iteration %% 100 == 1), aes(x = iteration, y = theta1, group = ichain, colour = factor(ichain))) + geom_line(alpha = 0.5)
wsmc.wasserstein.df <- wsmc_to_dataframe(results_wsmc)
last_step <- max(wsmc.wasserstein.df$step)
wsmc.wasserstein.constraints.df <- wsmc_to_dataframe(results_wsmc_constraints)
last_step_constraints <- max(wsmc.wasserstein.constraints.df$step)

results_wsmc$ncomputed %>% sum
results_wsmc_constraints$ncomputed %>% sum

thetas <- tail(results_wsmc$thetas_history, 1)[[1]]
burnin <- 10000
g1 <- ggplot(pmmh.df %>% filter(iteration > burnin), aes(x = theta1)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
# g1 <- g1 + geom_density(data = pmmh.alive.df %>% filter(iteration > burnin), aes(y = ..density.., fill = "Posterior 2", colour = "Posterior 2"),  alpha = 0.5)
g1 <- g1 + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors)
g1 <- g1 + geom_density(data = wsmc.wasserstein.df %>% filter(step == last_step), aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"),  alpha = 0.5)
g1 <- g1 + geom_density(data = wsmc.wasserstein.constraints.df %>% filter(step == last_step_constraints), aes(y = ..density.., fill = "W + constraint", colour = "W + constraint"),  alpha = 0.5)
g1 <- g1 + xlab(expression(theta[1]))
if (regime == "intermediate"){
  # g1 <- g1 + geom_vline(xintercept = 4)
  g1 <- g1 + geom_label(data = data.frame(x = c(4.15,4.2, 3.7), y = c(9, 3, 3.5),
                                          method = c("Posterior", "Wasserstein", "W + constraint")),
                        aes(x = x, y = y, colour = method, label = method), size = 8) + theme(legend.position = "none")
}
g1
ggsave(plot = g1, filename = paste0(prefix, "queue_", regime, "_theta1.pdf"), width = fig.width, height = fig.height)
ggsave(plot = g1, filename = paste0(prefix, "queue_", regime, "_theta1.png"), width = fig.width, height = fig.height, dpi = 150)

g2 <- ggplot(pmmh.df %>% filter(iteration > burnin), aes(x = theta1+theta2minus1)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
# g2 <- g2 + geom_density(data = pmmh.alive.df %>% filter(iteration > burnin), aes(y = ..density.., fill = "Posterior 2", colour = "Posterior 2"),  alpha = 0.5)
g2 <- g2 + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors)
g2 <- g2 + geom_density(data = wsmc.wasserstein.df %>% filter(step == last_step), aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"),  alpha = 0.5)
g2 <- g2 + geom_density(data = wsmc.wasserstein.constraints.df %>% filter(step == last_step_constraints), aes(y = ..density.., fill = "W + constraint", colour = "W + constraint"),  alpha = 0.5)
g2 <- g2 + xlab(expression(theta[2]))
if (regime == "intermediate"){
  # g2 <- g2 + geom_vline(xintercept = 7)
  g2 <- g2 + geom_label(data = data.frame(x = c(7, 6.7, 6.9), y = c(9, 2, 3),
                                          method = c("Posterior", "Wasserstein", "W + constraint")),
                        aes(x = x, y = y, colour = method, label = method), size = 8) + theme(legend.position = "none")
}
ggsave(plot = g2, filename = paste0(prefix, "queue_", regime, "_theta2.pdf"), width = fig.width, height = fig.height)
ggsave(plot = g2, filename = paste0(prefix, "queue_", regime, "_theta2.png"), width = fig.width, height = fig.height, dpi = 150)

g3 <- ggplot(pmmh.df %>% filter(iteration > burnin), aes(x = theta3)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
# g3 <- g3 + geom_density(data = pmmh.alive.df %>% filter(iteration > burnin), aes(y = ..density.., fill = "Posterior 2", colour = "Posterior 2"),  alpha = 0.5)
g3 <- g3 + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors)
g3 <- g3 + geom_density(data = wsmc.wasserstein.df %>% filter(step == last_step), aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"),  alpha = 0.5)
g3 <- g3 + geom_density(data = wsmc.wasserstein.constraints.df %>% filter(step == last_step_constraints), aes(y = ..density.., fill = "W + constraint", colour = "W + constraint"),  alpha = 0.5)
g3 <- g3 + xlab(expression(theta[3]))
if (regime == "intermediate"){
  # g3 <- g3 + geom_vline(xintercept = 0.15)
  g3 <- g3 + geom_label(data = data.frame(x = c(0.12, 0.125, 0.235), y = c(15, 13, 15),
                                          method = c("Posterior", "Wasserstein", "W + constraint")),
                        aes(x = x, y = y, colour = method, label = method), size = 8) + theme(legend.position = "none")
}
ggsave(plot = g3, filename = paste0(prefix, "queue_", regime, "_theta3.png"), width = fig.width, height = fig.height, dpi = 150)
ggsave(plot = g3, filename = paste0(prefix, "queue_", regime, "_theta3.pdf"), width = fig.width, height = fig.height)


thetas <- tail(results_wsmc$thetas_history, 1)[[1]]
colMeans(thetas)
sqrt(diag(cov(thetas)))

