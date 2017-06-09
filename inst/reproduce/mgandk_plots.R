library(winference)
registerDoParallel(cores = 6)
rm(list = ls())
setmytheme()
set.seed(11)
target <- get_mgandk()
fig.height <- 5
fig.width <- 5

my_colors <- get_my_colors()
prefix <- ""

nobservations <- 500
load(paste0(prefix, "mgandkdata.RData"))
obs <- obs[,1:nobservations]
target$simulate <- function(theta) target$robservation(nobservations, theta)

compute_d <- get_hilbert_to_y(obs)
y_sim <- target$simulate(target$rprior(1, target$parameters)[1,])
compute_d(y_sim)

param_algo <- list(nthetas = 1024, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 100000)

filename <- paste0(prefix, "mgandk.wsmc.n", nobservations, ".hilbert.RData")
load(filename)
results_hilbert <- results
wsmc.hilbert.df <- wsmc_to_dataframe(results_hilbert)
sum(results_hilbert$ncomputed)
step_hilbert <- length(results_hilbert$thetas_history)
step_hilbert <- which(cumsum(results_hilbert$ncomputed) > 2e6)[1]
#
filename <- paste0(prefix, "mgandk.wsmc.n", nobservations, ".mmd.RData")
load(filename)
results_mmd <- results
wsmc.mmd.df <- wsmc_to_dataframe(results_mmd)
sum(results_mmd$ncomputed)
step_mmd <- length(results_mmd$thetas_history)
step_mmd <- which(cumsum(results_mmd$ncomputed) > 2e6)[1]
#
filename <- paste0(prefix, "mgandk.wsmc.n", nobservations, ".swap.RData")
load(filename)
results_swap <- results
wsmc.swap.df <- wsmc_to_dataframe(results_swap)
sum(results_swap$ncomputed)
step_swap <- length(results_swap$thetas_history)
step_swap <- which(cumsum(results_swap$ncomputed) > 2e6)[1]
#
filename <- paste0(prefix, "mgandk.wsmc.n", nobservations, ".sinkhorn.RData")
load(filename)
results_sinkhorn <- results
wsmc.sinkhorn.df <- wsmc_to_dataframe(results_sinkhorn)
sum(results_sinkhorn$ncomputed)
step_sinkhorn <- length(results_sinkhorn$thetas_history)
step_sinkhorn <- which(cumsum(results_sinkhorn$ncomputed) > 2e6)[1]
# plot_marginal_time(results_sinkhorn, 7)

#
filename <- paste0(prefix, "mgandk.wsmc.n", nobservations, ".wasserstein.RData")
load(filename)
results_wasserstein <- results
wsmc.wasserstein.df <- wsmc_to_dataframe(results_wasserstein)
sum(results_wasserstein$ncomputed)
step_wasserstein <- length(results_wasserstein$thetas_history)
step_wasserstein <- which(cumsum(results_wasserstein$ncomputed) > 2e6)[1]

#
mhfile <- paste0(prefix, "mgandk.mcmc.n", nobservations, ".mh.initfromABC.RData")
load(mhfile)
if (exists("mh_results")){
  mcmc.df <- mhchainlist_to_dataframe(mh_results$chains)
  mcmc.df <- mcmc.df %>% filter(iteration < mh_results$iteration)
} else {
  mcmc.df <- mhchainlist_to_dataframe(mh$chains)
}
burnin <- 10000
names(mcmc.df) <- c("ichain", "iteration", target$parameter_names)
target$parameter_names

g <- ggplot(mcmc.df %>% filter(iteration > burnin), aes(x = a1)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
g <- g + geom_density(data = wsmc.wasserstein.df %>% filter(step == step_wasserstein), aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab(expression(a[1]))
g <- g + geom_label(data = data.frame(x = c(2.9, 3.1), y = c(4, 6), method = c("Wasserstein", "Posterior")),
                    aes(x = x, y = y, colour = method, label = method), size = 8) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "mgandk.n", nobservations, ".a1.pdf"), plot = g, width = fig.width, height = fig.height)

g <- ggplot(mcmc.df %>% filter(iteration > burnin), aes(x = b1)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
g <- g + geom_density(data = wsmc.wasserstein.df %>% filter(step == step_wasserstein), aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab(expression(b[1]))
g <- g + geom_label(data = data.frame(x = c(1.2, 0.9), y = c(3, 5), method = c("Wasserstein", "Posterior")),
                    aes(x = x, y = y, colour = method, label = method), size = 8) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "mgandk.n", nobservations, ".b1.pdf"), plot = g, width = fig.width, height = fig.height)


g <- ggplot(mcmc.df %>% filter(iteration > burnin), aes(x = g1)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
g <- g + geom_density(data = wsmc.wasserstein.df %>% filter(step == step_wasserstein), aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab(expression(g[1]))
g <- g + geom_label(data = data.frame(x = c(1.25, 0.9), y = c(2, 4), method = c("Wasserstein", "Posterior")),
                    aes(x = x, y = y, colour = method, label = method), size = 8) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "mgandk.n", nobservations, ".g1.pdf"), plot = g, width = fig.width, height = fig.height)

g <- ggplot(mcmc.df %>% filter(iteration > burnin), aes(x = k1)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
g <- g + geom_density(data = wsmc.wasserstein.df %>% filter(step == step_wasserstein), aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab(expression(k[1]))
g <- g + geom_label(data = data.frame(x = c(0.4, 0.6), y = c(3, 5), method = c("Wasserstein", "Posterior")),
                    aes(x = x, y = y, colour = method, label = method), size = 8) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "mgandk.n", nobservations, ".k1.pdf"), plot = g, width = fig.width, height = fig.height)

g <- ggplot(mcmc.df %>% filter(iteration > burnin), aes(x = a2)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
g <- g + geom_density(data = wsmc.wasserstein.df %>% filter(step == step_wasserstein), aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab(expression(a[2]))
g <- g + geom_label(data = data.frame(x = c(3.95, 4.05), y = c(5, 10), method = c("Wasserstein", "Posterior")),
aes(x = x, y = y, colour = method, label = method), size = 8) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "mgandk.n", nobservations, ".a2.pdf"), plot = g, width = fig.width, height = fig.height)

g <- ggplot(mcmc.df %>% filter(iteration > burnin), aes(x = b2)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
g <- g + geom_density(data = wsmc.wasserstein.df %>% filter(step == step_wasserstein), aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab(expression(b[2]))
g <- g + geom_label(data = data.frame(x = c(0.6, 0.4), y = c(3, 7), method = c("Wasserstein", "Posterior")),
                    aes(x = x, y = y, colour = method, label = method), size = 8) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "mgandk.n", nobservations, ".b2.pdf"), plot = g, width = fig.width, height = fig.height)


g <- ggplot(mcmc.df %>% filter(iteration > burnin), aes(x = g2)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
g <- g + geom_density(data = wsmc.wasserstein.df %>% filter(step == step_wasserstein), aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab(expression(g[2]))
g <- g + geom_label(data = data.frame(x = c(6, 2.5), y = c(0.5, 3), method = c("Wasserstein", "Posterior")),
                    aes(x = x, y = y, colour = method, label = method), size = 8) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "mgandk.n", nobservations, ".g2.pdf"), plot = g, width = fig.width, height = fig.height)
plot_marginal(results_wasserstein, 7, from_step = 50) + geom_vline(xintercept = true_theta[7])


g <- ggplot(mcmc.df %>% filter(iteration > burnin), aes(x = k2)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
g <- g + geom_density(data = wsmc.wasserstein.df %>% filter(step == step_wasserstein), aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab(expression(k[2]))
g <- g + geom_label(data = data.frame(x = c(0.3, 0.5), y = c(3, 7), method = c("Wasserstein", "Posterior")),
                    aes(x = x, y = y, colour = method, label = method), size = 8) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "mgandk.n", nobservations, ".k2.pdf"), plot = g, width = fig.width, height = fig.height)

g <- ggplot(mcmc.df %>% filter(iteration > burnin), aes(x = rho)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
g <- g + geom_density(data = wsmc.wasserstein.df %>% filter(step == step_wasserstein), aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g <- g + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors) + xlab(expression(rho))
g <- g + geom_label(data = data.frame(x = c(0.45, 0.65), y = c(3, 10), method = c("Wasserstein", "Posterior")),
                    aes(x = x, y = y, colour = method, label = method), size = 8) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "mgandk.n", nobservations, ".rho.pdf"), plot = g, width = fig.width, height = fig.height)

## plot for g2 only
g <- ggplot(wsmc.wasserstein.df %>% filter(step > 70), aes(x = g2, colour = step, group = step))
g <- g + geom_density(aes(y = ..density..))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue") + theme(legend.position = "none")
g <- g + xlab(expression(g[2]))
g <- g + geom_vline(xintercept = true_theta[7]) + geom_segment(aes(x = 0, xend = 10, y = 0.1, yend = 0.1), colour = rgb(1,0.5,0.5))
g
ggsave(filename = paste0(prefix, "mgandk.n", nobservations, ".g2evolution.pdf"), plot = g, width = fig.width, height = fig.height)



# compute W distance between approx and posterior
wasserstein_to_posterior <- function(results, posterior_sample, nmax){
  nsteps <- length(results$thetas_history)
  if (nsteps > nmax){
    steps <- floor(seq(from = 1, to = nsteps, length.out = nmax))
  } else {
    steps <- 1:nsteps
  }
  w_to_post_ <- as.numeric(foreach (istep = steps, .combine = c) %dorng% {
    thetas <- results$thetas_history[[istep]]
    x <- try(exact_transport_distance(t(thetas), t(posterior_sample), p = 1, ground_p = 2))
    if (inherits(x, "try-error")){
      x <- sinkhorn_distance(t(thetas), t(posterior_sample), p = 1, ground_p = 2, eps = 0.1, niterations = 1000)$uncorrected
    }
    x
  })
  return(data.frame(steps = steps, ncomputed = cumsum(results$ncomputed)[steps], times = results$compute_times[steps], w = w_to_post_))
}

npost <- 1024
mcmcpostburnin <- mcmc.df %>% filter(iteration > burnin)
posterior_sample <- mcmcpostburnin[sample(x = 1:nrow(mcmcpostburnin), size = npost, replace = TRUE),3:ncol(mcmcpostburnin)]

filename <- paste0(prefix, "mgandk.n", nobservations, ".w2posterior.RData")
n_w_comp <- 50
w_to_post_wasserstein <- wasserstein_to_posterior(results_wasserstein, posterior_sample, n_w_comp)
w_to_post_hilbert <- wasserstein_to_posterior(results_hilbert, posterior_sample, n_w_comp)
w_to_post_swap <- wasserstein_to_posterior(results_swap, posterior_sample, n_w_comp)
w_to_post_mmd <- wasserstein_to_posterior(results_mmd, posterior_sample, n_w_comp)
w_to_post_sinkhorn <- wasserstein_to_posterior(results_sinkhorn, posterior_sample, n_w_comp)
save(w_to_post_wasserstein, w_to_post_hilbert, w_to_post_swap, w_to_post_mmd, w_to_post_sinkhorn, file = filename)
load(file = filename)

w_to_post_wasserstein %>% tail

g <- ggplot(w_to_post_wasserstein, aes(x = ncomputed, y = w)) + geom_line(aes(colour = "Wasserstein")) + geom_point(aes(colour = "Wasserstein"))
g <- g + geom_line(data = w_to_post_hilbert, aes(colour = "Hilbert")) + geom_point(data = w_to_post_hilbert, aes(colour = "Hilbert"))
g <- g + geom_line(data = w_to_post_swap, aes(colour = "Swap")) + geom_point(data = w_to_post_swap, aes(colour = "Swap"))
g <- g + geom_line(data = w_to_post_mmd, aes(colour = "MMD")) + geom_point(data = w_to_post_mmd, aes(colour = "MMD"))
# g <- g + geom_line(data = w_to_post_sinkhorn, aes(colour = "Sinkhorn")) + geom_point(data = w_to_post_sinkhorn, aes(colour = "Sinkhorn"))
g <- g + scale_color_manual(name = "", values = my_colors) + xlab("# model simulations") + ylab("W-distance to posterior")
g <- g + scale_x_log10(breaks = c(1e4, 1e5, 1e6), limits = c(1e4, 1e7)) + scale_y_log10(breaks = c(3, 5,10))
# g <- g + geom_hline(yintercept = self.distances.qt[1], linetype = 2) + geom_hline(yintercept = self.distances.qt[2], linetype = 2)
g <- g + geom_label(data = data.frame(x = c(3e5, 7e6, 5e6, 3e6), y = c(2.75, 3.4, 2.75, 4.1), method = c("Wasserstein", "Hilbert", "Swap", "MMD")),
aes(x = x, y = y, colour = method, label = method), size = 8) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "mgandk.n", nobservations, ".w2post.pdf"), plot = g, width = 2*fig.width, height = fig.height)
