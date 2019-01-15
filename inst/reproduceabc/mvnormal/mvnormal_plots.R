library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()

my_colors <- get_my_colors()

set.seed(12)

doRun <- FALSE
max_time <- 30*60
d <- 2
target <- get_multivariate_normal(d)
target$parameters$tau <- 5
nobservations <- 100
nparticles = 2048
p <- 1
prefix <- ""

obsfile <- paste0(prefix, "mvnormaldata.d", d, ".n", nobservations, ".RData")
load(obsfile)

#exact posterior
mu_prior <- rep(target$parameters$mu_0, d)
sigma_prior <- diag(target$parameters$tau^2, d, d)
target$parameters$S
sigma_posterior <- solve(nobservations * solve(target$parameters$S) + solve(sigma_prior))
mu_posterior <- sigma_posterior %*% (solve(sigma_prior) %*% matrix(mu_prior, ncol = 1) + nobservations * solve(target$parameters$S) %*% matrix(rowMeans(obs), ncol = 1))

# exact draws from posterior
posterior_sample <- fast_rmvnorm(nparticles, mu_posterior, sigma_posterior)

# compute distance from posterior sample to posterior sample
# filename <- paste0(prefix, "mvnormal.d", d, ".n", nobservations, ".selfdistance.RData")
# self.distances <- as.numeric(foreach(i = 1:100) %dorng% {
#   x <- exact_transport_distance(t(fast_rmvnorm(nparticles, mu_posterior, sigma_posterior)), t(fast_rmvnorm(nparticles, mu_posterior, sigma_posterior)), p = 1, ground_p = 2)
#   x
# })
# save(self.distances, file = filename)
# load(filename)
# self.distances %>% summary
# self.distances.qt <- quantile(self.distances, prob = c(0.05, 0.95))

filename <- paste0(prefix, "mvnormalwsmc.d", d, ".n", nobservations, ".wasserstein.RData")
load(filename)
results_wasserstein <- results

filename <- paste0(prefix, "mvnormalwsmc.d", d, ".n", nobservations, ".euclidean.RData")
load(filename)
results_euclidean <- results

filename <- paste0(prefix, "mvnormalwsmc.d", d, ".n", nobservations, ".summary.RData")
load(filename)
results_summary <- results

filename <- paste0(prefix, "mvnormalrejection.d", d, ".n", nobservations, ".summary.RData")
load(filename)
results_rejection_summary <- results

filename <- paste0(prefix, "mvnormalrejection.d", d, ".n", nobservations, ".wasserstein.RData")
load(filename)
results_rejection_wasserstein <- results

sum(results_euclidean$ncomputed)
sum(results_summary$ncomputed)
sum(results_wasserstein$ncomputed)


wsmc.euclidean.df <- wsmc_to_dataframe(results_euclidean) %>% filter(step == length(results_euclidean$thetas_history))
wsmc.summary.df <- wsmc_to_dataframe(results_summary) %>% filter(step == length(results_summary$thetas_history))
wsmc.wasserstein.df <- wsmc_to_dataframe(results_wasserstein) %>% filter(step == length(results_wasserstein$thetas_history))
rej.summary.df = data.frame(X1 = results_rejection_summary[,1], X2 = results_rejection_summary[,2])
rej.wasserstein.df = data.frame(X1 = results_rejection_wasserstein[,1], X2 = results_rejection_wasserstein[,2])



g1 <- ggplot(data.frame(posterior_sample), aes(x = X1)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
g1 <- g1 + geom_density(data=wsmc.summary.df, aes(y = ..density.., fill = "Summary", colour = "Summary"), alpha = 0.5)
g1 <- g1 + geom_density(data=wsmc.wasserstein.df, aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g1 <- g1 + geom_density(data=rej.summary.df, aes(y = ..density.., fill = "Rej. Summary", colour = "Rej. Summary"), alpha = 0.5)
#g1 <- g1 + geom_density(data=rej.wasserstein.df, aes(y = ..density.., fill = "Rej. Summary", colour = "Rej. Summary"), alpha = 0.5)
g1 <- g1 + geom_density(data=wsmc.euclidean.df, aes(y = ..density.., fill = "Euclidean", colour = "Euclidean"), alpha = 0.5)
g1 <- g1 + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors)
g1 <- g1 + xlab(expression(theta[1]))
g1 <- g1 + geom_label(data = data.frame(x = c(-0.95,-1, -1.05, -0.31, -1.12), y = c(3.9, 3.4, 2.7, 1.25, 1.7),
                                  method = c("Posterior", "Summary", "Wasserstein", "Euclidean", "Rej. Summary")),
                aes(x = x, y = y, colour = method, label = method), size = 8) + theme(legend.position = "none")
g1
ggsave(filename = paste0(prefix, "mvnorm.posterior1.pdf"), plot = g1, width = 7.5, height = 7)
ggsave(filename = paste0(prefix, "mvnorm.posterior1.png"), plot = g1, width = 7.5, height = 7, dpi = 150)


g2 <- ggplot(data.frame(posterior_sample), aes(x = X2)) + geom_density(aes(y = ..density.., fill = "Posterior", colour = "Posterior"), alpha = 0.5)
g2 <- g2 + geom_density(data=wsmc.summary.df, aes(y = ..density.., fill = "Summary", colour = "Summary"), alpha = 0.5)
g2 <- g2 + geom_density(data=wsmc.wasserstein.df, aes(y = ..density.., fill = "Wasserstein", colour = "Wasserstein"), alpha = 0.5)
g2 <- g2 + geom_density(data=rej.summary.df, aes(y = ..density.., fill = "Rej. Summary", colour = "Rej. Summary"), alpha = 0.5)
g2 <- g2 + geom_density(data=wsmc.euclidean.df, aes(y = ..density.., fill = "Euclidean", colour = "Euclidean"), alpha = 0.5)
g2 <- g2 + scale_color_manual(name = "", values = my_colors) + scale_fill_manual(name = "", values = my_colors)
g2 <- g2 + geom_label(data = data.frame(x = c(-0.13,-0.15, -0.23, 0.55, 0.5), y = c(4, 3.6, 2.8, 1.1, 1.6),
                                        method = c("Summary", "Posterior", "Wasserstein", "Euclidean","Rej. Summary")),
                      aes(x = x, y = y, colour = method, label = method), size = 8) + theme(legend.position = "none")
g2 <- g2 + xlab(expression(theta[2]))
g2
ggsave(filename = paste0(prefix, "mvnorm.posterior2.pdf"), plot = g2, width = 7.5, height = 7)
ggsave(filename = paste0(prefix, "mvnorm.posterior2.png"), plot = g2, width = 7.5, height = 7, dpi = 150)


library(gridExtra)
grid.arrange(g1 + theme(legend.position = "none"), g2 + theme(legend.position = "none"), nrow = 1)


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
    x <- try(exact_transport_distance(t(thetas), t(fast_rmvnorm(nparticles, mu_posterior, sigma_posterior)), p = 1, ground_p = 2))
    if (inherits(x, "try-error")){
      x <- sinkhorn_distance(t(thetas), t(fast_rmvnorm(nparticles, mu_posterior, sigma_posterior)), p = 1, ground_p = 2, eps = 0.1, niterations = 1000)$uncorrected
    }
    x
  })
  return(data.frame(steps = steps, ncomputed = cumsum(results$ncomputed)[steps], times = results$compute_times[steps], w = w_to_post_))
}

filename <- paste0(prefix, "mvnormalwsmc.d", d, ".n", nobservations, ".w1posterior.RData")
n_w_comp <- 50
w_to_post_wasserstein <- wasserstein_to_posterior(results_wasserstein, posterior_sample, n_w_comp)
w_to_post_euclidean <- wasserstein_to_posterior(results_euclidean, posterior_sample, n_w_comp)
w_to_post_summary <- wasserstein_to_posterior(results_summary, posterior_sample, n_w_comp)
w_to_post_rejection_summary <- exact_transport_distance(t(results_rejection_summary), t(posterior_sample), p = 1, ground_p = 2)
save(w_to_post_wasserstein, w_to_post_euclidean, w_to_post_summary, w_to_post_rejection_summary, file = filename)
load(file = filename)
#
# w_to_post_wasserstein %>% head
#
g <- ggplot(w_to_post_wasserstein, aes(x = ncomputed, y = w)) + geom_line(aes(colour = "Wasserstein")) + geom_point(aes(colour = "Wasserstein"))
g <- g + geom_line(data = w_to_post_summary, aes(colour = "Summary")) + geom_point(data = w_to_post_summary, aes(colour = "Summary"))
g <- g + geom_line(data = w_to_post_euclidean, aes(colour = "Euclidean"))  + geom_point(data = w_to_post_euclidean, aes(colour = "Euclidean"))
g <- g + scale_color_manual(name = "", values = my_colors) + xlab("# model simulations") + ylab("W-distance to posterior")
g <- g + scale_x_log10(breaks = c(1e3, 1e4, 1e5, 1e6)) + scale_y_log10()
# g <- g + geom_hline(yintercept = w_to_post_rejection_summary, linetype = 2, aes(colour = "Swap"))
# g <- g + geom_hline(yintercept = self.distances.qt[1], linetype = 2) + geom_hline(yintercept = self.distances.qt[2], linetype = 2)
g <- g + geom_label(data = data.frame(x = c(5e5,5e5,5e5), y = c(0.05, 0.01, 0.5),
                                        method = c("Wasserstein", "Summary", "Euclidean")),
                      aes(x = x, y = y, colour = method, label = method), size = 8) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "mvnorm.w1post.distances.pdf"), plot = g, width = 7.5, height = 7)
ggsave(filename = paste0(prefix, "mvnorm.w1post.distances.png"), plot = g, width = 7.5, height = 7, dpi = 150)
#
#
# g <- ggplot(w_to_post_wasserstein, aes(x = times, y = w)) + geom_line(aes(colour = "Wasserstein")) + geom_point(aes(colour = "Wasserstein"))
# g <- g + geom_line(data = w_to_post_summary, aes(colour = "Summary")) + geom_point(data = w_to_post_summary, aes(colour = "Summary"))
# g <- g + geom_line(data = w_to_post_euclidean, aes(colour = "Euclidean"))  + geom_point(data = w_to_post_euclidean, aes(colour = "Euclidean"))
# g <- g + scale_color_manual(name = "", values = my_colors) + scale_y_log10() + xlab("time (seconds)") + ylab("W-distance to posterior")
# g <- g + scale_x_log10(breaks = c(10, 100, 1000))
# # g <- g + geom_hline(yintercept = self.distances.qt[1], linetype = 2) + geom_hline(yintercept = self.distances.qt[2], linetype = 2)
# g
# # ggsave(filename = paste0(prefix, "mvnorm.w1post.times.pdf"), plot = g, width = 9, height = 5)
