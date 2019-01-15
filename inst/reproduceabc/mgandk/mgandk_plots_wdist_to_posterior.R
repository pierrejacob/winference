library(winference)
registerDoParallel(cores = detectCores())
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


param_algo <- list(nthetas = 2048, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 100000)


filename <- paste0(prefix, "mgandk.wsmc.n", nobservations, ".hilbert.RData")
load(filename)
results_hilbert <- results
wsmc.hilbert.df <- wsmc_to_dataframe(results_hilbert)
sum(results_hilbert$ncomputed)
step_hilbert <- length(results_hilbert$thetas_history)

#
filename <- paste0(prefix, "mgandk.wsmc.n", nobservations, ".mmd.RData")
load(filename)
results_mmd <- results
wsmc.mmd.df <- wsmc_to_dataframe(results_mmd)
sum(results_mmd$ncomputed)
step_mmd <- length(results_mmd$thetas_history)

#
filename <- paste0(prefix, "mgandk.wsmc.n", nobservations, ".swap.RData")
load(filename)
results_swap <- results
wsmc.swap.df <- wsmc_to_dataframe(results_swap)
sum(results_swap$ncomputed)
step_swap <- length(results_swap$thetas_history)

#
filename <- paste0(prefix, "mgandk.wsmc.n", nobservations, ".wasserstein.RData")
load(filename)
results_wasserstein <- results
wsmc.wasserstein.df <- wsmc_to_dataframe(results_wasserstein)
sum(results_wasserstein$ncomputed)
step_wasserstein <- length(results_wasserstein$thetas_history)

#
mhfile <- paste0(prefix, "mgandk.mcmc.n", nobservations, ".mh.initfromABC.RData")
load(mhfile)
if (exists("mh_results")){
  mcmc.df <- mhchainlist_to_dataframe(mh_results$chains)
  mcmc.df <- mcmc.df %>% filter(iteration < mh_results$iteration)
} else {
  mcmc.df <- mhchainlist_to_dataframe(mh$chains)
}
burnin <- 50000
names(mcmc.df) <- c("ichain", "iteration", target$parameter_names)
target$parameter_names


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

npost <- 2048
mcmcpostburnin <- mcmc.df %>% filter(iteration > burnin)
posterior_sample <- mcmcpostburnin[sample(x = 1:nrow(mcmcpostburnin), size = npost, replace = TRUE),3:ncol(mcmcpostburnin)]
posterior_sample2 = mcmcpostburnin[sample(x = 1:nrow(mcmcpostburnin), size = npost, replace = TRUE),3:ncol(mcmcpostburnin)]

filename <- paste0(prefix, "mgandk.n", nobservations, ".w1posterior.RData")
n_w_comp <- 50
w_to_post_wasserstein <- wasserstein_to_posterior(results_wasserstein, posterior_sample, n_w_comp)
w_to_post_hilbert <- wasserstein_to_posterior(results_hilbert, posterior_sample, n_w_comp)
w_to_post_swap <- wasserstein_to_posterior(results_swap, posterior_sample, n_w_comp)
w_to_post_mmd <- wasserstein_to_posterior(results_mmd, posterior_sample, n_w_comp)
save(w_to_post_wasserstein, w_to_post_hilbert, w_to_post_swap, w_to_post_mmd, file = filename)
load(file = filename)

w_to_post_wasserstein %>% tail

g <- ggplot(w_to_post_wasserstein, aes(x = ncomputed, y = w)) + geom_line(aes(colour = "Wasserstein")) + geom_point(aes(colour = "Wasserstein"))
g <- g + geom_line(data = w_to_post_hilbert, aes(colour = "Hilbert")) + geom_point(data = w_to_post_hilbert, aes(colour = "Hilbert"))
g <- g + geom_line(data = w_to_post_swap, aes(colour = "Swap")) + geom_point(data = w_to_post_swap, aes(colour = "Swap"))
g <- g + geom_line(data = w_to_post_mmd, aes(colour = "MMD")) + geom_point(data = w_to_post_mmd, aes(colour = "MMD"))
g <- g + scale_color_manual(name = "", values = my_colors) + xlab("# model simulations") + ylab("W-distance to posterior")
g <- g + scale_x_log10(breaks = c(5*1e4,1e5,2*1e5,5*1e5,1e6,2*1e6), limits = c(5*1e4,3*1e6)) + scale_y_log10(breaks = c(3,4,5,7,10))
g <- g + geom_label(data = data.frame(x = c(5e5, 2.85e6, 2.8e6, 2.7e6), y = c(3.1, 3.5, 3, 4.9), method = c("Wasserstein", "Hilbert", "Swap", "MMD")),
                    aes(x = x, y = y, colour = method, label = method), size = 7) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "mgandk.n", nobservations, ".w1post.pdf"), plot = g, width = 2*fig.width, height = fig.height)
#ggsave(filename = paste0(prefix, "mgandk.n", nobservations, ".w1post.png"), plot = g, width = 2*fig.width, height = fig.height, dpi = 150)
