library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()
set.seed(11)
target <- get_gandk()

fig.height <- 5
fig.width <- 5

my_colors <- get_my_colors()

prefix = ""

nobservations = 250
load(paste0(prefix, "gandkdata.RData"))
obs <- obs[1:nobservations]

# Wasserstein SMC
filename = paste0(prefix, "gandkwsmc.n", nobservations, ".RData")
load(filename)
results_was = results


# MCMC
mhfile <- paste0(prefix, "gandkmcmc.n", nobservations, "mh.RData")
load(mhfile)
mcmc.df <- mhchainlist_to_dataframe(mh$chains)
names(mcmc.df) <- c("ichain", "iteration", target$parameter_names)
burnin = 50000
mcmc.df = mcmc.df[mcmc.df$iteration>burnin,]


#Calculate Wasserstein of WABC to posterior
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

filename <- paste0(prefix, "gandk.n", nobservations, ".w1posterior.RData")
n_w_comp <- 50
w_to_post_wasserstein <- wasserstein_to_posterior(results_was, posterior_sample, n_w_comp)
save(w_to_post_wasserstein, file = filename)
load(file = filename)

g <- ggplot(w_to_post_wasserstein, aes(x = ncomputed, y = w)) + geom_line(aes(colour = "Wasserstein")) + geom_point(aes(colour = "Wasserstein"))
g <- g + scale_color_manual(name = "", values = my_colors) + xlab("# model simulations") + ylab("W-distance to posterior")
g <- g + scale_x_log10(breaks = c(1e4,1e5,1e6,1e7,1e8)) + scale_y_log10(breaks = c(10,5,2,1,0.5,0.25,0.1,0.06))
g <- g + geom_label(data = data.frame(x = c(8e6), y = c(0.1), method = c("Wasserstein")),
                    aes(x = x, y = y, colour = method, label = method), size = 7) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "gandk.n", nobservations, ".w1post.pdf"), plot = g, width = 2*fig.width, height = fig.height)



