library(winference)
registerDoParallel(cores = 6)
rm(list = ls())
setmytheme()
set.seed(11)
# model
target <- get_autoregressive()
#
# number of observations
nobservations <- 1000
prefix <- ""
load(file = paste0(prefix, "ar1data.RData"))
obs <- obs[1:nobservations]

## distribution obtained with Wasserstein on the marginal
filename <- paste0(prefix, "ar1.n", nobservations, ".wsmc_marginal.RData")
load(filename)
results_marginal <- results
wsmc.df <- wsmc_to_dataframe(results_marginal)
nsteps_marginal <- length(results_marginal$thetas_history)
# evolution of distance threshold
plot_threshold_time(results_marginal) + scale_y_log10()
plot_ncomputed(results_marginal)
# distribution of distances
qplot(x = results_marginal$distances_history[[nsteps_marginal]], geom = "histogram") + xlab("distances")
# bivariate ABC posterior
plot_bivariate(results_marginal, 1, 2)
g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = rho, y = logsigma, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps_marginal/2)) + theme(legend.position = "none")
g <- g + xlab(expression(rho)) + ylab(expression(log(sigma)))
g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
print(g)

ggsave(filename = paste0(prefix, "ar1_identification_withoutlag.png"), plot = g, width = 7, height = 5, dpi = 150)

## distribution obtained with Wasserstein on the bivariate marginal of (y_t,y_t+1)
## subsetted to have nobservations / 2 points
filename <- paste0(prefix, "ar1.n", nobservations, ".wsmc_delay1.RData")
load(filename)
results_delay1 <- results
wsmc.df <- wsmc_to_dataframe(results_delay1)
nsteps_delay1 <- length(results_delay1$thetas_history)
# evolution of distance threshold
plot_threshold_time(results_delay1) + scale_y_log10() + geom_point()
# evolution of cost per step
plot_ncomputed(results_delay1) + geom_point()
# distribution of distances
qplot(x = results_delay1$distances_history[[nsteps_delay1]], geom = "histogram") + xlab("distances")
# bivariate ABC posterior
# plot_bivariate(results_delay1, 1, 2)
wsmc.df <- wsmc_to_dataframe(results_delay1)
g <- ggplot(wsmc.df %>% filter(step > 0), aes(x = rho, y = logsigma, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps_marginal/2)) + theme(legend.position = "none")
g <- g + xlab(expression(rho)) + ylab(expression(log(sigma)))
g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
print(g)

ggsave(filename = paste0(prefix, "ar1_identification_withlag1.png"), plot = g, width = 7, height = 5, dpi = 150)

wabc_thetas <- results_delay1$thetas_history[[nsteps_delay1]]
## compare ABC posterior with exact posterior, and posterior in alternate, bivariate Normal model
filename <- paste0(prefix, "ar1.n", nobservations, ".mcmc.RData")
load(filename)
mh_posterior <- mh
chains.df <- mhchainlist_to_dataframe(mh_posterior$chains)
niterations <- max(chains.df$iteration)
# ggplot(chains.df, aes(x = iteration, y = X.1, group = ichain)) + geom_line()
burnin <- niterations / 2
mcmc_thetas <- as.matrix((chains.df %>% filter(iteration > burnin))[,3:4])
#
filename <- paste0(prefix, "ar1.n", nobservations, ".mcmc_alternate.RData")
load(filename)
mh_alternate <- mh
chains_alternate.df <- mhchainlist_to_dataframe(mh_alternate$chains)
niterations <- max(chains_alternate.df$iteration)
burnin <- niterations / 2
mcmc_alternate_thetas <- as.matrix((chains_alternate.df %>% filter(iteration > burnin, iteration %% 10 == 1))[,3:4])

grho <- qplot(x = mcmc_thetas[,1], geom = "blank") + geom_histogram(aes(y = ..density..))
grho <- grho + geom_histogram(aes(x = mcmc_alternate_thetas[,1], y = ..density..), fill = "red", alpha = 0.4)
grho <- grho + geom_histogram(aes(x = wabc_thetas[,1], y = ..density..), fill = "blue", alpha = 0.4)
grho

gsigma <- qplot(x = mcmc_thetas[,2], geom = "blank") + geom_histogram(aes(y = ..density..))
gsigma <- gsigma + geom_histogram(aes(x = mcmc_alternate_thetas[,2], y = ..density..), fill = "red", alpha = 0.4)
gsigma <- gsigma + geom_histogram(aes(x = wabc_thetas[,2], y = ..density..), fill = "blue", alpha = 0.4)
gsigma

