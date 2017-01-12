# load package
library(winference)
# register parallel cores
registerDoMC(cores = 8)
# remove all
rm(list = ls())
# apply preferences for ggplotting
setmytheme()
# set RNG seed
set.seed(11)
#
# This script applies an SMC sampler to sample
# from the Wasserstein quasi-posterior distribution,
# in an auto-regressive model with delay reconstruction.
# We first use the Hilbert distance, and later the regularized Wasserstein
# calculated with Sinkhorn's algorithm.

# The model is autoregressive of order 1.
# The model specifies Y_1 ~ Normal(0, sigma^2/(1-rho^2)),
# and Y_t ~ Normal(rho Y_{t-1}, sigma^2) for all t >= 2.
# The parameters are (rho, sigma).
# The prior is Uniform on [-1,1] on rho,
# and Normal(0,1) on log(sigma).

# load autoregressive model
target <- get_autoregressive()
# number of observations
nobservations <- 100
# data-generating parameter (parametrized as rho and log(sigma))
true_theta <- c(0.7, 0.9)
# observations generated from model
obs <- target$robservation(nobservations, true_theta,
                           target$parameters, target$generate_randomness(nobservations))
# plot data
plot(obs, type = "l")

# now we introduce delay reconstructions, with a lag of 1
lagvalue <- 1
# we construct a matrix of (y_t, y_t-1) for t >= 2 (t indexes the columns)
lag_obs <- create_lagmatrix(matrix(obs, nrow = 1), lagvalue)
lag_obs <- lag_obs[,(lagvalue+1):ncol(lag_obs)]
# we can sort this matrix according to the Hilbert space filling curve as follows:
order_obs <- hilbert_order(lag_obs)
orderded_obs <- lag_obs[,order_obs]
#
# compute Hilbert pseudo-distance (of order 1) between data generated given theta and the actual data
compute_hilbert_distance <- function(theta){
  # generate fake data
  fake_rand <- target$generate_randomness(nobservations)
  fake_obs <- target$robservation(nobservations, theta, target$parameters, fake_rand)
  # create delay reconstruction
  fake_obs <- create_lagmatrix(matrix(fake_obs, nrow = 1), lagvalue)
  fake_obs <- fake_obs[,(lagvalue+1):ncol(fake_obs)]
  # Hilbert sort
  order_fake <- hilbert_order(fake_obs)
  # compute distance of sorted samples
  distance <- nrow(fake_obs) * mean(abs(orderded_obs - fake_obs[,order_fake]))
  return(distance)
}

# let's try to compute the Hilbert distance for draws with rho = 0.5, log(sigma) = 1
compute_hilbert_distance(c(0.5,1))
compute_hilbert_distance(c(0.5,1))
compute_hilbert_distance(c(0.5,1))
# each time we call the distance function, we get a different value, due to the randomness
# in the sampling of the synthetic dataset

# we now call the SMC sampler to explore the quasi-posterior as the tolerance/threshold
# epsilon goes to zero

# we use a mixture of Normal distributions for the proposal
proposal <- mixture_proposal()
# we specify algorithmic parameters,
# such as the number of particles, the number of moves per rejuvenation step,
# the proposal mechanism, the number of steps to perform, the diversity parameter
# used in the threshold adaptation, the number of hits to use in the r-hit kernel,
# and the maximum number of trials to use in the r-hit kernel before rejecting.
param_algo <- list(nthetas = 1024, nmoves = 1, proposal = proposal,
                   nsteps = 5, minimum_diversity = 0.5, R = 2, maxtrials = 1e5)
# we run the algorithm; on my laptop it takes about 20 seconds
wsmcresults <- wsmc_rhit(compute_hilbert_distance, target, param_algo)
# we have access to all the generated particles, distances, and thresholds
names(wsmcresults)
# we transform the output for easy plotting
wsmc.df <- wsmc_to_dataframe(wsmcresults, target$parameter_names)
nsteps <- max(wsmc.df$step)
# let's see the evolution of the thresholds
gthreshold <- qplot(x = 1:(length(wsmcresults$threshold_history)), y = wsmcresults$threshold_history, geom = "line")
gthreshold <- gthreshold + xlab("step") + ylab("threshold")
gthreshold
# the thresholds are still going down, so we should do more steps

# let's plot rho and log(sigma) on a scatter plot, colored by step
g <- ggplot(wsmc.df, aes(x = rho, y = logsigma, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
g <- g + xlab(expression(rho)) + ylab(expression(log(sigma)))
g

# we can also look at the marginals
g <- ggplot(wsmc.df, aes(x = rho, colour = step, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(rho))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue") + geom_vline(xintercept = true_theta[1])
g
#
g <- ggplot(wsmc.df, aes(x = logsigma, colour = step, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(log(sigma)))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue") + geom_vline(xintercept = true_theta[2])
g

####
# now we can do the same experiments, but calculating the distance
# with regularized Wasserstein techniques instead of Hilbert

# compute regularized transport distance between delay reconstructions
compute_wasserstein_distance <- function(theta, transportiterations = 100){
  # generate fake data
  fake_rand <- target$generate_randomness(nobservations)
  fake_obs <- target$robservation(nobservations, theta, target$parameters, fake_rand)
  # create delay reconstruction
  lag_fake_obs <- create_lagmatrix(matrix(fake_obs, nrow = target$ydim), lagvalue)
  # compute matrix of pairwise distances
  C <- cost_matrix_L1(lag_obs, lag_fake_obs[,(lagvalue+1):nobservations,drop=FALSE])
  # regularization parameter
  epsilon <- 0.05 * median(C)
  # dummy vector of equal weights attributed to the empirical distributions
  equalw <- rep(1/(nobservations-lagvalue), (nobservations-lagvalue))
  # call Sinkhorn's algorithm
  wass <- wasserstein(equalw, equalw, C, epsilon, transportiterations)
  return(as.numeric(wass$distances))
}

# let's try to compute the Wasserstein distance for draws with rho = 0.5, log(sigma) = 1
compute_wasserstein_distance(c(0.5,1))
compute_wasserstein_distance(c(0.5,1))
compute_wasserstein_distance(c(0.5,1))

# we now call the SMC sampler to explore the quasi-posterior as the tolerance/threshold
# epsilon goes to zero; this takes longer than with the Hilbert curve, especially
# for large data sets
wsmcresults <- wsmc_rhit(compute_wasserstein_distance, target, param_algo)
wsmc.df <- wsmc_to_dataframe(wsmcresults, target$parameter_names)
nsteps <- max(wsmc.df$step)

# let's plot rho and log(sigma) on a scatter plot, colored by step
g <- ggplot(wsmc.df, aes(x = rho, y = logsigma, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
g <- g + xlab(expression(rho)) + ylab(expression(log(sigma)))
g

# we can also look at the marginals
g <- ggplot(wsmc.df, aes(x = rho, colour = step, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(rho))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue") + geom_vline(xintercept = true_theta[1])
g
#
g <- ggplot(wsmc.df, aes(x = logsigma, colour = step, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(log(sigma)))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue") + geom_vline(xintercept = true_theta[2])
g
