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
# in a Normal location model.

# The parameters are (mu, sigma).
# The prior is standard Normal on mu,
# and Gamma(2,1) on sigma.

# The model specifies Y ~ Normal(mu, sigma).
# We are going to fit it on a sample obtained for a Gamma distribution,
# so we are dealing with a misspecified setting.

# number of observations
nobservations <- 1000
# observations from a Gamma model, with mean 2
obs <- rgamma(nobservations, shape = 10, rate = 5)
# sort the observations, to facilitate later computations
obs_sorted <- sort(obs)

# function to generate from prior distribution
# first argument is number of desired samples
# second argument contains hyper-parameters
rprior <- function(nparticles, parameters){
  particles <- matrix(nrow = nparticles, ncol = 2)
  particles[,1] <- rnorm(nparticles, mean = parameters$mu_0, sd = 1/sqrt(parameters$nu))
  particles[,2] <- rgamma(nparticles, shape = parameters$alpha, rate = parameters$beta)
  return(particles)
}
# function to evaluate prior log-density
# first argument is a matrix of parameters (one per row)
# second argument contains hyper-parameters
dprior <- function(thetaparticles, parameters){
  logdensities <- dnorm(thetaparticles[,1], mean = parameters$mu_0, sd = 1/sqrt(parameters$nu), log = TRUE)
  logdensities <- logdensities + dgamma(thetaparticles[,2], shape = parameters$alpha, rate = parameters$beta, log = TRUE)
  return(logdensities)
}

# we collect the prior and other aspects of the model
# in a list; the "parameters" list contains the hyper-parameters
target <- list(rprior = rprior,
              dprior = dprior,
              parameter_names = c("mu", "sigma"),
              parameters = list(mu_0 = 0, nu = 1, alpha = 2, beta = 1),
              thetadim = 2, ydim = 1)

# we specify a data-generating mechanism, given a parameter theta
robservation <- function(nobservations, theta){
  observations <- theta[1] + rnorm(nobservations) * theta[2]
  return(observations)
}

# function to compute 1-Wasserstein distance between observed data and data generated given theta
compute_distance <- function(theta){
  fake_obs <- robservation(nobservations, theta)
  fake_obs_sorted <- sort(fake_obs)
  return(sum(abs(obs_sorted - fake_obs_sorted))/length(obs_sorted))
}

# let's try to compute the Wasserstein distance for a draw with mu = 0, sigma = 1
compute_distance(c(0,1))
compute_distance(c(0,1))
compute_distance(c(0,1))
# each time we call the distance function, we get a different value, due to the randomness
# in the sampling of the synthetic dataset

# We now call the SMC sampler to explore the quasi-posterior as the tolerance/threshold
# epsilon goes to zero

# we could use a Normal proposal fitted on the data, by specifying the proposal as
proposal <- independent_proposal()
# or we could use a mixture of Normal distributions, with the following line
proposal <- mixture_proposal()

# we specify algorithmic parameters,
# such as the number of particles, the number of moves per rejuvenation step,
# the proposal mechanism, the number of steps to perform, the diversity parameter
# used in the threshold adaptation, the number of hits to use in the r-hit kernel,
# and the maximum number of trials to use in the r-hit kernel before rejecting.
param_algo <- list(nthetas = 1024, nmoves = 1, proposal = proposal,
                   nsteps = 5, minimum_diversity = 0.5, R = 2, maxtrials = 1e5)
# we run the algorithm; on my laptop it takes about 20 seconds
wsmcresults <- wsmc_rhit(compute_distance, target, param_algo)
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

# let's plot mu and sigma on a scatter plot, colored by step
g <- ggplot(wsmc.df, aes(x = mu, y = sigma, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + geom_vline(xintercept = 10/5) + geom_hline(yintercept = sqrt(10/25))
g <- g + xlab(expression(mu)) + ylab(expression(sigma))
g

# we can also look at the marginals
g <- ggplot(wsmc.df, aes(x = mu, colour = step, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(mu))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue") + geom_vline(xintercept = 10/5)
g
#
g <- ggplot(wsmc.df, aes(x = sigma, colour = step, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(sigma))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue") + geom_vline(xintercept = sqrt(10/25))
g
