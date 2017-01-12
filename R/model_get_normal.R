# Normal Gamma model
# theta = (mu, tau)
#'@rdname get_normal
#'@title Normal model
#'@description This function returns a list representing
#' a Normal location model.
#' The prior is mu ~ Normal(mu_0, nu^{-1}), where nu is precision
#' and  tau ~ Gamma(alpha, beta), where beta is rate (1/scale).
#'  The likelihood is Y ~ Normal(mu, tau^2) where tau^2 is the variance.
#'@return The list contains rprior, dprior (generate and evaluate the density of prior distribution),
#' generate_randomness (generate data-generating variables), robservation (create synthetic
#' data sets), parameter_names (useful for plotting), thetadim (dimension of parameter),
#' ydim (dimension of observations), parameters (list of hyperparameters,
#' to be passed to rprior,dprior,robservation)
#'@export
get_normal <- function(){
  rprior <- function(nparticles, parameters){
    particles <- matrix(nrow = nparticles, ncol = 2)
    particles[,1] <- rnorm(nparticles, mean = parameters$mu_0, sd = 1/sqrt(parameters$nu))
    particles[,2] <- rgamma(nparticles, shape = parameters$alpha, rate = parameters$beta)
    return(particles)
  }
  # evaluate the log-density of the prior, for each particle
  # parameters is a list containing mu_0, nu, alpha, beta
  dprior <- function(thetaparticles, parameters){
    logdensities <- dnorm(thetaparticles[,1], mean = parameters$mu_0, sd = 1/sqrt(parameters$nu), log = TRUE)
    logdensities <- logdensities + dgamma(thetaparticles[,2], shape = parameters$alpha, rate = parameters$beta, log = TRUE)
    return(logdensities)
  }
  # generate random variables used to compute a synthetic dataset
  generate_randomness <- function(nobservations){
    return(rnorm(nobservations))
  }
  # function to compute a dataset for each theta value
  robservation <- function(nobservations, theta, parameters, randomness){
    observations <- theta[1] + randomness * theta[2]
    return(observations)
  }
  parameters <- list(mu_0 = 0, nu = 1, alpha = 2, beta = 1)
  #
  model <- list(rprior = rprior,
                dprior = dprior,
                generate_randomness = generate_randomness,
                robservation = robservation,
                parameter_names = c("mu", "sigma"),
                parameters = parameters,
                thetadim = 2, ydim = 1)
  return(model)
}
