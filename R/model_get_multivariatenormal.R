# The observations are i.i.d. Normal with mean mu (parameters)
# and covariance matrix Sigma, defined
# as the identity matrix, and +0.5 on the upper and lower diagonals.
#'@rdname get_multivariate_normal
#'@title Multivariate Normal model
#'@description This function returns a list representing
#' a Normal location model, in a dimension specified by the user as an argument.
#'@return The list contains rprior, dprior (generate and evaluate the density of prior distribution),
#' generate_randomness (generate data-generating variables), robservation (create synthetic
#' data sets), parameter_names (useful for plotting), thetadim (dimension of parameter),
#' ydim (dimension of observations), parameters (list of hyperparameters,
#' to be passed to rprior,dprior,robservation)
#'@export
get_multivariate_normal <- function(dimension){
  target <- list()
  target$rprior <- function(nparticles, parameters){
    return(fast_rmvnorm(nparticles, rep(parameters$mu_0, dimension), diag(parameters$tau^2, dimension, dimension)))
    # particles <- matrix(nrow = nparticles, ncol = dimension)
    # for (id in 1:dimension){
    #   particles[,id] <- rnorm(nparticles, mean = parameters$mu_0, sd = parameters$tau)
    # }
    # return(particles)
  }
  # evaluate the log-density of the prior, for each particle
  # parameters is a list containing mu_0, nu, alpha, beta
  target$dprior <- function(thetas, parameters){
    # logdensities <- rep(0, nrow(thetas))
    return(fast_dmvnorm(thetas, rep(parameters$mu_0, dimension), diag(parameters$tau^2, dimension, dimension)))
    # for (id in 1:dimension){
      # logdensities <- logdensities + dnorm(thetas[,id],  mean = parameters$mu_0, sd = parameters$tau, log = TRUE)
    # }
    # return(logdensities)
  }

  # generate random variables used to compute a synthetic dataset
  target$generate_randomness <- function(nobservations){
    return(fast_rmvnorm(nobservations, rep(0, dimension), diag(1, dimension, dimension)))
  }

  S <- diag(1, dimension, dimension)
  if (dimension > 1){
    for (i in 1:(dimension-1)){
      S[i,i+1] <- S[i+1,i] <- 0.5
    }
  }
  target$parameters <- list(S = S, mu_0 = 0, tau = 10)
  # function to compute a dataset for each theta value
  target$robservation <- function(nobservations, theta, parameters, ...){
    return(t(fast_rmvnorm(nobservations, theta, parameters$S)))
  }
  target$loglikelihood <- function(thetaparticles, observations, parameters){
    logdensities <- rep(0, nrow(thetaparticles))
    for (i in 1:nrow(thetaparticles)){
      logdensities[i] <- sum(fast_dmvnorm(t(observations), thetaparticles[i,], parameters$S))
    }
    return(logdensities)
  }
  target$thetadim <- dimension
  target$ydim <- dimension
  target$parameter_names <- paste0("X", 1:dimension)
  return(target)
}
