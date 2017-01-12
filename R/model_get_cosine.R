# cosine model
# theta = (phi, logsigma)
#'@rdname get_cosine
#'@title Cosine model
#'@description This function returns a list representing a cosine trend model.
#'@return The list contains rprior, dprior (generate and evaluate the density of prior distribution),
#' generate_randomness (generate data-generating variables), robservation (create synthetic
#' data sets), parameter_names (useful for plotting), thetadim (dimension of parameter),
#' ydim (dimension of observations), parameters (list of hyperparameters,
#' to be passed to rprior,dprior,robservation)
#'@export
get_cosine <- function(){
  # generate phi ~ Unif(0,2pi) and log-sigma is normally distributed
  rprior <- function(nparticles, parameters){
    omegas <- runif(nparticles, min = 0, max = 1/10)
    phis <- runif(nparticles, min = 0, max = 2*pi)
    logsigma <- rnorm(nparticles)
    logA <- rnorm(nparticles)
    return(cbind(omegas, phis, logsigma, logA))
  }

  # evaluate the log-density of the prior, for each particle
  dprior <- function(thetas, parameters){
    logdensities <- dnorm(thetas[,3], 0, 1, log = TRUE)
    logdensities <- logdensities + dnorm(thetas[,4], 0, 1, log = TRUE)
    logdensities[thetas[,1] > 1/10] <- -Inf
    logdensities[thetas[,1] < 0] <- -Inf
    logdensities[thetas[,2] > 2*pi] <- -Inf
    logdensities[thetas[,2] < 0] <- -Inf
    return(logdensities)
  }
  #
  generate_randomness <- function(nobservations){
    return(rnorm(nobservations))
  }
  # function to generate a dataset for each theta value
  robservation <- function(nobservations, theta, parameters, randomness){
    observations <- exp(theta[4]) * cos(2*pi*theta[1]*(1:nobservations) + theta[2]) + exp(theta[3]) * randomness
    return(observations)
  }
  #
  model <- list(rprior = rprior,
                dprior = dprior,
                generate_randomness = generate_randomness,
                robservation = robservation,
                parameter_names = c("omega", "phi", "logsigma", "logA"),
                thetadim = 4, ydim = 1,
                parameters = list())
  return(model)
}
