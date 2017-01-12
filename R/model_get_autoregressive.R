# Autoregressive model
# theta = (phi, logsigma)
#'@rdname get_autoregressive
#'@title Autoregressive model
#'@description This function returns a list representing an auto-regressive
#'model of order 1.
#'@return The list contains rprior, dprior (generate and evaluate the density of prior distribution),
#' generate_randomness (generate data-generating variables), robservation (create synthetic
#' data sets), parameter_names (useful for plotting), thetadim (dimension of parameter),
#' ydim (dimension of observations), parameters (list of hyperparameters,
#' to be passed to rprior,dprior,robservation)
#'@export
get_autoregressive <- function(){
  # generate phi ~ Unif(-1,1) and log-sigma is normally distributed
  rprior <- function(nparticles, parameters){
    phis <- 2*runif(nparticles) - 1
    logsigmas <- rnorm(nparticles, 0, 1)
    return(cbind(phis, logsigmas))
  }

  # evaluate the log-density of the prior, for each particle
  dprior <- function(thetas, parameters){
    logdensities <- dnorm(thetas[,2], 0, 1, log = TRUE)
    logdensities[thetas[,1] > 1] <- -Inf
    logdensities[thetas[,1] < -1] <- -Inf
    return(logdensities)
  }
  #
  generate_randomness <- function(nobservations){
    return(rnorm(nobservations))
  }
  # function to generate a dataset for each theta value
  # X_0 is Normal(0,sigma^2/(1-phi^2))
  # X_t is Normal(phi X_t-1,sigma^2)
  robservation <- function(nobservations, theta, parameters, randomness){
    observations <- rep(0, nobservations)
    observations[1] <- randomness[1] * exp(theta[2]) / sqrt(1 - theta[1]^2)
    for (idata in 2:nobservations){
      observations[idata] <- theta[1] * observations[idata-1] + randomness[idata] * exp(theta[2])
    }
    return(observations)
  }
  #
  model <- list(rprior = rprior,
                dprior = dprior,
                generate_randomness = generate_randomness,
                robservation = robservation,
                parameter_names = c("rho", "logsigma"),
                thetadim = 2, ydim = 1,
                parameters = list())
  return(model)
}
