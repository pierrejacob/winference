#'@rdname get_gandk
#'@title G and k model
#'@description This function returns a list representing the g-and-k
#' quantile distribution.
#'@return The list contains rprior, dprior (generate and evaluate the density of prior distribution),
#' generate_randomness (generate data-generating variables), robservation (create synthetic
#' data sets), parameter_names (useful for plotting), thetadim (dimension of parameter),
#' ydim (dimension of observations), parameters (list of hyperparameters,
#' to be passed to rprior,dprior,robservation)
#'@export
get_gandk <- function(){
  rprior <- function(nparticles, parameters){
    return(matrix(runif(nparticles*4, min = 0, max = 10), ncol = 4))
  }
  # evaluate the log-density of the prior, for each particle
  dprior <- function(thetaparticles, parameters){
    densities <- rep(0, nrow(thetaparticles))
    for (i in 1:nrow(thetaparticles)){
      if (any(thetaparticles[i,] > 10) || any(thetaparticles[i,] < 0)){
        densities[i] <- -Inf
      }
    }
    return(densities)
  }
  # generate random variables used to compute a synthetic dataset
  generate_randomness <- function(nobservations){
    return(rnorm(nobservations))
  }
  # function to compute a dataset for each theta value
  robservation <- function(nobservations, theta, parameters, randomness){
    observations <- gandkinversecdf_givennormals(randomness, theta)
    return(observations)
  }
  parameters <- list()
  #
  model <- list(rprior = rprior,
                dprior = dprior,
                generate_randomness = generate_randomness,
                robservation = robservation,
                parameter_names = c("A", "B", "g", "k"),
                parameters = parameters,
                thetadim = 4, ydim = 1)
  return(model)
}
