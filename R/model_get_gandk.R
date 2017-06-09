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
  loglikelihood <- function(thetas, ys, ...){
    n <- length(ys)
    evals <- rep(0, nrow(thetas))
    for (itheta in 1:nrow(thetas)){
      ll <- function(ys, h = 1e-5, tolerance = 1e-10){
        all_ys <- c(ys-h, ys+h)
        o <- order(all_ys)
        x <- rep(0, length(all_ys))
        x[o[1]] <- gandkcdf(y = all_ys[o[1]], theta = thetas[itheta,], tolerance = tolerance)
        for (i in 2:length(all_ys)){
          x[o[i]] <- gandkcdf(y = all_ys[o[i]], theta = thetas[itheta,], tolerance = tolerance, lower = x[o[i-1]])
        }
        return(sum(log((x[(n+1):(2*n)] - x[1:n])/(2*h))))
      }
      evals[itheta] <- ll(ys)
    }
    return(evals)
  }
  parameters <- list()
  #
  model <- list(rprior = rprior,
                dprior = dprior,
                generate_randomness = generate_randomness,
                robservation = robservation,
                loglikelihood = loglikelihood,
                parameter_names = c("A", "B", "g", "k"),
                parameters = parameters,
                thetadim = 4, ydim = 1)
  return(model)
}
