#'@rdname get_ricker
#'@title Ricker model
#'@description This function returns a list representing
#' the Ricker model in
#' Wood (2010) Statistical inference for noisy nonlinear ecological dynamic systems
#'@return The list contains rprior, dprior (generate and evaluate the density of prior distribution),
#' generate_randomness (generate data-generating variables), robservation (create synthetic
#' data sets), parameter_names (useful for plotting), thetadim (dimension of parameter),
#' ydim (dimension of observations), parameters (list of hyperparameters,
#' to be passed to rprior,dprior,robservation)
#'@export
get_ricker <- function(){
  rprior <- function(nparticles, ...){
    theta1 <- runif(nparticles, min = 0, max = 10)
    theta2 <- runif(nparticles, min = 0, max = 20)
    theta3 <- runif(nparticles, min = 0, max = 2)
    return(cbind(theta1, theta2, theta3))
  }
  dprior <- function(thetas, ...){
    if (is.null(dim(thetas))) thetas <- matrix(thetas, nrow = 1)
    density_evals <- dunif(thetas[,1], min = 0, max = 10, log = TRUE)
    density_evals <- density_evals + dunif(thetas[,2], min = 0, max = 20, log = TRUE)
    density_evals <- density_evals + dunif(thetas[,3], min = 0, max = 2, log = TRUE)
    return(density_evals)
  }
  #
  generate_randomness <- function(nobservations){
    return(list())
  }
  robservation <- function(nobservations, theta, parameters, randomness){
    obs <- rep(0, nobservations)
    r <- exp(theta[1])
    phi <- theta[2]
    sigma_e <- theta[3]
    state <- 1
    for (t in 1:nobservations){
      state = r * state * exp(-state + sigma_e * rnorm(1))
      obs[t] <- rpois(1, phi*state)
    }
    return(obs)
  }
  #
  model <- list(rprior = rprior,
                dprior = dprior,
                generate_randomness = generate_randomness,
                robservation = robservation,
                parameters = NULL,
                parameter_names = c("logr", "phi", "sigma_e"),
                thetadim = 3, ydim = 1)
  return(model)
}
