#'@rdname get_pz_6param
#'@title Phytoplankton-zooplankton model
#'@description This function returns a list representing
#' a Lotka-Volterra type model for plankton. See
#' Jones, E., Parslow, J., and Murray, L. (2010). A Bayesian approach to state and parameter estimation in a phytoplankton-zooplankton model. Australian Meteorological and Oceanographic Journal, 59:7–16.
#'@return The list contains rprior, dprior (generate and evaluate the density of prior distribution),
#' generate_randomness (generate data-generating variables), robservation (create synthetic
#' data sets), parameter_names (useful for plotting), thetadim (dimension of parameter),
#' ydim (dimension of observations), parameters (list of hyperparameters,
#' to be passed to rprior,dprior,robservation)
#'@export
get_pz_6param <- function(){
  rprior <- function(nparticles, ...){
    ## evaluate prior density on the transformed parameter
    theta1 <- runif(nparticles)
    theta2 <- runif(nparticles)
    # set the other parameters deterministically
    theta3 <- runif(nparticles)
    theta4 <- runif(nparticles)
    theta5 <- runif(nparticles)
    theta6 <- runif(nparticles)

    return(cbind(theta1, theta2, theta3, theta4, theta5, theta6))
  }
  # prior, on the transformed parameters
  dprior <- function(thetas, ...){
    ## evaluate prior density on the transformed parameter
    if (is.null(dim(thetas))) thetas <- matrix(thetas, nrow = 1)
    density_evals <- dunif(thetas[,1], log = TRUE) + dunif(thetas[,2], log = TRUE)
    density_evals <- density_evals + dunif(thetas[,3], log = TRUE) + dunif(thetas[,4], log = TRUE)
    density_evals <- density_evals + dunif(thetas[,5], log = TRUE) + dunif(thetas[,6], log = TRUE)
    return(density_evals)
  }

  #
  generate_randomness <- function(nobservations){
    return(list(x_0 = rnorm(2), x = rnorm(nobservations), y = rnorm(nobservations)))
  }
  robservation <- function(nobservations, theta, parameters, randomness){
    # untr_theta <- untransform_theta(theta)
    state <- matrix(exp(log(2) + randomness$x_0), nrow = 2)
    states <- matrix(nrow = 2, ncol = nobservations+1)
    states[,1] <- state
    log_obs <- rep(0, nobservations)
    for (t in 1:nobservations){
      alpha <- theta[2] * randomness$x[t] + theta[1]
      state <- pz_transition(state, alpha, t-1, c(theta[3:6]))
      states[,t+1] <- state
      log_obs[t] <- 0.2*randomness$y[t] + log(state[1,1])
    }
    return(log_obs)
  }
  ### additional functions to run PMMH
  rinit <- function(nparticles, theta, rand, ...){
    return(exp(matrix(log(2) + rand[1:(2*nparticles)], nrow = 2)))
  }
  rtransition <- function(xparticles, theta, time, rand, ...){
    nparticles <- ncol(xparticles)
    ra <- rand[((2+time-1)*nparticles + 1):((2+time)*nparticles)]
    alphas <- theta[1] + theta[2] * ra
    xparticles <- pz_transition(xparticles, alphas, time-1, c(theta[3:6]))
    return(xparticles)
  }
  dmeasurement <- function(xparticles, theta, observation, ...) {
    return(dnorm(x = observation, mean = log(xparticles[1,]), sd = 0.2, log = TRUE))
  }
  generate_randomness_pmmh <- function(nparticles, datalength){
    return(pz_generate_randomness_cpp(nparticles, datalength))
  }
  precompute <- function(...){
    return(list())
  }

  #
  model <- list(rprior = rprior,
                dprior = dprior,
                generate_randomness = generate_randomness,
                robservation = robservation,
                parameter_names = c("mu_alpha", "sigma_alpha", "c", "e", "ml", "mq"),
                thetadim = 6, ydim = 1,
                rinit = rinit, rtransition = rtransition, dmeasurement = dmeasurement, generate_randomness_pmmh = generate_randomness_pmmh, precompute = precompute)
  return(model)
}
