#' @rdname get_toggleswitch
#' @title Toggle switch model
#' @description This function returns a list representing the toggle switch model
#' of Bonassi, F. V., West, M., et al. (2015).
#' Sequential Monte Carlo with adaptive weights for approximate Bayesian computation. Bayesian Analysis, 10(1):171–187.
#' @return The list contains rprior, dprior (generate and evaluate the density of prior distribution),
#' generate_randomness (generate data-generating variables), robservation (create synthetic
#' data sets), parameter_names (useful for plotting), thetadim (dimension of parameter),
#' ydim (dimension of observations), parameters (list of hyperparameters,
#' to be passed to rprior,dprior,robservation)
#' @export
get_toggleswitch <- function(){
  library(truncnorm)
  rprior <- function(nparticles, parameters){
    particles <- matrix(nrow = nparticles, ncol = 7)
    particles[,1] <- runif(nparticles, 0, 50)
    particles[,2] <- runif(nparticles, 0, 50)
    particles[,3] <- runif(nparticles, 0, 5)
    particles[,4] <- runif(nparticles, 0, 5)
    particles[,5] <- runif(nparticles, 250, 450)
    particles[,6] <- runif(nparticles, 0, 0.5)
    particles[,7] <- runif(nparticles, 0, 0.4)
    return(particles)
  }
  # evaluate the log-density of the prior, for each particle
  dprior <- function(thetaparticles, parameters){
    logdensities <- rep(0, nrow(thetaparticles))
    logdensities <- dunif(thetaparticles[,1], min = 0, max = 50, log = TRUE)
    logdensities <- logdensities + dunif(thetaparticles[,2], min = 0, max = 50, log = TRUE)
    logdensities <- logdensities + dunif(thetaparticles[,3], min = 0, max = 5, log = TRUE)
    logdensities <- logdensities + dunif(thetaparticles[,4], min = 0, max = 5, log = TRUE)
    logdensities <- logdensities + dunif(thetaparticles[,5], min = 250, max = 450, log = TRUE)
    logdensities <- logdensities + dunif(thetaparticles[,6], min = 0, max = 0.5, log = TRUE)
    logdensities <- logdensities + dunif(thetaparticles[,7], min = 0, max = 0.4, log = TRUE)
    return(logdensities)
  }
  # generate random variables used to compute a synthetic dataset
  generate_randomness <- function(nobservations){
    return(list())
  }
  # function to compute a dataset for each theta value
  robservation <- function(nobservations, theta, parameters, randomness){
    # constants used in the data generating process
    h <- 1
    tau <- 300
    u0 <- 10
    v0 <- 10
    #
    u <- matrix(0, nrow = nobservations, ncol = tau+1)
    v <- matrix(0, nrow = nobservations, ncol = tau+1)
    u[,1] <- u0
    v[,1] <- v0
    # noise <- array(randomness[(nobservations+1):length(randomness)], dim = c(nobservations, tau, 2))
    for (t in 1:tau){
      u[,t+1] <- u[,t] + h * theta[1] / (1 + v[,t]^theta[3]) - h * (1 + 0.03 * u[,t])
      u[,t+1] <- u[,t+1] + h * 0.5 * rtruncnorm(nobservations, a = - u[,t+1]/(h*0.5))
      v[,t+1] <- v[,t] + h * theta[2] / (1 + u[,t]^theta[4]) - h * (1 + 0.03 * v[,t])
      v[,t+1] <- v[,t+1] + h * 0.5 * rtruncnorm(nobservations, a = - v[,t+1]/(h*0.5))
    }
    lb = -(u[,tau+1] + theta[5]) / (theta[5] * theta[6]) * (u[,tau+1]^theta[7])
    y <- u[,tau+1] + theta[5] + theta[5] * theta[6] * rtruncnorm(nobservations, a = lb) / (u[,tau+1]^theta[7])
    return(y)
  }
  parameters <- list()
  #
  model <- list(rprior = rprior,
                dprior = dprior,
                generate_randomness = generate_randomness,
                robservation = robservation,
                parameter_names = c("alpha_1", "alpha_2", "beta_1", "beta_2", "mu", "sigma", "gamma"),
                parameters = parameters,
                thetadim = 7, ydim = 1)
  return(model)
}
