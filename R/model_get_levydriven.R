#'@rdname get_levydriven
#'@title Levy driven stochastic volatility model
#'@description This function returns a list representing a Levy driven stochastic volatility model.
#'@return The list contains rprior, dprior (generate and evaluate the density of prior distribution),
#' generate_randomness (generate data-generating variables), robservation (create synthetic
#' data sets), parameter_names (useful for plotting), thetadim (dimension of parameter),
#' ydim (dimension of observations), parameters (list of hyperparameters,
#' to be passed to rprior,dprior,robservation)
#'@export
get_levydriven <- function(){
  ## Stochastic volatility : one-factor model
  #
  # Y_t = mu + beta * v_t + v_t^(0.5) * epsilon_t
  # X_t = (v_t, z_t)
  # v_t+1 = lambda^(-1) ( z_t - z_{t+1} + sum_j=1^k e_j )
  # z_t+1 = e^(-lambda) * z_t + sum_j=1^k exp(-lambda(t + 1 - c_j)) e_j
  # k ~ Poisson(lambda * xi^2 / omega^2)
  # c_{1:k} ~ Uniform(t, t + 1)
  # e_{1:k} ~ Exp(xi / omega^2) (rate parameter)
  #
  # v_0 does not matter
  # z_0 ~ Gamma(xi^2 / omega^2, xi / omega^2) (rate parameter)
  #
  # theta <- c(0, 0, 0.5, 0.0625, 0.01)
  rprior <- function(nparticles, ...){
    theta1 <- rnorm(nparticles, 0, sd = sqrt(2))
    theta2 <- rnorm(nparticles, 0, sd = sqrt(2))
    theta3 <- rexp(nparticles, rate = 0.2)
    theta4 <- rexp(nparticles, rate = 0.2)
    theta5 <- rexp(nparticles, rate = 1)
    return(cbind(theta1, theta2, theta3, theta4, theta5))
  }
  dprior <- function(thetas, ...){
    if (is.null(dim(thetas))) thetas <- matrix(thetas, nrow = 1)
    density_evals <- dnorm(thetas[,1], mean = 0, sd = sqrt(2), log = TRUE)
    density_evals <- density_evals + dnorm(thetas[,2], mean = 0, sd = sqrt(2), log = TRUE)
    density_evals <- density_evals + dexp(thetas[,3], rate = 0.2, log = TRUE)
    density_evals <- density_evals + dexp(thetas[,4], rate = 0.2, log = TRUE)
    density_evals <- density_evals + dexp(thetas[,5], rate = 1, log = TRUE)
    return(density_evals)
  }

  generate_randomness <- function(nobservations){
    return(list())
  }
  robservation <- function(nobservations, theta, parameters, randomness){
    obs <- rep(0, nobservations)
    state <- rgamma(2, shape = theta[3] * theta[3] / theta[4], scale = theta[4]/theta[3])
    for (t in 1:nobservations){
      rtransition_r <- levydriven_rtransition_rand(1, theta)
      new_z <- exp(-theta[5]) * state[2] + rtransition_r$sum_weighted_e
      new_v <- (1/theta[5]) * (state[2] - new_z + rtransition_r$sum_e)
      state[1] <- new_v
      state[2] <- new_z
      obs[t] <- rnorm(1, mean = theta[1] + theta[2] * state[1], sd = sqrt(state[1]))
    }
    return(obs)
  }
  model <- list(rprior = rprior,
                dprior = dprior,
                generate_randomness = generate_randomness,
                robservation = robservation,
                parameter_names = c("mu", "beta", "xi", "omega2", "lambda"),
                thetadim = 5, ydim = 1)
  return(model)
}
