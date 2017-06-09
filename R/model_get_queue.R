#'@export
get_queue <- function(){
  rprior <- function(ntheta, parameters){
    theta1 <- runif(n = ntheta, min = 0, max = 10)
    theta2minus1 <- runif(n = ntheta, min = 0, max = 10)
    theta3 <- runif(n = ntheta, min = 0, max = 1/3)
    return(cbind(theta1, theta2minus1, theta3))
  }

  dprior <- function(thetas, parameters){
    evals <- dunif(thetas[,1], min = 0, max = 10, log = TRUE)
    evals <- evals + dunif(thetas[,2], min = 0, max = 10, log = TRUE)
    evals <- evals + dunif(thetas[,3], min = 0, max = 1/3, log = TRUE)
    return(evals)
  }

  robservation <- function(nobservations, theta, parameters, ...){
    theta1 <- theta[1]
    theta2 <- theta[2] + theta[1]
    theta3 <- theta[3]
    u <- runif(nobservations, theta1, theta2)
    v <- rep(0, nobservations)
    y <- rep(0, nobservations)
    x <- rep(0, nobservations)
    v[1] <- rexp(n = 1, rate = theta3)
    y[1] <- u[1] + v[1]
    x[1] <- y[1]
    for (t in 2:nobservations){
      v[t] <- v[t-1] + rexp(n = 1, rate = theta3)
      y[t] <- u[t] + max(0, v[t] - x[t-1])
      x[t] <- x[t-1] + y[t]
    }
    return(y)
  }
  model <- list(rprior = rprior,
                dprior = dprior,
                robservation = robservation,
                parameter_names = c("theta1", "theta2minus1", "theta3"),
                parameters = list(),
                thetadim = 3, ydim = 1)
  return(model)
}
