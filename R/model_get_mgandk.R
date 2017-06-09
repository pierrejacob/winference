#'@rdname get_mgandk
#'@title multivariate G and k model
#'@description This function returns
#'@return The list contains rprior, dprior (generate and evaluate the density of prior distribution),
#' generate_randomness (generate data-generating variables), robservation (create synthetic
#' data sets), parameter_names (useful for plotting), thetadim (dimension of parameter),
#' ydim (dimension of observations), parameters (list of hyperparameters,
#' to be passed to rprior,dprior,robservation)
#'@export
get_mgandk <- function(){
  rprior <- function(nparticles, parameters){
    thetas <- matrix(runif(nparticles*8, min = 0, max = 10), ncol = 8)
    thetas <- cbind(thetas, runif(nparticles, min = -1, max = 1))
    return(thetas)
  }
  # evaluate the log-density of the prior, for each particle
  dprior <- function(thetas, parameters){
    densities <- rep(0, nrow(thetas))
    for (i in 1:nrow(thetas)){
      if (any(thetas[i,1:8] > 10) || any(thetas[i,1:8] < 0)){
        densities[i] <- -Inf
      }
      densities[i] <- densities[i] + dunif(thetas[i,9], min = -1, max = 1, log = TRUE)
    }
    return(densities)
  }
  # function to compute a dataset for each theta value
  robservation <- function(nobservations, theta, ...){
    theta_1 <- theta[1:4]
    theta_2 <- theta[5:8]
    rho <- theta[9]
    covariance <- matrix(c(1, rho, rho, 1), ncol = 2, nrow = 2)
    normals <- fast_rmvnorm(nobservations, rep(0, 2), covariance)
    y_1 <- gandkinversecdf_givennormals(normals[,1], theta_1)
    y_2 <- gandkinversecdf_givennormals(normals[,2], theta_2)
    return(rbind(y_1, y_2))
  }
  gandk_loglikelihood <- function(thetas, ys, ...){
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
  #
  loglikelihood <- function(thetas, ys, ...){
    lls <- rep(0, nrow(thetas))
    for (itheta in 1:nrow(thetas)){
      Fy1 <- sapply(ys[1,], function(v) gandkcdf(v, thetas[itheta,1:4]))
      Fy2 <- sapply(ys[2,], function(v) gandkcdf(v, thetas[itheta,5:8]))
      x1 <- qnorm(Fy1)
      x2 <- qnorm(Fy2)
      covariance <- matrix(c(1, thetas[itheta,9], thetas[itheta,9], 1), ncol = 2, nrow = 2)
      lls[itheta] <- sum(fast_dmvnorm(cbind(x1, x2), rep(0, 2), covariance))
      lls[itheta] <- lls[itheta] - sum(dnorm(qnorm(Fy1), log = T)) - sum(dnorm(qnorm(Fy2), log = T))
    }
    lls <- lls + gandk_loglikelihood(thetas[,1:4,drop=F], ys[1,])
    lls <- lls + gandk_loglikelihood(thetas[,5:8,drop=F], ys[2,])
    return(lls)
  }
  parameters <- list()
  #
  model <- list(rprior = rprior,
                dprior = dprior,
                robservation = robservation,
                loglikelihood = loglikelihood,
                parameter_names = c("a1", "b1", "g1", "k1", "a2", "b2", "g2", "k2", "rho"),
                parameters = parameters,
                thetadim = 9, ydim = 2)
  return(model)
}
