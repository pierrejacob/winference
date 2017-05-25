#'@export
independent_proposal <- function(){
  proposal <- list(r = function(thetas, param_prop) fast_rmvnorm(nrow(thetas), param_prop$mean, param_prop$cov),
                   d = function(thetas, param_prop) fast_dmvnorm(thetas, param_prop$mean, param_prop$cov),
                   param_update = function(thetas, weights){
                     m <- wmean(thetas, weights)
                     c <- wcovariance(thetas, weights, m)
                     return(list(mean = m, cov = c))
                   })
  return(proposal)
}

#'@export
randomwalk_proposal <- function(){
  proposal <- list(r = function(thetas, param_prop) thetas + fast_rmvnorm(nrow(thetas), rep(0, ncol(thetas)), param_prop$cov),
                   d = function(thetas, param_prop) rep(0, nrow(thetas)),
                   param_update = function(thetas, weights){
                     m <- wmean(thetas, weights)
                     c <- wcovariance(thetas, weights, m)
                     return(list(mean = m, cov = c))
                   })
  return(proposal)
}

#'@export
mixture_mclust <- function(nclust = 5){
  library(mclust)
  proposal <- list(r = function(thetas, param_prop){
    fit <- param_prop$fit
    x <- sim(modelName = fit$modelName, parameters = fit$parameters, n = nrow(thetas))
    x <- x[,2:ncol(x)]
    return(matrix(x, ncol = ncol(thetas)))
  },
  d = function(thetas, param_prop){
    fit <- param_prop$fit
    z <- dens(modelName = fit$modelName, data = thetas, logarithm = TRUE, parameters = fit$parameters)
    return(z)
  },
  param_update = function(thetas, weights){
    ancestors <- systematic_resampling(weights)
    thetas_check <- matrix(thetas[ancestors,], ncol = ncol(thetas))
    fit <- Mclust(thetas_check, G = nclust)
    return(list(fit = fit))
  })
  return(proposal)
}

#'@export
mixture_rmixmod <- function(nclust = 5){
  library(Rmixmod)
  proposal <- list(r = function(thetas, param_prop){
    proportions <- param_prop$fit@bestResult@parameters@proportions
    means <- param_prop$fit@bestResult@parameters@mean
    variances <- param_prop$fit@bestResult@parameters@variance
    K <- nrow(means)
    n <- nrow(thetas)
    X <- matrix(0, nrow = n, ncol = ncol(means))
    # sample allocations
    allocations <- systematic_resampling_n(proportions, n)
    for (k in 1:K){
      which.k <- which(allocations == k)
      nk <- length(which.k)
      if (nk > 0){
        xk <- fast_rmvnorm(nk, means[k,], variances[[k]])
        X[which.k,] <- xk
      }
    }
    # random shuffling
    X <- X[sample(x = 1:n, size = n, replace = FALSE),,drop=FALSE]
    return(X)
  },
  d = function(thetas, param_prop){
    proportions <- param_prop$fit@bestResult@parameters@proportions
    means <- param_prop$fit@bestResult@parameters@mean
    variances <- param_prop$fit@bestResult@parameters@variance
    n <- nrow(thetas)
    d <- ncol(thetas)
    K <- nrow(means)
    evals <- matrix(0, nrow = n, ncol = K)
    for (k in 1:K){
      evals[,k] <- fast_dmvnorm(x = thetas, mean = means[k,], covariance = variances[[k]]) + log(proportions[k])
    }
    g <- function(row){
      m <- max(row)
      return(m + log(sum(exp(row - m))))
    }
    results <- apply(X = evals, MARGIN = 1, FUN = g)
    return(results)
  },
  param_update = function(thetas, weights){
    # print("fitting mixture...")
    ancestors <- systematic_resampling(weights)
    thetas_check <- thetas[ancestors,,drop=FALSE]
    fit <- mixmodCluster(data = data.frame(thetas_check), nbCluster = nclust, dataType = "quantitative")
    # print("mixture fitted")
    return(list(fit = fit))
  })
  return(proposal)
}
