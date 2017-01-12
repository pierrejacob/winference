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
mixture_proposal <- function(nclust = 5){
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
