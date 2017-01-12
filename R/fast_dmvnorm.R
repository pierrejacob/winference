#'@rdname fast_dmvnorm
#'@title fast_dmvnorm
#'@description evaluate multivariate Normal log-densities
#'@export

fast_dmvnorm <- function(x, mean, covariance){
  return(dmvnorm(x, mean, covariance))
}
