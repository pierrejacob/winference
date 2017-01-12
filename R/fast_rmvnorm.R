#'@rdname fast_rmvnorm
#'@title fast_rmvnorm
#'@description generate multivariate Normal draws
#'@export

fast_rmvnorm <- function(nparticles, mean, covariance){
  return(rmvnorm(nparticles, mean, covariance))
}
