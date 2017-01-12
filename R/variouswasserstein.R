#' #'@export
#' wasserstein <- function(p, qs, cost_matrix, epsilon, niterations){
#'   return(wasserstein_(p, qs, cost_matrix, epsilon, niterations))
#' }

#'@export
wasserstein_SAG <- function(p, q, cost_matrix, epsilon, niterations, stepsize){
  return(wasserstein_SAG_(p, q, cost_matrix, epsilon, niterations, stepsize))
}

#'@export
wasserstein_SAG_auto <- function(p, q, cost_matrix, epsilon, tolerance, stepsize, maxiterations){
  return(wasserstein_SAG_auto_(p, q, cost_matrix, epsilon, tolerance, stepsize, maxiterations))
}
