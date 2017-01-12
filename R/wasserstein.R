#'@rdname wasserstein
#'@title wasserstein
#'@description Compute regularized Wasserstein distance between two empirical distributions,
#' p and q, specified as vector of probabilities summing to one.
#' The third argument is the cost matrix, i.e. a matrix of pair-wise distances,
#' the fourth argument is the regularization parameter, e.g. 0.05*median(cost_matrix),
#' and the last argument is the number of Sinkhorn iterations to perform, e.g. 100.
#' Important references are
#'
#' - Cuturi, M. (2013). Sinkhorn distances: Lightspeed computation of optimal transport. In Advances in Neural Information Processing Systems (NIPS), pages 2292–2300.
#'
#' - Cuturi, M. and Doucet, A. (2014). Fast computation of Wasserstein barycenters. In Proceedings of the 31st International Conference on Machine Learning (ICML), pages 685–693.
#'
#'@return a list with "distances", "transportmatrix", "u" and "v"
#'@export
wasserstein <- function(p, q, cost_matrix, epsilon, niterations){
  return(wasserstein_(p, q, cost_matrix, epsilon, niterations))
}
