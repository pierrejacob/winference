#'@export
#'
wasserstein_auto <- function(p, q, cost_matrix, epsilon, desired_alpha){
  return(wasserstein_auto_(p, q, cost_matrix, epsilon, desired_alpha))
}
