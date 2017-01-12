#'@rdname cost_matrix_L1
#'@title cost_matrix_L1
#'@description Compute cost matrix L1 between two matrices of dimension d x N
#'@export
#'
cost_matrix_L1 <- function(x, y){
  return(cost_matrix_L1_(x,y))
}

#'@rdname cost_matrix_L2
#'@title cost_matrix_L2
#'@description Compute cost matrix L2 between two matrices of dimension d x N
#'@export
#'
cost_matrix_L2 <- function(x, y){
  return(cost_matrix_L2_(x,y))
}
