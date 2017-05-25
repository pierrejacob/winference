# returns a function that computes MMD between y and z, for fixed y
#'@export
get_mmd_to_y <- function(y){
  nobs <- ncol(y)
  Cy1 <- cost_matrix_L1(y, y)
  Cy2 <- cost_matrix_L2(y, y)^2
  eps <- median(as.numeric(Cy1))
  k_y <- exp(-Cy2/(2*(eps^2)))
  first_term <- sum(k_y) / (nobs*nobs)
  f <- function(z){
    return(mmd_c(first_term, eps, z, y))
  }
  return(f)
}
# compute MMD between y and z
#'@export
mmd <- function(y, z){
  nobs <- ncol(y)
  Cy1 <- cost_matrix_L1(y, y)
  Cy2 <- cost_matrix_L2(y, y)^2
  eps <- median(as.numeric(Cy1))
  k_y <- exp(-Cy2/(2*(eps^2)))
  first_term <- sum(k_y) / (nobs*nobs)
  return(mmd_c(first_term, eps, z, y))
}
