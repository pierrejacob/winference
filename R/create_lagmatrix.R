#'@rdname create_lagmatrix
#'@title create_lagmatrix
#'@description This function creates the delay reconstruction, i.e. a matrix
#' where the first row contains y_1, ..., y_T (the given univariate time series)
#' the second row contains NA, y_2, ..., y_T .... etc
#' the k-th row contains NA, NA, ..., y_{k+1}, ..., y_T
#'@return a matrix with k rows and n-k columns, where n is the length of the provided time series
#'@export
create_lagmatrix <- function(timeseries, k){
  if (k == 0){
    return(matrix(timeseries, nrow = 1))
  }
  res <- matrix(NA, nrow = k+1, ncol = ncol(timeseries))
  res[1,] <- timeseries
  for (lagvalue in 1:k){
    res[lagvalue+1,] <- lag(timeseries[1,], lagvalue)
  }
  return(res[,(k+1):ncol(res),drop=F])
}
