#'@rdname hilbert_order
#'@title hilbert_order
#'@description This function returns the "Hilbert order" of a sample of n
#'points of dimension d, stored in a matrix with d rows and n columns
#' where d is the dimension of each sample and n the number of samples.
#' The function essentially calls Hilbert_Sort_CGAL, from the CGAL library.
#'@return a vector of index corresponding to the ordered samples
#'@export
hilbert_order <- function(x){
  return(hilbert_order_(x) + 1)
}
