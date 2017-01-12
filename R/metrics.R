# compute distance between vectors of same size and dimension 1
#'@export
metricL1 <- function(xvec,yvec)  sum(abs(xvec - yvec)) / length(xvec)

#'@export
metricL2 <- function(xvec,yvec)  sqrt(sum((xvec - yvec)^2)) / length(xvec)
