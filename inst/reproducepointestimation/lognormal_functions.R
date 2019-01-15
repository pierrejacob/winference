
L <- 10

target <- list()

target$generate_randomness <- function(nobservations){
  return(rnorm(nobservations*L))
}

target$robservation <- function(theta, randomness){
  normals_ <- theta[1] + theta[2] * randomness
  lognormals_ <- exp(normals_)
  return(rowSums(matrix(lognormals_, ncol = L)))
}

metricL1 = function(xvec, yvec)  mean(abs(xvec - yvec))

true_theta <- c(0,1)

thetadim = 2 #Inference on all parameters

target$thetadim = 2
