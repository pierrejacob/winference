target <- list()

target$generate_randomness <- function(nobservations){
  return(rnorm(nobservations))
}

target$thetadim = 2

metricL1 = function(xvec,yvec)  mean(abs(xvec - yvec))
