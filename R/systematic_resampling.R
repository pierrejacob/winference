#'@rdname systematic_resampling
#'@title systematic_resampling
#'@description systematic_resampling
#'@export
systematic_resampling <- function(normalized_weights){
  return(systematic_resampling_(normalized_weights))
}

#'@export
systematic_resampling_given_u <- function(normalized_weights, u){
  return(systematic_resampling_n_(normalized_weights, length(normalized_weights), u))
}

#'@export
systematic_resampling_n <- function(normalized_weights, ndraws){
  return(systematic_resampling_n_(normalized_weights, ndraws, runif(1)))
}

# systematic_resampling <- function(normalized_weights){
#   N <- length(normalized_weights)
#   indices <- rep(0, N)
#   normalized_weights <- N * normalized_weights
#   j <- 1
#   csw <- normalized_weights[1]
#   u <- runif(1, min = 0, max = 1)
#   for (k in 1:N){
#     while (csw < u){
#       j <- j + 1
#       csw <- csw + normalized_weights[j]
#     }
#     indices[k] <- j
#     u <- u + 1
#   }
#   return(indices)
# }
