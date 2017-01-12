#'@rdname wcovariance
#'@title wcovariance
#'@description wcovariance
#'@export
wcovariance <- function(xparticles, normweights, mean){
  return(wcovariance_(xparticles, normweights, mean))
}
