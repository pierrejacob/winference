#'@export
levydriven_rtransition_rand <- function(nparticles, theta){
  levydriven_rtransition_rand_cpp(nparticles, theta)
}
