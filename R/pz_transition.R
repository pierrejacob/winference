#'@rdname pz_transition
#'@title pz_transition
#'@description Solve PZ ODE for each particle, given each alpha, from time to time + 1,
#' and given the parameters (c, e, ml, mq).
#'@export
#'
pz_transition <- function(xparticles, alphas, time, parameters){
  return(one_step_pz_vector(xparticles, alphas, time, parameters))
}
