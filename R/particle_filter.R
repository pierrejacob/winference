# requires model$rinit, model$rtransition, and model$dobs
# rinit specifies x_1 and the first observation is y_1
#'@export
particle_filter <- function(nparticles, model, theta, observations, storex = FALSE){
  datalength <- ncol(observations)
  # initialization
  xparticles <- model$rinit(nparticles, theta)
  logw <- model$dobs(observations, 1, xparticles, theta)
  if (all(is.infinite(logw))){
    return(NA)
  }
  maxlw <- max(logw)
  w <- exp(logw - maxlw)
  # update log likelihood estimate
  ll <- maxlw + log(mean(w))
  normweights <- w / sum(w)
  #
  if (storex){
    xhistory <- rep(list(matrix(ncol = nparticles, nrow = nrow(xparticles))), datalength)
    whistory <- rep(list(rep(0, nparticles)), datalength + 1)
    ahistory <- rep(list(rep(0, nparticles)), datalength)
    xhistory[[1]] <- xparticles
    whistory[[1]] <- normweights
  }
  # step t > 1
  for (time in 2:datalength){
    ancestors <- systematic_resampling_given_u(normweights, runif(1))
    xparticles <- xparticles[ancestors,,drop=F]
    xparticles <- model$rtransition(xparticles, time, theta)
    logw <- model$dobs(observations, time, xparticles, theta)
    if (all(is.infinite(logw))){
      return(NA)
    }
    maxlw <- max(logw)
    w <- exp(logw - maxlw)
    # update log likelihood estimate
    ll <- ll + maxlw + log(mean(w))
    normweights <- w / sum(w)

    #
    if (storex){
      ahistory[[time]] <- ancestors
      xhistory[[time+1]] <- xparticles
      whistory[[time+1]] <- normweights
    }
  }
  if (storex){
    return(list(ll = ll, xhistory = xhistory, whistory = whistory, ahistory = ahistory))
  } else {
    return(ll)
  }
}
