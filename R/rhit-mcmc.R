#'@export
mhrhit_step_onetheta <- function(theta, distance, compute_d, target, param_algo){
  threshold <- param_algo$threshold
  proposal <- param_algo$proposal
  R <- param_algo$R
  maxtrials <- param_algo$maxtrials
  if (is.null(maxtrials)){
    maxtrials <- Inf
  }
  #
  theta <- matrix(theta, nrow = 1)
  #
  prior_current <- target$dprior(theta, target$parameters)
  prop_current <- proposal$d(theta, proposal$param_prop)
  #
  nhits <- 0
  nproposals <- 0
  hit_indices <- c()
  # sample thetas and associated distances until R hits
  thetas_prop <- c()
  associated_dists <- c()
  # also compute prior densities as it will be required in case of
  # hit, and also allows to not sample from parameters outside the admissible range
  associated_priors <- c()
  #
  ncomputed <- 0
  while (nhits < R && nproposals < maxtrials){
    theta_prop <- proposal$r(theta, proposal$param_prop)
    prior_prop <- target$dprior(theta_prop, target$parameters)
    dproposed <- Inf
    if (!is.infinite(prior_prop)){
      dproposed <- compute_d(theta_prop)
      ncomputed <- ncomputed + 1
    }
    if (is.na(dproposed)){
      dproposed <- Inf
    }
    # if there's a hit
    if (dproposed < threshold){
      nhits <- nhits + 1
      # remember which parameter made the hit
      hit_indices <- c(hit_indices, nproposals+1)
    }
    nproposals <- nproposals + 1
    thetas_prop <- rbind(thetas_prop, theta_prop)
    associated_dists <- c(associated_dists, dproposed)
    associated_priors <- c(associated_priors, prior_prop)
  }
  # nproposals is N^prime in the notation of the paper
  # choose one of the succesful proposed parameter, except not the last one
  if (nproposals == maxtrials){
    cat("first while loop takes too long, let's reject\n")
    return(list(accepted = FALSE, theta = theta, distance = distance, nproposals = nproposals, ncurrent = NA, ncomputed = ncomputed))
  } else {
    L <- sample(x = hit_indices[1:(length(hit_indices)-1)], size = 1)
    if (length(hit_indices) == 2){
      # it took me ages to find this bug: if L is of size 1, then sample(x = L) will actually
      # sample from 1 to L[1] instead of returning L[1] with probability 1. Damn!
      L <- hit_indices[1]
    }
    theta_L <- matrix(thetas_prop[L,], nrow = 1)
    associated_prior_L <- associated_priors[L]
    associated_distance_L <- associated_dists[L]
    # now do it again from the selected parameter, until R - 1 hits
    nhits <- 0
    # ncurrent will be N in the notation of the paper
    ncurrent <- 0
    #
    while (nhits < (R-1) && ncurrent < maxtrials){
      theta_prop <- proposal$r(theta_L, proposal$param_prop)
      prior_prop <- target$dprior(theta_prop, target$parameters)
      if (!is.infinite(prior_prop[1])){
        dproposed <- compute_d(theta_prop)
        ncomputed <- ncomputed + 1
      } else {
        dproposed <- Inf
      }
      if (is.na(dproposed)) dproposed <- Inf
      if (dproposed < threshold){
        nhits <- nhits + 1
        hit_indices <- c(hit_indices, ncurrent+1)
      }
      ncurrent <- ncurrent + 1
    }
    if (ncurrent == maxtrials){
      cat("second while loop takes too long, let's reject\n")
      return(list(accepted = FALSE, theta = theta, distance = distance, nproposals = nproposals, ncurrent = ncurrent, ncomputed = ncomputed))
    } else {
      logratio <- associated_prior_L - (prior_current) +
        (prop_current - proposal$d(theta_L, proposal$param_prop)) + log(ncurrent/(nproposals-1))
      accepted <- (log(runif(1)) < logratio)
      if (accepted){
        theta <- theta_L
        distance <- associated_distance_L
        return(list(accepted = accepted, theta = theta, distance = distance, nproposals = nproposals, ncurrent = ncurrent, ncomputed = ncomputed))
      } else {
        return(list(accepted = FALSE, theta = theta, distance = distance, nproposals = nproposals, ncurrent = ncurrent, ncomputed = ncomputed))
      }
    }
  }
}

#'@export
mh_rhit_step <- function(thetas, distances, compute_d, target, param_algo){
  res_foreach <- foreach(i = 1:nrow(thetas), .combine = rbind) %dorng% {
    theta <- thetas[i,]
    distance <- distances[i]
    res <- mhrhit_step_onetheta(theta, distance, compute_d, target, param_algo)
    matrix(c(res$theta, res$distance, res$accepted, res$ncomputed), nrow = 1)
  }
  # print("finished")
  ncomputed <- sum(res_foreach[,(target$thetadim+3)])
  accepts <- res_foreach[,(target$thetadim+2)]
  thetas <- res_foreach[,(1:target$thetadim)]
  distances <- res_foreach[,(target$thetadim+1)]
  # accepts <- res_foreach[,(target$thetadim+2)]
  return(list(acceptrate = mean(accepts), thetas = thetas, distances = distances, ncomputed = ncomputed))
}


# multiple steps of MH r-hit kernel, for multiple chains (where each chain is parallelized)
#'@export
mhrhit_parallelchains <- function(nmoves, nchains, target, thetas_, distances_, param_algo){
  results <- foreach(ichain = 1:nchains, .combine = rbind) %dorng% {
    MCMC_results <- matrix(nrow = nmoves, ncol = target$thetadim+6)
    #
    theta_ <- thetas_[ichain,]
    distance_ <- distances_[ichain]
    #
    for (imove in 1:nmoves){
      mh_step_results <- mhrhit_step_onetheta(theta_, distance_, target, param_algo)
      theta_ <- mh_step_results$theta
      distance_ <- mh_step_results$distance
      MCMC_results[imove,1:target$thetadim] <- theta_
      MCMC_results[imove,(target$thetadim+1):(target$thetadim+6)] <- c(distance_, mh_step_results$accepted, ichain, imove,
                                                                       mh_step_results$nproposals, mh_step_results$ncurrent)
    }
    MCMC_results <- data.frame(MCMC_results)
    names(MCMC_results) <- c(paste0("X", 1:target$thetadim), "distance", "accept", "chain", "iteration", "nproposals", "ncurrent")
    MCMC_results
  }
  cat("average acceptance rate: ", 100*mean(results$accept), "%\n")
  return(results)
}
