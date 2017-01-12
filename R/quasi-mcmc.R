# One step of pseudo-marginal MH using the quasilikelihood provided, for each given theta
# uses a Gaussian random walk proposal
#'@export
quasimh_step <- function(thetas, priorvalues, qvalues, distances, threshold, qloglikelihood, target, param_algo){
  proposal <- param_algo$proposal
  nchains <- nrow(thetas)
  if (is.null(param_algo$parallel)) param_algo$parallel <- FALSE
  # if (ncol(qvalues) != param_algo$M+1){
  #   qq <- matrix(0, nrow = nrow(thetas), ncol = param_algo$M+1)
  #   qq[,1:ncol(qvalues)] <- qvalues
  #   qvalues <- qq
  # }
  # proposal
  thetas_prop <- proposal$r(thetas, proposal$param_prop)
  prior_proposed <- target$dprior(thetas_prop, target$parameters)
  prop_proposed <- proposal$d(thetas_prop, proposal$param_prop)
  # proposal density evaluated at current values
  prop_current <- proposal$d(thetas, proposal$param_prop)
  #
  logratios <- (prior_proposed - priorvalues) + (prop_current - prop_proposed)
  # find proposal that are admissible with respect to prior
  potential_indices <- which(!is.infinite(logratios))
  naccepts <- 0
  if (length(potential_indices) > 0){
    for (ipot in potential_indices){
      qresult <- qloglikelihood(thetas_prop[ipot,], threshold = threshold, M = param_algo$M)
      logratio <- logratios[ipot] + qresult$qvalue - qvalues[ipot]
      accepted <- (log(runif(1)) < logratio)
      if (accepted){
        naccepts <- naccepts + 1
        thetas[ipot,] <- thetas_prop[ipot,]
        priorvalues[ipot] <- prior_proposed[ipot]
        qvalues[ipot] <- qresult$qvalue
        distances[ipot,] <- qresult$distances
      }
    }
  }
  return(list(acceptrate = naccepts/nchains, thetas = thetas, qvalues = qvalues, distances = distances, priorvalues = priorvalues))
}

# A number nmoves of steps of pseudo-marginal MH using the quasilikelihood provided, for each given theta
# potentially adapts the covariance matrix of the Gaussian random walk proposal
#'@export
quasimh_steps <- function(thetas, priorvalues, qvalues, distances, threshold, qloglikelihood, target, param_algo){
  nmoves <- param_algo$nmoves
  nchains <- nrow(thetas)
  if (is.null(param_algo$adaptation)) param_algo$adaptation <- 0
  # if (is.null(param_algo$parallel)) param_algo$parallel <- FALSE
  # store histories
  chains_ <- rep(list(matrix(nrow = nmoves, ncol = target$thetadim)), nchains)
  qchains_ <- matrix(0, nrow = nmoves, ncol = nchains)
  distances_ <- rep(list(matrix(nrow = nmoves, ncol = param_algo$M)), nchains)
  priorchains_ <- matrix(nrow = nmoves, ncol = nchains)
  #
  accepts <- 0
  for (imove in 1:nmoves){
    if (imove %% 100 == 1){
      cat("iteration ", imove, "/", nmoves, "\n")
      cat("average acceptance:", accepts / imove * 100, "%\n")
    }
    if (imove > 100 && param_algo$adaptation > 0  && (imove %% param_algo$adaptation) == 0){
      # adapt the proposal covariance matrix based on the last < 50,000 samples of all chains
      mcmc_samples <- foreach(ichain = 1:nchains, .combine = rbind) %do% {
        matrix(chains_[[ichain]][max(1, imove - 50000):(imove-1),], ncol = target$thetadim)
      }
      new_cov <- cov(mcmc_samples) / target$thetadim
      if (sum(new_cov) > 1e-5){
        param_algo$proposal$param_prop$cov <- new_cov
      }
    }
    res <- quasimh_step(thetas, priorvalues, qvalues, distances, threshold, qloglikelihood, target, param_algo)
    thetas <- res$thetas
    qvalues <- res$qvalues
    distances <- res$distances
    priorvalues <- res$priorvalues
    accepts <- accepts + res$acceptrate
    # storing things
    priorchains_[imove,] <- priorvalues
    qchains_[imove,] <- qvalues
    for (ichain in 1:nchains){
      chains_[[ichain]][imove,] <- thetas[ichain,]
      distances_[[ichain]][imove,] <- distances[ichain,]
    }
  }
  chains_df <- foreach (i = 1:nchains, .combine = rbind) %do% {
    data.frame(ichain = rep(i, nmoves), iteration = 1:nmoves, X = chains_[[i]])
  }
  return(list(chains_list = chains_, chains_df = chains_df, qchains = qchains_, priorchains = priorchains_,
              distancechains = distances_, param_algo = param_algo))
}

# draw thetas from prior until quasi likelihood is non-zero for the given threshold
# (only makes "maxattempts" attempts per theta)
#'@export
initialization <- function(threshold, qloglikelihood, target, param_algo, maxattempts = 1000){
  nchains <- param_algo$nchains
  thetas <- target$rprior(nchains, parameters = target$parameters)
  priorvalues <- target$dprior(thetas, target$parameters)
  qvalues <- rep(Inf, nchains)
  distances <- matrix(0, nrow = nchains, ncol = param_algo$M)
  for (i in 1:nchains){
    attempts <- 0
    while(is.infinite(qvalues[i]) && attempts < maxattempts){
      attempts <- attempts + 1
      thetas[i,] <- target$rprior(1, parameters = target$parameters)
      qresult <- qloglikelihood(thetas[i,], threshold = threshold, M = param_algo$M)
      qvalues[i] <- qresult$qvalue
      distances[i,] <- qresult$distances
    }
    if (attempts == 1000){
      print("too hard to initialize")
    }
  }
  return(list(thetas = thetas, priorvalues = priorvalues, qvalues = qvalues, distances = distances))
}


# Performs pseudo-marginal MH using the quasilikelihood provided
# potentially adapts the covariance matrix of the Gaussian random walk proposal
# initialize using the initialization function
#'@export
quasimh <- function(threshold, qloglikelihood, target, param_algo, maxinitattempts = 1000){
  # first initialize from the prior distribution
  init <- initialization(threshold, qloglikelihood, target, param_algo, maxattempts = maxinitattempts)
  if (any(is.infinite(init$qvalues))){
    print("initialization failed, threshold might be too small")
    return(NULL)
  }
  # then run mh_steps
  res <- quasimh_steps(init$thetas, init$priorvalues, init$qvalues, init$distances, threshold, qloglikelihood, target, param_algo)
  return(res)
}

# take results of quasimh or quasimh_steps, and find theta values in there such that their
# quasi likelihood given a new threshold would still be non-zero
#'@export
find_init_values <- function(mh_results, threshold){
  nmoves <- nrow(mh_results$priorchains)
  nchains <- ncol(mh_results$priorchains)
  # load distances
  distances <- mh_results$distancechains
  M <- ncol(distances[[1]])
  ds <- foreach (i = 1:nchains, .combine = rbind) %do% {
    data.frame(ichain = rep(i, nmoves), iteration = 1:nmoves, q = distances[[i]])
  }
  # find index of points corresponding to nonzero quasi likelihood under new threshold
  potential_indices <- c()
  qliks <- apply(ds, 1, function(v) mean(v[3:(2+M)] < threshold))
  potential_indices <- which(qliks > 0)
  # hopefully there are some points
  if (length(potential_indices) == 0){
    print("no matching index for this threshold")
    return(NULL)
  }
  if (length(potential_indices) > nchains){
    selected_indices <- as.numeric(sample(x = potential_indices, size = nchains, replace = FALSE))
  } else {
    selected_indices <- as.numeric(sample(x = potential_indices, size = nchains, replace = TRUE))
  }
  thetas <- matrix(0, nrow = nchains, ncol = target$thetadim)
  priorvalues <- rep(0, nchains)
  qvalues <- rep(0, nchains)
  distances <- matrix(0, nrow = nchains, ncol = M)
  for (i in 1:nchains){
    ichain <- ds$ichain[selected_indices[i]]
    iteration <- ds$iteration[selected_indices[i]]
    qvalues[i] <- log(qliks[selected_indices[i]])
    priorvalues[i] <- mh_results$priorchains[iteration, ichain]
    thetas[i,] <- mh_results$chains_list[[ichain]][iteration,]
    distances[i,] <- mh_results$distancechains[[ichain]][iteration,]
  }
  return(list(thetas = thetas, priorvalues = priorvalues, qvalues = qvalues, distances = distances))
}

# function to continue the MCMC sampling from a previous result (obtained e.g. with quasimh)
# but with a new threshold, and potentially a new value for M, for the number of chains, etc
#'@export
quasimh_continue <- function(mh_results, threshold, qloglikelihood, target, param_algo){
  # find nchains with appropriate init value given the new threshold
  init <- find_init_values(mh_results, threshold)
  thetas <- init$thetas
  priorvalues <- init$priorvalues
  qvalues <- init$qvalues
  distances <- init$distances

  if (is.null(init)){
    print("cannot initialize properly the chains for that threshold; using latest states instead")
    nmoves <- mh_results$param_algo$nmoves
    nchains <- mh_results$param_algo$nchains
    thetas <- matrix(0, nrow = nchains, ncol = target$thetadim)
    distances <- matrix(0, nrow = nchains, ncol = mh_results$param_algo$M)
    for (ichain in 1:nchains){
      thetas[ichain,] <- mh_results$chains_list[[ichain]][nmoves,]
      distances[ichain,] <- mh_results$distancechains[[ichain]][nmoves,]
    }
    qvalues <- mh_results$qchains[nmoves,]
    priorvalues <- mh_results$priorchains[nmoves,]
  }
  # if in fact there's a new value of nchains...
  if (param_algo$nchains != mh_results$param_algo$nchains){
    ind <- sample(x = 1:mh_results$param_algo$nchains, size = param_algo$nchains, replace = TRUE)
    thetas <- thetas[ind,]
    priorvalues <- priorvalues[ind]
    qvalues <- qvalues[ind]
    distances <- distances[ind,]
  }
  if (param_algo$M != mh_results$param_algo$M){
    ind <- sample(x = 1:mh_results$param_algo$M, size = param_algo$M, replace = TRUE)
    distances <- distances[,ind]
  }
  # then run mh_steps
  res <- quasimh_steps(thetas, priorvalues, qvalues, distances, threshold, qloglikelihood, target, param_algo)
  return(res)
}
