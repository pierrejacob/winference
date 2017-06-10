# Function to perform Metropolis-Hastings
#'@export
metropolishastings <- function(observations, target, tuning_parameters, savefile = NULL, verbose = FALSE){
  # Posterior density function (log)
  posterior <- function(thetas){
    logdens <- target$dprior(thetas, target$parameters)
    which.ok <- which(is.finite(logdens))
    if (length(which.ok) > 0){
      theta.ok <- thetas[which.ok,,drop=FALSE]
      logdens[which.ok] <- logdens[which.ok] + target$loglikelihood(theta.ok, observations, target$parameters)
    }
    return(logdens)
  }

  niterations <- tuning_parameters$niterations
  nchains <- tuning_parameters$nchains
  cov_proposal <- tuning_parameters$cov_proposal
  p <- ncol(tuning_parameters$init_chains)

  # store whole chains
  chains <- rep(list(matrix(nrow = niterations, ncol = p)), nchains)
  # current states of the chains
  current_chains <- matrix(nrow = nchains, ncol = p)
  # initialization of the chains
  current_chains <- matrix(tuning_parameters$init_chains, nrow = nchains, ncol = p)
  for (ichain in 1:nchains){
    chains[[ichain]][1,] <- current_chains[ichain,]
  }
  # log target density values associated with the current states of the chains
  current_dtarget <-  posterior(current_chains)
  #
  naccepts <- 0
  # run the chains
  for (iteration in 2:niterations){
    if ((iteration %% max(1, floor(niterations/100)) == 1) && (verbose)){
      cat("iteration ", iteration, "/", niterations, "\n")
      cat("average acceptance:", naccepts / (iteration*nchains) * 100, "%\n")
    }
    if (iteration > 250 && tuning_parameters$adaptation > 0  && (iteration %% tuning_parameters$adaptation) == 0){
      # adapt the proposal covariance matrix based on the last < 50,000 samples of all chains
      mcmc_samples <- foreach(ichain = 1:nchains, .combine = rbind) %do% {
        matrix(chains[[ichain]][max(1, iteration - tuning_parameters$adaptation):(iteration-1),], ncol = p)
      }
      cov_proposal <- cov(mcmc_samples) / p
    }
    # proposals
    proposals <- current_chains + fast_rmvnorm(nchains, rep(0, p), cov_proposal)
    # proposals' target density
    proposal_dtarget <- posterior(proposals)
    # log Metropolis Hastings ratio
    acceptance_ratios <- (proposal_dtarget - current_dtarget)
    # uniforms for the acceptance decisions
    uniforms <- runif(n = nchains)
    # acceptance decisions
    accepts <- (log(uniforms) < acceptance_ratios)
    naccepts <- naccepts + sum(accepts)
    # make the appropriate replacements
    current_chains[accepts,] <- proposals[accepts,]
    if (is.null(dim(current_chains))) current_chains <- matrix(current_chains, ncol = p)
    current_dtarget[accepts] <- proposal_dtarget[accepts]
    # book keeping
    for (ichain in 1:nchains){
      chains[[ichain]][iteration,] <- current_chains[ichain,]
    }
    if (!is.null(savefile) && iteration %% 1000 == 1){
      mh_results <- list(chains = chains, naccepts = naccepts, cov_proposal = cov_proposal, iteration = iteration)
      save(mh_results, file = savefile)
    }
  }
  cat("average acceptance:", naccepts / (niterations*nchains) * 100, "%\n")
  return(list(chains = chains, naccepts = naccepts, cov_proposal = cov_proposal))
}

#'@export
mhchainlist_to_dataframe <- function(chains_list){
  nchains <- length(chains_list)
  niterations <- nrow(chains_list[[1]])
  chaindf <- foreach (i = 1:nchains, .combine = rbind) %do% {
    data.frame(ichain = rep(i, niterations), iteration = 1:niterations, X = chains_list[[i]])
  }
  return(chaindf)
}
