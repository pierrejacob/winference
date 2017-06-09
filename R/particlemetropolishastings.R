# Particle marginal Metropolis-Hastings
#'@export
pmmh <- function(observations, target, tuning_parameters, savefile = NULL){
  niterations <- tuning_parameters$niterations
  nchains <- tuning_parameters$nchains
  cov_proposal <- tuning_parameters$cov_proposal
  p <- ncol(tuning_parameters$init_chains)
  nparticles <- tuning_parameters$nparticles
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
  current_dprior <-  target$dprior(current_chains, target$parameters)
  current_lls <- rep(0, nchains)
  for (ichain in 1:nchains){
    if (!is.infinite(current_dprior[ichain])){
      current_lls[ichain] <- particle_filter(nparticles, target, current_chains[ichain,], observations)
    }
  }
  current_lls[is.na(current_lls)] <- -10^50 # so that
  #
  naccepts <- 0
  # run the chains
  for (iteration in 2:niterations){
    if (iteration %% 10 == 1){
      cat("iteration ", iteration, "/", niterations, "\n")
      cat("average acceptance:", naccepts / (iteration*nchains) * 100, "%\n")
    }
    if (iteration > 50 && tuning_parameters$adaptation > 0  && (iteration %% tuning_parameters$adaptation) == 0){
      # adapt the proposal covariance matrix based on the last < 50,000 samples of all chains
      mcmc_samples <- foreach(ichain = 1:nchains, .combine = rbind) %do% {
        matrix(chains[[ichain]][max(1, iteration - tuning_parameters$adaptation):(iteration-1),], ncol = p)
      }
      cov_proposal <- cov(mcmc_samples) / p
    }
    # proposals
    proposals <- current_chains + fast_rmvnorm(nchains, rep(0, p), cov_proposal)
    # proposals' prior density
    proposal_dprior <- target$dprior(proposals, target$parameters)
    proposal_lls <- rep(0, nchains)
    for (ichain in 1:nchains){
      if (!is.infinite(proposal_dprior[ichain])){
        proposal_lls[ichain] <- particle_filter(nparticles, target, proposals[ichain,], observations)
      }
    }
    proposal_dtarget <- proposal_lls + proposal_dprior
    proposal_dtarget[is.na(proposal_dtarget)] <- -Inf
    # log Metropolis Hastings ratio
    acceptance_ratios <- (proposal_dtarget - current_lls - current_dprior)
    # uniforms for the acceptance decisions
    uniforms <- runif(n = nchains)
    # acceptance decisions
    accepts <- rep(FALSE, nchains)
    for (ichain in 1:nchains){
      if (is.finite(proposal_dtarget[ichain])){
        if (log(uniforms[ichain]) < acceptance_ratios[ichain]){
          accepts[ichain] <- TRUE
        }
      }
    }
    naccepts <- naccepts + sum(accepts)
    # make the appropriate replacements
    current_chains[accepts,] <- proposals[accepts,]
    current_lls[accepts] <- proposal_lls[accepts]
    if (is.null(dim(current_chains))) current_chains <- matrix(current_chains, ncol = p)
    # current_dtarget[accepts] <- proposal_dtarget[accepts]
    # book keeping
    for (ichain in 1:nchains){
      chains[[ichain]][iteration,] <- current_chains[ichain,]
    }
    if (!is.null(savefile) && iteration %% 1000 == 1){
      pmmh_results <- list(chains = chains, naccepts = naccepts, cov_proposal = cov_proposal, iteration = iteration)
      save(pmmh_results, file = savefile)
    }
  }
  cat("average acceptance:", naccepts / (niterations*nchains) * 100, "%\n")
  return(list(chains = chains, naccepts = naccepts, cov_proposal = cov_proposal))
}
