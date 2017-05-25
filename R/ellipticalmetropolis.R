#'@export
ellipticalmetropolis <- function(niterations, nobservations, target, proposal, compute_d, epsilon, thetas, distances, us, maxangle = 6.283185, savefile = NULL){
  # store the chains in big matrix, with
  # parameters, then index of chain, then index of iteration
  nchains <- nrow(thetas)
  thetadim <- target$thetadim
  # store whole chains
  chains <- rep(list(matrix(nrow = niterations, ncol = thetadim)), nchains)
  # current states of the chains
  # current_chains <- thetas
  # initialization of the chains
  for (ichain in 1:nchains){
    chains[[ichain]][1,] <- thetas[ichain,]
  }
  # also store the associated iistances
  chain_distances <- matrix(nrow = niterations, ncol = nchains)
  chain_distances[1,] <- distances
  # and the acceptance rates
  acceptance_rates <- rep(0, niterations)
  nevals <- rep(0, nchains)
  widths <- matrix(0, nrow = niterations, ncol = nchains)
  for (iteration in 2:niterations){
    if (iteration %% 100 == 1){
      cat("iteration", iteration, "/", niterations, ", accceptance rate =", 100*mean(acceptance_rates[1:iteration]), "%\n")
    }
    ############################
    ############## update theta given u
    thetas_prop <- proposal$r(thetas, proposal$param_prop)
    prior_proposed <- target$dprior(thetas_prop, target$parameters)
    prop_proposed <- proposal$d(thetas_prop, proposal$param_prop)
    #
    prior_current <- target$dprior(thetas, target$parameters)
    prop_current <- proposal$d(thetas, proposal$param_prop)
    logratios <- (prior_proposed - prior_current) + (prop_current - prop_proposed)
    dist_proposed <- rep(Inf, nchains)
    potential_indices <- which(!is.infinite(logratios))
    if (length(potential_indices) > 0){
      # res_foreach <- foreach(ipot = potential_indices, .combine = c) %dopar% {
        # compute_d(thetas_prop[ipot,], us[ipot,])
      # }
      # dist_proposed[potential_indices] <- as.numeric(res_foreach[,1])
      for (ipot in potential_indices){
        dist_proposed[ipot] <- compute_d(thetas_prop[ipot,], us[ipot,])
      }
    }
    logratios <- logratios + log(dist_proposed < epsilon)
    accepted <- (log(runif(nchains)) < logratios)
    #
    thetas[accepted,] <- thetas_prop[accepted,,drop=F]
    distances[accepted] <- dist_proposed[accepted]
    ############################
    ############## update u given theta
    for (ichain in 1:nchains){
      theta <- thetas[ichain,]
      f <- us[ichain,]
      logLf <- 0 # otherwise wouldn't have been accepted before
      ellipt_result <- elliptical_update(f, logLf, nobservations, theta, compute_d, epsilon, maxangle)
      us[ichain,] <- ellipt_result$fprime
      nevals[ichain] <- nevals[ichain]  + ellipt_result$nevals
      widths[iteration,ichain] <- ellipt_result$width
    }
    ##############
    for (ichain in 1:nchains){
      chains[[ichain]][iteration,] <- thetas[ichain,]
    }
    chain_distances[iteration,] <- distances
    acceptance_rates[iteration] <- mean(accepted)
    if (!is.null(savefile) && (iteration %% 100 == 1)){
      save_chains <- chains
      for (ichain in 1:nchains){
        save_chains[[ichain]] <- save_chains[[ichain]][1:iteration,,drop=F]
      }
      mh_elliptical_res <- list(chains = save_chains, chain_distances = chain_distances[1:iteration,], nevals = nevals,
                                widths  = widths[1:iteration,], acceptance_rate = sum(acceptance_rates)/iteration)
      save(mh_elliptical_res, file = savefile)
    }
  }
  cat("average acceptance rate: ", 100*mean(acceptance_rates), "%\n")
  return(list(chains = chains, chain_distances = chain_distances, nevals = nevals,
              widths  = widths, acceptance_rate = mean(acceptance_rates)))
}

#'@export
elliptical_update <- function(f, logLf, nobservations, theta, compute_d, epsilon, maxangle = 6.283185){
  nu <- rnorm(nobservations)
  uniform <- runif(1, min = 0, max = 1)
  logy <- logLf + log(uniform)
  angle <- runif(1, min = 0, max = maxangle)
  angle_min <- angle - maxangle
  angle_max <- angle
  accept <- FALSE
  nevals <- 0
  fprime <- rep(0, nobservations)
  while (!accept && nevals < 1e3){
    fprime <- f * cos(angle) + nu * sin(angle)
    cd <- compute_d(theta, fprime)
    logLfprime <- log(cd < epsilon)
    # logLfprime <- log_target_given_theta(fprime, theta)
    nevals <- nevals + 1
    if (logLfprime > logy){
      accept <- TRUE
    } else {
      if (angle < 0){
        angle_min <- angle
      } else {
        angle_max <- angle
      }
      angle <- runif(1, min = angle_min, max = angle_max)
    }
  }
  if (nevals == 1e3){
    print("elliptical update did 1000 evaluations")
  }
  return(list(nevals = nevals, logLfprime = logLfprime, fprime = fprime, width = angle_max - angle_min))
}
