#'@rdname wsmc
#'@title Adaptive SMC sampler with r-hit kernels
#'@description Runs an adaptive SMC sampler with r-hit kernels
#'
#' References for the r-hit kernels:
#' \itemize{
#' \item Lee, A. (2012). On the choice of MCMC kernels for approximate Bayesian computation with SMC samplers.
#'  In Proceedings of the 2012 Winter Simulation Conference, pages 304–315.
#' \item Lee, A. and Łatuszyński, K. (2014). Variance bounding and geometric ergodicity of Markov chain Monte Carlo kernels for approximate Bayesian computation. Biometrika, 101(3):655–671.
#'}
#'@return a list containing
#' \itemize{
#' \item thetas_history (all theta particles at all steps),
#'  \item distances_history, threshold_history, param_algo
#'  }
#'@export
wsmc <- function(compute_d, target, param_algo, savefile = NULL, parallel = TRUE, debug = FALSE, maxstep = Inf, maxtime = Inf, maxsimulation = Inf){
  if (is.infinite(maxstep) && is.infinite(maxtime) && is.infinite(maxsimulation)){
    print("error: you have to provide a finite value for 'maxstep' OR 'maxtime' OR 'maxsimulation'")
    return(NULL)
  }
  init_time <- proc.time()[3]
  param_algo$original_proposal <- param_algo$proposal
  if (is.null(param_algo$minimum_diversity)){
    param_algo$minimum_diversity <- 0.5
  }
  # sample from prior
  thetas <- target$rprior(param_algo$nthetas, target$parameters)
  # compute distances and associated datasets
  dy <- list()
  if (debug){
    print("compute first distances")
  }
  if (parallel){
    dy <- foreach(itheta = 1:param_algo$nthetas) %dorng% {
      y_fake <- target$simulate(thetas[itheta,])
      distance <- compute_d(y_fake)
      list(d = distance, y = y_fake)
    }
  } else {
    for(itheta in 1:param_algo$nthetas) {
      y_fake <- target$simulate(thetas[itheta,])
      distance <- compute_d(y_fake)
      dy[[itheta]] <- list(d = distance, y = y_fake)
    }
  }
  distances <- sapply(X = dy, FUN = function(x) x$d)
  latest_y <- lapply(X = dy, FUN = function(x) x$y)
  if (debug){
    print("first distances computed")
    print(summary(distances))
  }
  # distances[1:10,1:3]
  if (target$ydim == 1){
    nobs <- length(latest_y[[1]])
  } else {
    nobs <- ncol(latest_y[[1]])
  }
  distances[is.na(distances)] <- Inf
  # count of the number of distance calculations
  ncomputed <- c(param_algo$nthetas)
  normcst <- c()
  #
  param_algo$threshold <- as.numeric(quantile(distances, probs = param_algo$minimum_diversity))
  threshold_history <- param_algo$threshold
  if (debug){
    cat("first threshold: ", param_algo$threshold, "\n")
  }
  #
  thetas_history <- list()
  thetas_history[[1]] <- thetas
  distances_history <- list()
  distances_history[[1]] <- distances
  #
  step_time <- proc.time()[3]
  total_time <- step_time - init_time
  compute_times <- c(total_time)
  istep <- 1
  cat("step 1 completed, in", total_time, "seconds,", ncomputed, "distances calculated\n")
  while (total_time < maxtime && sum(ncomputed) < maxsimulation && istep < maxstep){
    istep <- istep + 1
    # then continue
    cat("step ", istep, "... running ", sep = "")
    if (is.finite(maxstep))
      cat("for", maxstep, "steps\n")
    if (is.finite(maxtime))
      cat("until time >=", maxtime, "seconds\n")
    if (is.finite(maxsimulation)){
      cat("until # distances calculated >=", maxsimulation, "\n")
    }
    one_step_results <- wsmc_one_step(thetas, distances, latest_y, compute_d, target, param_algo, debug = debug, parallel = parallel)
    #
    thetas <- one_step_results$thetas
    distances <- one_step_results$distances
    latest_y <- one_step_results$latest_y
    param_algo <- one_step_results$param_algo
    ncomputed <- c(ncomputed, one_step_results$ncomputed_thisstep)
    normcst <- c(normcst, one_step_results$normcst)
    cat("  total # distances calculated: ", sum(ncomputed), " (for this step: ", one_step_results$ncomputed_thisstep, ")\n", sep = "")

    # store
    thetas_history[[istep]] <- thetas
    distances_history[[istep]] <- distances
    threshold_history <- c(threshold_history, param_algo$threshold)
    #
    new_step_time <- proc.time()[3]
    total_time <- total_time + new_step_time - step_time
    cat("  total time spent:", total_time, "seconds (for this step:", new_step_time - step_time, "seconds)\n")
    step_time <- new_step_time
    compute_times <- c(compute_times, total_time)
    if (!is.null(savefile)){
      results <- list(thetas_history = thetas_history, distances_history = distances_history,
                      threshold_history = threshold_history,
                      latest_y = latest_y,
                      ncomputed = ncomputed, param_algo = param_algo, compute_d = compute_d,
                      target = target, normcst = normcst, compute_times = compute_times)
      save(results, file = savefile)
    }
    if (total_time > maxtime){
      break
    }
  }
  return(list(thetas_history = thetas_history, distances_history = distances_history,
              threshold_history = threshold_history, param_algo = param_algo,
              latest_y = latest_y, ncomputed = ncomputed, compute_d = compute_d, target = target, normcst = normcst,
              compute_times = compute_times))
}

#'@export
wsmc_continue <- function(results, savefile = NULL, maxstep = Inf, maxtime = Inf, maxsimulation = Inf){
  if (is.infinite(maxstep) && is.infinite(maxtime) && is.infinite(maxsimulation)){
    print("error: you have to provide a finite value for 'maxstep' OR 'maxtime' OR 'maxsimulation'")
    return(NULL)
  }
  target <- results$target
  compute_d <- results$compute_d
  param_algo <- results$param_algo
  #
  thetas_history <- results$thetas_history
  distances_history <- results$distances_history
  threshold_history <- results$threshold_history
  ncomputed <- results$ncomputed
  normcst <- results$normcst
  compute_times <- results$compute_times
  #
  performed_steps <- length(results$thetas_history)
  cat("result file contains the result of", performed_steps, "steps\n")
  # if (is.finite(maxstep)){}
  # nsteps <- performed_steps - 1 + nsteps
  # cat("now aiming for ", param_algo$nsteps, "steps in total\n")
  # otherwise...
  thetas <- results$thetas_history[[performed_steps]]
  distances <- results$distances_history[[performed_steps]]
  latest_y <- results$latest_y
  #
  if (is.null(param_algo$original_proposal)){
    print("no 'original_proposal' in param_algo: something's weird, as the wsmc function should have created it")
    print("current proposal:")
    print(param_algo$proposal)
    print("using this as 'original_proposal' from now on")
    param_algo$original_proposal <- param_algo$proposal
  }
  param_algo$proposal <- param_algo$original_proposal
  #
  step_time <- proc.time()[3]
  total_time <- compute_times[length(compute_times)]
  total_time_continue <- 0
  #
  performed_calculations <- sum(ncomputed)
  istep <- performed_steps
  while (total_time_continue < maxtime && istep < performed_steps+maxstep && sum(ncomputed)-performed_calculations < maxsimulation){
    istep <- istep + 1
    # for (istep in performed_steps:param_algo$nsteps){
    cat("step ", istep, "... running ", sep = "")
    if (is.finite(maxstep))
      cat("for", maxstep, "more steps\n")
    if (is.finite(maxtime))
      cat("until time >=", maxtime, "seconds (since continue)\n")
    if (is.finite(maxsimulation)){
      cat("until # distances calculated >=", maxsimulation, "(since continue)\n")
    }
    one_step_results <- wsmc_one_step(thetas, distances, latest_y, compute_d, target, param_algo)
    #
    thetas <- one_step_results$thetas
    distances <- one_step_results$distances
    latest_y <- one_step_results$latest_y
    param_algo <- one_step_results$param_algo
    ncomputed <- c(ncomputed, one_step_results$ncomputed_thisstep)
    normcst <- c(normcst, one_step_results$normcst)
    cat("  total # distances calculated: ", sum(ncomputed), " (for this step: ", one_step_results$ncomputed_thisstep, ")\n", sep = "")
    # store
    thetas_history[[istep]] <- thetas
    distances_history[[istep]] <- distances
    threshold_history <- c(threshold_history, param_algo$threshold)
    #
    new_step_time <- proc.time()[3]
    total_time <- total_time + new_step_time - step_time
    total_time_continue <- total_time_continue + new_step_time - step_time
    cat("  total time spent:", total_time, "seconds (since continue:", total_time_continue, "s; this step:", new_step_time - step_time, "s)\n")
    step_time <- new_step_time
    compute_times <- c(compute_times, total_time)
    if (!is.null(savefile)){
      results <- list(thetas_history = thetas_history, distances_history = distances_history,
                      threshold_history = threshold_history, param_algo = param_algo,
                      latest_y = latest_y, ncomputed = ncomputed, compute_d = compute_d, target = target,
                      normcst = normcst, compute_times = compute_times)
      save(results, file = savefile)
    }
    if (total_time_continue > maxtime){
      break
    }
  }
  return(list(thetas_history = thetas_history, distances_history = distances_history,
              threshold_history = threshold_history, latest_y = latest_y, param_algo = param_algo, ncomputed = ncomputed,
              compute_d = compute_d, target = target, normcst = normcst, compute_times = compute_times))
}

# takes as argument
# thetas
# distances: associated with each theta
# compute_d, to give to mh_rhit_step
# target, likewise
# param_algo, including a new $threshold and original_proposal
#'@export
wsmc_one_step <- function(thetas, distances, latest_y, compute_d, target, param_algo, debug = FALSE, parallel = TRUE){
  # compute weights
  weights <- as.numeric(distances <= param_algo$threshold)
  normcst <- mean(weights)
  weights <- weights / sum(weights)
  # fit proposal
  param_algo <- update_proposal(param_algo, thetas, weights)
  # systematic resampling
  ancestors <- systematic_resampling(weights)
  thetas <- thetas[ancestors,,drop=F]
  distances <- distances[ancestors]
  old_latest_y <- latest_y
  for (i in 1:nrow(thetas)){
    latest_y[[i]] <- old_latest_y[[ancestors[i]]]
  }
  weights <- rep(1/param_algo$nthetas, param_algo$nthetas)
  # rejuvenation moves
  acceptrates <- rep(0, param_algo$nmoves)
  ncomputed_thisstep <- 0
  for (i in 1:param_algo$nmoves){
    mh_res <- move_step(thetas, distances, latest_y, compute_d, target, param_algo, debug = debug, parallel = parallel)
    thetas <- mh_res$thetas
    distances <- mh_res$distances
    latest_y <- mh_res$latest_y
    ncomputed_thisstep <- ncomputed_thisstep + mh_res$ncomputed
    # re-fit proposal
    param_algo <- update_proposal(param_algo, thetas, weights)
    acceptrates[i] <- mh_res$acceptrate
    cat("  acceptance rates:", 100 * acceptrates[i], "%, threshold =", param_algo$threshold,
        ", min. dist. =", min(distances), "\n")
  }
  #
  nunique <- length(unique(thetas[,1])) # number of unique thetas
  if (debug){
    print("distances before finding new threshold:")
    print(summary(distances))
  }
  # TODO: code robust binary search
  if (nunique/param_algo$nthetas > param_algo$minimum_diversity){
    # if diversity is large enough, we decrease the threshold so that the
    # expected diversity after resampling is still above the threshold
    # that expected diversity is
    one_uniform <- runif(1)
    g <- function(epsilon){
      logw <- log(as.numeric(distances <= epsilon))
      w <- exp(logw - max(logw))
      w <- w / sum(w)
      a <- systematic_resampling_given_u(w, one_uniform)
      th <- thetas[a,1]
      return((length(unique(th))-1)/param_algo$nthetas)
    }
    lower_threshold <- 0
    upper_threshold <- param_algo$threshold
    if (is.infinite(upper_threshold)){
      upper_threshold <- max(distances)
    }
    opt <- try(optimize(f = function(e) (g(e) - param_algo$minimum_diversity)^2, interval = c(lower_threshold, upper_threshold)))
    if (inherits(opt, "try-error")){
      print("error trying to find new threshold: keeping previous threshold")
    } else {
      param_algo$threshold <- opt$minimum
    }
    if (debug){
      cat("new threshold: ", param_algo$threshold, "\n")
    }
  }
  return(list(param_algo = param_algo, thetas = thetas, distances = distances, latest_y = latest_y,
              ncomputed_thisstep = ncomputed_thisstep, normcst = normcst, acceptrates = acceptrates))
}

#'@export
move_step <- function(thetas, distances, latest_y, compute_d, target, param_algo, debug = FALSE, parallel = TRUE){
  res_foreach <- list()
  if (param_algo$R == 0){
    # perform standard MH move
    if (parallel){
      res_foreach <- foreach(i = 1:nrow(thetas)) %dorng% {
        theta <- thetas[i,]
        res <- std_move_step_onetheta(theta, compute_d, target, param_algo)
        res
      }
    } else {
      for (i in 1:nrow(thetas)){
        theta <- thetas[i,]
        res_foreach[[i]] <- std_move_step_onetheta(theta, compute_d, target, param_algo)
      }
    }
  } else {
    # perform R-hit move
    if (parallel){
      res_foreach <- foreach(i = 1:nrow(thetas)) %dorng% {
        theta <- thetas[i,]
        res <- move_step_onetheta(theta, compute_d, target, param_algo)
        res
      }
    } else {
      for (i in 1:nrow(thetas)){
        theta <- thetas[i,]
        res_foreach[[i]] <- move_step_onetheta(theta, compute_d, target, param_algo)
      }
    }
  }
  # thetas_accepted <- sapply(res_foreach, function(x) x$theta)
  ncomputed <- sum(sapply(res_foreach, function(x) x$ncomputed))
  if (debug){
    ncurrents <- sapply(res_foreach, function(x) x$ncurrent)
    nproposals <- sapply(res_foreach, function(x) x$nproposals)
    print("ncurrents:")
    print(summary(ncurrents))
    print("nproposals:")
    print(summary(nproposals))
  }
  accepts <- sapply(res_foreach, function(x) x$accepted)
  for (i in 1:nrow(thetas)){
    if (accepts[i]){
      thetas[i,] <- res_foreach[[i]]$theta
      distances[i] <- res_foreach[[i]]$distance
      latest_y[[i]] <- res_foreach[[i]]$y
    }
  }
  return(list(acceptrate = mean(accepts), thetas = thetas, distances = distances, latest_y = latest_y, ncomputed = ncomputed))
}

# std kernel
#'@export
std_move_step_onetheta <- function(theta, compute_d, target, param_algo){
  threshold <- param_algo$threshold
  proposal <- param_algo$proposal
  theta <- matrix(theta, nrow = 1)
  prior_current <- target$dprior(theta, target$parameters)
  prop_current <- proposal$d(theta, proposal$param_prop)
  theta_prop <- proposal$r(theta, proposal$param_prop)
  prior_proposed <- target$dprior(theta_prop, target$parameters)
  prop_proposed <- proposal$d(theta_prop, proposal$param_prop)
  dproposed <- Inf
  y_prop <- 0
  if (!is.infinite(prior_proposed)){
    y_prop <- target$simulate(theta_prop)
    dproposed <- compute_d(y_prop)
  }
  logratio <- (prior_proposed - prior_current) + (prop_current - prop_proposed)
  logratio <- logratio + log(dproposed < threshold)
  accepted <- (log(runif(1)) < logratio)
  if (accepted){
    return(list(accepted = TRUE, theta = theta_prop, distance = dproposed, y = y_prop,
                nproposals = 1, ncurrent = 0, ncomputed = 1))
  } else {
    return(list(accepted = FALSE, nproposals = 1, ncurrent = 0, ncomputed = 1))
  }
}

# R-hit kernel
#'@export
move_step_onetheta <- function(theta, compute_d, target, param_algo){
  threshold <- param_algo$threshold;   proposal <- param_algo$proposal;  R <- param_algo$R;   maxtrials <- param_algo$maxtrials
  if (is.null(maxtrials)){ maxtrials <- Inf }
  theta <- matrix(theta, nrow = 1)
  #
  prior_current <- target$dprior(theta, target$parameters)
  prop_current <- proposal$d(theta, proposal$param_prop)
  #
  nhits <- 0;  nproposals <- 0;  hit_indices <- c()
  # sample thetas and associated distances until R hits
  thetas_prop <- c();  associated_dists <- c()
  # also compute prior densities as it will be required in case of
  # hit, and also allows to not sample data from parameters outside an admissible range
  associated_priors <- c();  associated_y_fake <- list()
  #
  ncomputed <- 0
  while (nhits < R && nproposals < maxtrials){
    theta_prop <- proposal$r(theta, proposal$param_prop)
    prior_prop <- target$dprior(theta_prop, target$parameters)
    dproposed <- Inf
    y_prop <- 0
    if (!is.infinite(prior_prop)){
      y_prop <- target$simulate(theta_prop)
      dproposed <- compute_d(y_prop)
      ncomputed <- ncomputed + 1
    }
    if (is.na(dproposed)){ dproposed <- Inf }
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
    # store y only if hit, otherwise memory intensive
    if (dproposed < threshold){
      associated_y_fake[[length(associated_y_fake)+1]] <- y_prop
    } else {
      associated_y_fake[[length(associated_y_fake)+1]] <- 0
    }
  }
  # nproposals is N^prime in the notation of the paper
  # choose one of the succesful proposed parameter, except not the last one
  if (nproposals == maxtrials){
    cat("first while loop takes too long, let's reject\n")
    return(list(accepted = FALSE, nproposals = nproposals, ncurrent = NA, ncomputed = ncomputed))
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
    associated_y_fake_L <- associated_y_fake[[L]]
    # now do it again from the selected parameter, until R - 1 hits
    nhits <- 0
    # ncurrent will be N in the notation of the paper
    ncurrent <- 0
    #
    while (nhits < (R-1) && ncurrent < maxtrials){
      theta_prop <- proposal$r(theta_L, proposal$param_prop)
      prior_prop <- target$dprior(theta_prop, target$parameters)
      if (!is.infinite(prior_prop[1])){
        y_prop <- target$simulate(theta_prop)
        dproposed <- compute_d(y_prop)
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
      return(list(accepted = FALSE, nproposals = nproposals, ncurrent = ncurrent, ncomputed = ncomputed))
    } else {
      logratio <- associated_prior_L - (prior_current) +
        (prop_current - proposal$d(theta_L, proposal$param_prop)) +
        log(ncurrent/(nproposals-1))
      accepted <- (log(runif(1)) < logratio)
      if (accepted){
        return(list(accepted = TRUE, theta = theta_L, distance = associated_distance_L, y = associated_y_fake_L,
                    nproposals = nproposals, ncurrent = ncurrent, ncomputed = ncomputed))
      } else {
        return(list(accepted = FALSE, nproposals = nproposals, ncurrent = ncurrent, ncomputed = ncomputed))
      }
    }
  }
}

#'@export
update_proposal <- function(param_algo, thetas, weights){
  # reverse back to original proposal, in case it's been replaced by something else
  param_algo$proposal <- param_algo$original_proposal
  param_algo$proposal$param_prop <- param_algo$proposal$param_update(thetas, weights)
  a <- try(param_algo$proposal$r(thetas[1:2,,drop=F], param_prop = param_algo$proposal$param_prop))
  if (inherits(a, "try-error")){
    print("error in fitting MCMC proposal, trying Rmixmod 5 times..")
    param_algo$proposal <- mixture_rmixmod()
    param_algo$proposal$param_prop <- param_algo$proposal$param_update(thetas, weights)
    a <- try(param_algo$proposal$r(thetas[1:2,,drop=F], param_prop = param_algo$proposal$param_prop))
    nattempts <- 5
    attempt <- 0
    while ((inherits(a, "try-error")) && (attempt < nattempts)){
      attempt <- attempt + 1
      param_algo$proposal$param_prop <- param_algo$proposal$param_update(thetas, weights)
      a <- try(param_algo$proposal$r(thetas[1:2,,drop=F], param_prop = param_algo$proposal$param_prop))
    }
    if (inherits(a, "try-error")){
      print("error in fitting MCMC proposal, switching to Gaussian Mixture with mclust proposal")
      # tmpfilename <- tempfile(fileext = ".RData")
      # save(thetas, weights, param_algo, file = tmpfilename)
      # cat("dumping thetas, weights, param_algo in ", tmpfilename,
      #     "; try executing 'load('", tmpfilename, "'); param_algo$proposal$param_update(thetas, weights)'\n")
      param_algo$proposal <- mixture_mclust()
      param_algo$proposal$param_prop <- param_algo$proposal$param_update(thetas, weights)
      a <- try(param_algo$proposal$r(thetas[1:2,,drop=F], param_prop = param_algo$proposal$param_prop))
      if (inherits(a, "try-error")){
        print("error in fitting MCMC proposal, switching to independent Gaussian")
        # tmpfilename <- tempfile(fileext = ".RData")
        # save(thetas, weights, param_algo, file = tmpfilename)
        # cat("dumping thetas, weights, param_algo in ", tmpfilename,
        # "; try executing 'load('", tmpfilename, "'); param_algo$proposal$param_update(thetas, weights)'\n")
        param_algo$proposal <- independent_proposal()
        param_algo$proposal$param_prop <- param_algo$proposal$param_update(thetas, weights)
        a <- try(param_algo$proposal$r(thetas[1:2,,drop=F], param_prop = param_algo$proposal$param_prop))
        if (inherits(a, "try-error")){
          print("independent Gaussian proposal failed")
          tmpfilename <- "~/tmpfitfail.RData"
          save(thetas, weights, param_algo, file = tmpfilename)
          cat("dumping thetas, weights, param_algo in ", tmpfilename,
              "; try executing 'load('", tmpfilename, "'); param_algo$proposal$param_update(thetas, weights)'\n")
        }
      }
    }
  }
  return(param_algo)
}


