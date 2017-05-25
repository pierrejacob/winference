#' #'@export
#' wsmc_fixedthresholds <- function(thresholds, compute_d, target, param_algo, savefile = NULL){
#'   nthetas <- param_algo$nthetas
#'   nsteps <- length(thresholds)
#'   nmoves <- param_algo$nmoves
#'   # sample from prior
#'   thetas <- target$rprior(nthetas, target$parameters)
#'   # compute distances
#'   distances <- foreach(itheta = 1:nthetas, .combine = c) %dorng% {
#'     compute_d(thetas[itheta,])
#'   }
#'   #
#'   thetas_history <- list()
#'   thetas_history[[1]] <- thetas
#'   distances_history <- list()
#'   distances_history[[1]] <- distances
#'   #
#'   for (istep in 1:nsteps){
#'     param_algo$threshold <- thresholds[istep]
#'     # compute weights
#'     weights <- as.numeric(distances <= param_algo$threshold)
#'     weights <- weights / sum(weights)
#'     #
#'     param_algo$proposal$param_prop <- param_algo$proposal$param_update(thetas, weights)
#'     # systematic resampling
#'     ancestors <- systematic_resampling(weights)
#'     thetas <- thetas[ancestors,]
#'     distances <- distances[ancestors]
#'     weights <- rep(1/nthetas, nthetas)
#'     # rejuvenation moves
#'     acceptrates <- rep(0, nmoves)
#'     for (i in 1:nmoves){
#'       #
#'       mh_res <- mh_step_within_smc(thetas, distances, compute_d, target, param_algo)
#'       thetas <- mh_res$thetas
#'       distances <- mh_res$distances
#'       param_algo$proposal$param_prop <- param_algo$proposal$param_update(thetas, weights)
#'       acceptrates[i] <- mh_res$acceptrate
#'     }
#'     cat("step", istep, "/", nsteps, ", ")
#'     cat(" acceptance rates:", 100 * acceptrates, "%, threshold =", param_algo$threshold,
#'         ", min. dist. =", min(distances), "\n")
#'     #
#'     thetas_history[[istep+1]] <- thetas
#'     distances_history[[istep+1]] <- distances
#'     #
#'     if (!is.null(savefile)){
#'       results <- list(thetas_history = thetas_history, distances_history = distances_history, param_algo = param_algo)
#'       save(results, file = savefile)
#'     }
#'   }
#'   return(list(thetas_history = thetas_history, distances_history = distances_history, param_algo = param_algo))
#' }
#'
#' #'@export
#' wsmc <- function(compute_d, target, param_algo, savefile = NULL){
#'   nthetas <- param_algo$nthetas
#'   nsteps <- param_algo$nsteps
#'   nmoves <- param_algo$nmoves
#'   minimum_diversity <- param_algo$minimum_diversity
#'   if (is.null(minimum_diversity)){
#'     minimum_diversity <- 0.5
#'   }
#'   # sample from prior
#'   thetas <- target$rprior(nthetas, target$parameters)
#'   # compute distances
#'   distances <- foreach(itheta = 1:nthetas, .combine = c) %dorng% {
#'     compute_d(thetas[itheta,])
#'   }
#'   #
#'   param_algo$threshold <- as.numeric(quantile(distances, probs = minimum_diversity))
#'   threshold_history <- param_algo$threshold
#'   #
#'   thetas_history <- list()
#'   thetas_history[[1]] <- thetas
#'   distances_history <- list()
#'   distances_history[[1]] <- distances
#'   #
#'   for (istep in 1:nsteps){
#'     # compute weights
#'     weights <- as.numeric(distances <= param_algo$threshold)
#'     weights <- weights / sum(weights)
#'     #
#'     param_algo$proposal$param_prop <- param_algo$proposal$param_update(thetas, weights)
#'     # systematic resampling
#'     ancestors <- systematic_resampling(weights)
#'     thetas <- thetas[ancestors,]
#'     distances <- distances[ancestors]
#'     weights <- rep(1/nthetas, nthetas)
#'     # rejuvenation moves
#'     acceptrates <- rep(0, nmoves)
#'     for (i in 1:nmoves){
#'       #
#'       mh_res <- mh_step_within_smc(thetas, distances, compute_d, target, param_algo)
#'       thetas <- mh_res$thetas
#'       distances <- mh_res$distances
#'       param_algo$proposal$param_prop <- param_algo$proposal$param_update(thetas, weights)
#'       acceptrates[i] <- mh_res$acceptrate
#'     }
#'     cat("step", istep, "/", nsteps, ", ")
#'     cat(" acceptance rates:", 100 * acceptrates, "%, threshold =", param_algo$threshold,
#'         ", min. dist. =", min(distances), "\n")
#'     #
#'     nunique <- length(unique(thetas[,1])) # number of unique thetas
#'     #
#'     if (nunique/nthetas > minimum_diversity){
#'       # if diversity is large enough, we decrease the threshold so that the
#'       # expected diversity after resampling is still above the threshold
#'       # that expected diversity is
#'       one_uniform <- runif(1)
#'       g <- function(epsilon){
#'         logw <- log(as.numeric(distances <= epsilon))
#'         w <- exp(logw - max(logw))
#'         w <- w / sum(w)
#'         a <- systematic_resampling_given_u(w, one_uniform)
#'         th <- thetas[a,1]
#'         return((length(unique(th))-1)/nthetas)
#'       }
#'       opt <- optimize(f = function(e) (g(e) - minimum_diversity)^2, interval = c(0, param_algo$threshold))
#'       param_algo$threshold <- opt$minimum
#'     }
#'     thetas_history[[istep+1]] <- thetas
#'     distances_history[[istep+1]] <- distances
#'     threshold_history <- c(threshold_history, param_algo$threshold)
#'     #
#'     if (!is.null(savefile)){
#'       results <- list(thetas_history = thetas_history, distances_history = distances_history,
#'                       threshold_history = threshold_history,
#'                       param_algo = param_algo)
#'       save(results, file = savefile)
#'     }
#'   }
#'   return(list(thetas_history = thetas_history, distances_history = distances_history,
#'               threshold_history = threshold_history,
#'               param_algo = param_algo))
#' }
#'
#'
#' #'@export
#' mh_step_within_smc <- function(thetas, distances, compute_d, target, param_algo){
#'   threshold <- param_algo$threshold
#'   proposal <- param_algo$proposal
#'   nthetas <- nrow(thetas)
#'   #
#'   thetas_prop <- proposal$r(thetas, proposal$param_prop)
#'   prior_proposed <- target$dprior(thetas_prop, target$parameters)
#'   prop_proposed <- proposal$d(thetas_prop, proposal$param_prop)
#'   #
#'   prior_current <- target$dprior(thetas, target$parameters)
#'   prop_current <- proposal$d(thetas, proposal$param_prop)
#'   #
#'   logratios <- (prior_proposed - prior_current) + (prop_current - prop_proposed)
#'   dist_proposed <- rep(Inf, nthetas)
#'   potential_indices <- which(!is.infinite(logratios))
#'   res_foreach <- foreach(ipot = potential_indices, .combine = rbind) %dorng% {
#'     matrix(compute_d(thetas_prop[ipot,]), nrow = 1)
#'   }
#'   dist_proposed[potential_indices] <- as.numeric(res_foreach[,1])
#'   logratios <- logratios + log(dist_proposed < threshold)
#'   accepted <- (log(runif(nthetas)) < logratios)
#'   #
#'   thetas[accepted] <- thetas_prop[accepted]
#'   distances[accepted] <- dist_proposed[accepted]
#'   return(list(acceptrate = mean(accepted), thetas = thetas, distances = distances))
#' }
#'
#'
#' #'@rdname wsmc_rhit
#' #'@title Adaptive SMC sampler with r-hit kernels
#' #'@description Runs an adaptive SMC sampler with r-hit kernels
#' #'
#' #' References for the r-hit kernels:
#' #' \itemize{
#' #' \item Lee, A. (2012). On the choice of MCMC kernels for approximate Bayesian computation with SMC samplers.
#' #'  In Proceedings of the 2012 Winter Simulation Conference, pages 304–315.
#' #' \item Lee, A. and Łatuszyński, K. (2014). Variance bounding and geometric ergodicity of Markov chain Monte Carlo kernels for approximate Bayesian computation. Biometrika, 101(3):655–671.
#' #'}
#' #'@return a list containing
#' #' \itemize{
#' #' \item thetas_history (all theta particles at all steps),
#' #'  \item distances_history, threshold_history, param_algo
#' #'  }
#' #'@export
#' wsmc_rhit <- function(compute_d, target, param_algo, savefile = NULL){
#'   original_proposal <- param_algo$proposal
#'   nthetas <- param_algo$nthetas
#'   nsteps <- param_algo$nsteps
#'   nmoves <- param_algo$nmoves
#'   minimum_diversity <- param_algo$minimum_diversity
#'   if (is.null(minimum_diversity)){
#'     minimum_diversity <- 0.5
#'   }
#'   # sample from prior
#'   thetas <- target$rprior(nthetas, target$parameters)
#'   # compute distances
#'   distances <- foreach(itheta = 1:nthetas, .combine = c) %dorng% {
#'     compute_d(thetas[itheta,])
#'   }
#'   # count of the number of distance calculations
#'   ncomputed <- c(nthetas)
#'   #
#'   param_algo$threshold <- as.numeric(quantile(distances, probs = minimum_diversity))
#'   threshold_history <- param_algo$threshold
#'   #
#'   thetas_history <- list()
#'   thetas_history[[1]] <- thetas
#'   distances_history <- list()
#'   distances_history[[1]] <- distances
#'   #
#'   for (istep in 1:nsteps){
#'     # compute weights
#'     weights <- as.numeric(distances <= param_algo$threshold)
#'     weights <- weights / sum(weights)
#'     # reverse back to original proposal, in case it's been replaced by something else
#'     param_algo$proposal <- original_proposal
#'     param_algo$proposal$param_prop <- param_algo$proposal$param_update(thetas, weights)
#'     a <- try(param_algo$proposal$r(thetas[1:2,,drop=F], param_prop = param_algo$proposal$param_prop))
#'     if (inherits(a, "try-error")){
#'       print("error in the MCMC proposal, switching to independent Gaussian proposal")
#'       iprop <- independent_proposal()
#'       param_algo$proposal <- iprop
#'       param_algo$proposal$param_prop <- param_algo$proposal$param_update(thetas, weights)
#'     }
#'     # systematic resampling
#'     ancestors <- systematic_resampling(weights)
#'     thetas <- thetas[ancestors,]
#'     distances <- distances[ancestors]
#'     weights <- rep(1/nthetas, nthetas)
#'     # rejuvenation moves
#'     acceptrates <- rep(0, nmoves)
#'     ncomputed_thisstep <- 0
#'     for (i in 1:nmoves){
#'       #
#'       mh_res <- mh_rhit_step(thetas, distances, compute_d, target, param_algo)
#'       thetas <- mh_res$thetas
#'       distances <- mh_res$distances
#'       ncomputed_thisstep <- ncomputed_thisstep + mh_res$ncomputed
#'       param_algo$proposal$param_prop <- param_algo$proposal$param_update(thetas, weights)
#'       acceptrates[i] <- mh_res$acceptrate
#'     }
#'     ncomputed <- c(ncomputed, ncomputed_thisstep)
#'     cat("step", istep, "/", nsteps, ", ")
#'     cat(" acceptance rates:", 100 * acceptrates, "%, threshold =", param_algo$threshold,
#'         ", min. dist. =", min(distances), "\n")
#'     cat("number of distance calculations at this step: ", ncomputed_thisstep, "\n")
#'     #
#'     nunique <- length(unique(thetas[,1])) # number of unique thetas
#'     #
#'     if (nunique/nthetas > minimum_diversity){
#'       # if diversity is large enough, we decrease the threshold so that the
#'       # expected diversity after resampling is still above the threshold
#'       # that expected diversity is
#'       one_uniform <- runif(1)
#'       g <- function(epsilon){
#'         logw <- log(as.numeric(distances <= epsilon))
#'         w <- exp(logw - max(logw))
#'         w <- w / sum(w)
#'         a <- systematic_resampling_given_u(w, one_uniform)
#'         th <- thetas[a,1]
#'         return((length(unique(th))-1)/nthetas)
#'       }
#'       opt <- optimize(f = function(e) (g(e) - minimum_diversity)^2, interval = c(0, param_algo$threshold))
#'       param_algo$threshold <- opt$minimum
#'     }
#'     thetas_history[[istep+1]] <- thetas
#'     distances_history[[istep+1]] <- distances
#'     threshold_history <- c(threshold_history, param_algo$threshold)
#'     #
#'     if (!is.null(savefile)){
#'       results <- list(thetas_history = thetas_history, distances_history = distances_history,
#'                       threshold_history = threshold_history,
#'                       param_algo = param_algo, ncomputed = ncomputed)
#'       save(results, file = savefile)
#'     }
#'   }
#'   return(list(thetas_history = thetas_history, distances_history = distances_history,
#'               threshold_history = threshold_history,
#'               param_algo = param_algo, ncomputed = ncomputed))
#' }
#'
#' #'@export
#' wsmc_rhit_continue <- function(results, compute_d, target, savefile = NULL){
#'   threshold_history <- results$threshold_history
#'   thetas_history <- results$threshold_history
#'   distances_history <- results$distances_history
#'   param_algo <- results$param_algo
#'   original_proposal <- param_algo$proposal
#'   #
#'   nsteps <- param_algo$nsteps
#'   donesteps <- length(thetas_history)
#'   thetas <- results$thetas_history[[donesteps]]
#'   distances <- results$distances_history[[donesteps]]
#'   #
#'   nthetas <- param_algo$nthetas
#'   nmoves <- param_algo$nmoves
#'   minimum_diversity <- param_algo$minimum_diversity
#'   #
#'   #
#'   for (istep in (donesteps+1):nsteps){
#'     # compute weights
#'     weights <- as.numeric(distances <= param_algo$threshold)
#'     weights <- weights / sum(weights)
#'     #
#'     param_algo$proposal <- original_proposal
#'     param_algo$proposal$param_prop <- param_algo$proposal$param_update(thetas, weights)
#'     a <- try(param_algo$proposal$r(thetas[1:2,], param_prop = param_algo$proposal$param_prop))
#'     if (inherits(a, "try-error")){
#'       iprop <- randomwalk_proposal()
#'       param_algo$proposal <- iprop
#'       param_algo$proposal$param_prop <- param_algo$proposal$param_update(thetas, weights)
#'     }
#'     # systematic resampling
#'     ancestors <- systematic_resampling(weights)
#'     thetas <- thetas[ancestors,]
#'     distances <- distances[ancestors]
#'     weights <- rep(1/nthetas, nthetas)
#'     # rejuvenation moves
#'     acceptrates <- rep(0, nmoves)
#'     for (i in 1:nmoves){
#'       #
#'       mh_res <- mh_rhit_step(thetas, distances, compute_d, target, param_algo)
#'       thetas <- mh_res$thetas
#'       distances <- mh_res$distances
#'       param_algo$proposal$param_prop <- param_algo$proposal$param_update(thetas, weights)
#'       acceptrates[i] <- mh_res$acceptrate
#'     }
#'     cat("step", istep, "/", nsteps, ", ")
#'     cat(" acceptance rates:", 100 * acceptrates, "%, threshold =", param_algo$threshold,
#'         ", min. dist. =", min(distances), "\n")
#'     #
#'     nunique <- length(unique(thetas[,1])) # number of unique thetas
#'     #
#'     if (nunique/nthetas > minimum_diversity){
#'       # if diversity is large enough, we decrease the threshold so that the
#'       # expected diversity after resampling is still above the threshold
#'       # that expected diversity is
#'       one_uniform <- runif(1)
#'       g <- function(epsilon){
#'         logw <- log(as.numeric(distances <= epsilon))
#'         w <- exp(logw - max(logw))
#'         w <- w / sum(w)
#'         a <- systematic_resampling_given_u(w, one_uniform)
#'         th <- thetas[a,1]
#'         return((length(unique(th))-1)/nthetas)
#'       }
#'       opt <- optimize(f = function(e) (g(e) - minimum_diversity)^2, interval = c(0, param_algo$threshold))
#'       param_algo$threshold <- opt$minimum
#'     }
#'     thetas_history[[length(thetas_history)+1]] <- thetas
#'     distances_history[[length(thetas_history)+1]] <- distances
#'     threshold_history <- c(threshold_history, param_algo$threshold)
#'     #
#'     if (!is.null(savefile)){
#'       results <- list(thetas_history = thetas_history, distances_history = distances_history,
#'                       threshold_history = threshold_history,
#'                       param_algo = param_algo)
#'       save(results, file = savefile)
#'     }
#'   }
#'   return(list(thetas_history = thetas_history, distances_history = distances_history,
#'               threshold_history = threshold_history,
#'               param_algo = param_algo))
#' }
