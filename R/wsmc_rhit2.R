#'
#' #'@rdname wsmc_rhit2
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
#' wsmc_rhit2 <- function(compute_d, target, param_algo, savefile = NULL){
#'   param_algo$original_proposal <- param_algo$proposal
#'   if (is.null(param_algo$minimum_diversity)){
#'     param_algo$minimum_diversity <- 0.5
#'   }
#'   # sample from prior
#'   thetas <- target$rprior(param_algo$nthetas, target$parameters)
#'   # compute distances
#'   distances <- foreach(itheta = 1:param_algo$nthetas, .combine = c) %dorng% {
#'     compute_d(thetas[itheta,])
#'   }
#'   distances[is.na(distances)] <- Inf
#'   # count of the number of distance calculations
#'   ncomputed <- c(param_algo$nthetas)
#'   normcst <- c()
#'   #
#'   param_algo$threshold <- as.numeric(quantile(distances, probs = param_algo$minimum_diversity))
#'   threshold_history <- param_algo$threshold
#'   #
#'   thetas_history <- list()
#'   thetas_history[[1]] <- thetas
#'   distances_history <- list()
#'   distances_history[[1]] <- distances
#'   #
#'   for (istep in 1:param_algo$nsteps){
#'     cat("step", istep, "/", param_algo$nsteps, "\n")
#'     one_step_results <- wsmc_rhit_one_step_(thetas, distances, compute_d, target, param_algo)
#'     #
#'     thetas <- one_step_results$thetas
#'     distances <- one_step_results$distances
#'     param_algo <- one_step_results$param_algo
#'     ncomputed <- c(ncomputed, one_step_results$ncomputed_thisstep)
#'     normcst <- c(normcst, one_step_results$normcst)
#'     # store
#'     thetas_history[[istep+1]] <- thetas
#'     distances_history[[istep+1]] <- distances
#'     threshold_history <- c(threshold_history, param_algo$threshold)
#'     #
#'     if (!is.null(savefile)){
#'       results <- list(thetas_history = thetas_history, distances_history = distances_history,
#'                       threshold_history = threshold_history,
#'                       ncomputed = ncomputed, param_algo = param_algo, compute_d = compute_d,
#'                       target = target, normcst = normcst)
#'       save(results, file = savefile)
#'     }
#'   }
#'   return(list(thetas_history = thetas_history, distances_history = distances_history,
#'               threshold_history = threshold_history, param_algo = param_algo,
#'               ncomputed = ncomputed, compute_d = compute_d, target = target, normcst = normcst))
#' }
#'
#' #'@export
#' wsmc_rhit2_continue <- function(results, nsteps, savefile = NULL){
#'   if (missing(nsteps) || nsteps <= 0){
#'     cat("Need to specify nsteps, an integer >= 1.\n")
#'     cat("It represents the number of new desired steps. Exiting now.\n")
#'     return(NULL)
#'   }
#'   target <- results$target
#'   compute_d <- results$compute_d
#'   param_algo <- results$param_algo
#'   #
#'   thetas_history <- results$thetas_history
#'   distances_history <- results$distances_history
#'   threshold_history <- results$threshold_history
#'   ncomputed <- results$ncomputed
#'   normcst <- results$normcst
#'   #
#'   performed_steps <- length(results$thetas_history)
#'   cat("result file contains the result of ", performed_steps-1, "steps\n")
#'   param_algo$nsteps <- performed_steps - 1 + nsteps
#'   cat("now aiming for ", param_algo$nsteps, "steps in total\n")
#'   # otherwise...
#'   thetas <- results$thetas_history[[performed_steps]]
#'   distances <- results$distances_history[[performed_steps]]
#'   #
#'   param_algo$original_proposal <- param_algo$proposal
#'   #
#'   for (istep in performed_steps:param_algo$nsteps){
#'     cat("step", istep, "/", param_algo$nsteps, "\n")
#'     one_step_results <- wsmc_rhit_one_step_(thetas, distances, compute_d, target, param_algo)
#'     #
#'     thetas <- one_step_results$thetas
#'     distances <- one_step_results$distances
#'     param_algo <- one_step_results$param_algo
#'     ncomputed <- c(ncomputed, one_step_results$ncomputed_thisstep)
#'     normcst <- c(normcst, one_step_results$normcst)
#'     # store
#'     thetas_history[[istep+1]] <- thetas
#'     distances_history[[istep+1]] <- distances
#'     threshold_history <- c(threshold_history, param_algo$threshold)
#'     #
#'     if (!is.null(savefile)){
#'       results <- list(thetas_history = thetas_history, distances_history = distances_history,
#'                       threshold_history = threshold_history, param_algo = param_algo,
#'                       ncomputed = ncomputed, compute_d = compute_d, target = target,
#'                       normcst = normcst)
#'       save(results, file = savefile)
#'     }
#'   }
#'   return(list(thetas_history = thetas_history, distances_history = distances_history,
#'               threshold_history = threshold_history, param_algo = param_algo, ncomputed = ncomputed,
#'               compute_d = compute_d, target = target, normcst = normcst))
#' }
#'
#'
#' # takes as argument
#' # thetas
#' # distances: associated with each theta
#' # compute_d, to give to mh_rhit_step
#' # target, likewise
#' # param_algo, including a new $threshold and original_proposal
#' wsmc_rhit_one_step_ <- function(thetas, distances, compute_d, target, param_algo){
#'   # compute weights
#'   weights <- as.numeric(distances <= param_algo$threshold)
#'   normcst <- mean(weights)
#'   weights <- weights / sum(weights)
#'   # fit proposal
#'   param_algo <- update_proposal(param_algo, thetas, weights)
#'   # systematic resampling
#'   ancestors <- systematic_resampling(weights)
#'   thetas <- thetas[ancestors,]
#'   distances <- distances[ancestors]
#'   weights <- rep(1/param_algo$nthetas, param_algo$nthetas)
#'   # rejuvenation moves
#'   acceptrates <- rep(0, param_algo$nmoves)
#'   ncomputed_thisstep <- 0
#'   for (i in 1:param_algo$nmoves){
#'     mh_res <- mh_rhit_step(thetas, distances, compute_d, target, param_algo)
#'     thetas <- mh_res$thetas
#'     distances <- mh_res$distances
#'     ncomputed_thisstep <- ncomputed_thisstep + mh_res$ncomputed
#'     # re-fit proposal
#'     param_algo <- update_proposal(param_algo, thetas, weights)
#'     acceptrates[i] <- mh_res$acceptrate
#'   }
#'   # ncomputed <- c(ncomputed, ncomputed_thisstep)
#'   cat(" acceptance rates:", 100 * acceptrates, "%, threshold =", param_algo$threshold,
#'       ", min. dist. =", min(distances), "\n")
#'   cat("number of distance calculations at this step: ", ncomputed_thisstep, "\n")
#'   #
#'   nunique <- length(unique(thetas[,1])) # number of unique thetas
#'   # TODO: code robust binary search
#'   if (nunique/param_algo$nthetas > param_algo$minimum_diversity){
#'     # if diversity is large enough, we decrease the threshold so that the
#'     # expected diversity after resampling is still above the threshold
#'     # that expected diversity is
#'     one_uniform <- runif(1)
#'     g <- function(epsilon){
#'       logw <- log(as.numeric(distances <= epsilon))
#'       w <- exp(logw - max(logw))
#'       w <- w / sum(w)
#'       a <- systematic_resampling_given_u(w, one_uniform)
#'       th <- thetas[a,1]
#'       return((length(unique(th))-1)/param_algo$nthetas)
#'     }
#'     opt <- optimize(f = function(e) (g(e) - param_algo$minimum_diversity)^2, interval = c(0, param_algo$threshold))
#'     param_algo$threshold <- opt$minimum
#'   }
#'   return(list(param_algo = param_algo, thetas = thetas, distances = distances,
#'               ncomputed_thisstep = ncomputed_thisstep, normcst = normcst))
#' }
