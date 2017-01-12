#+ presets, echo = FALSE, warning = FALSE, message = FALSE
library(winference)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(foreach)
library(doMC)
library(doRNG)
library(reshape2)
registerDoMC(cores = 10)
rm(list = ls())
setmytheme()
set.seed(11)
target <- get_normal()
# number of observations
nobservations <- 10
load(file = "~/Dropbox/ABCD/Results/data/normaldata.RData")
obs <- obs[1:nobservations]
obs_sorted <- sort(obs)

# function to compute distance between observed data and data generated given theta
compute_d <- function(theta, metric = metricL2){
  fake_rand <- target$generate_randomness(nobservations)
  fake_obs <- target$robservation(nobservations, theta, target$parameters, fake_rand)
  fake_obs_sorted <- sort(fake_obs)
  return(metric(obs_sorted, fake_obs_sorted))
}

proposal <- independent_proposal()

param_algo <- list(nthetas = 1024, nmoves = 1, proposal = proposal,
                   nsteps = 10, minimum_diversity = 0.5)

# param_algo$hit <- FALSE
# param_algo$R <- 2
# param_algo$maxtrials <- 1000

smc_step <- function(thetas, qvalues, qloglikelihood, target, param_algo){
  threshold <- param_algo$threshold
  nmoves <- param_algo$nmoves
  # compute weights
  weights <- as.numeric(qvalues[,1])
  weights <- weights / sum(weights)
  #
  param_algo$proposal$param_prop <- param_algo$proposal$param_update(thetas, weights)
  # systematic resampling
  ancestors <- systematic_resampling(weights)
  thetas <- thetas[ancestors,]
  priorvalues <- priorvalues[ancestors]
  qvalues <- qvalues[ancestors,]

  weights <- rep(1/nrow(thetas), nrow(thetas))
  # rejuvenation moves
  acceptrates <- rep(0, nmoves)
  for (i in 1:nmoves){
    #
    mh_res <- mh_step(thetas, priorvalues, qvalues, threshold = threshold, qloglikelihood, target, param_algo)
    thetas <- mh_res$thetas
    qvalues <- mh_res$qvalues
    priorvalues <- mh_res$priorvalues
    param_algo$proposal$param_prop <- param_algo$proposal$param_update(thetas, weights)
    acceptrates[i] <- mh_res$acceptrate
  }
  return(list(thetas = thetas, distances = distances, acceptrates = acceptrates,
              param_algo = param_algo))
}

smcsampler <- function(qloglikelihood, target, param_algo, savefile = NULL){
  nthetas <- param_algo$nthetas
  nsteps <- param_algo$nsteps
  if (is.null(param_algo$M)) param_algo$M <- 1
  minimum_diversity <- param_algo$minimum_diversity
  if (is.null(minimum_diversity)){
    minimum_diversity <- 0.5
  }
  # default to not using r-hit kernels
  if (is.null(param_algo$hit)){
    param_algo$hit <- FALSE
  }
  # sample from prior
  thetas <- target$rprior(nthetas, target$parameters)
  # compute W-distances
  qvalues <- matrix(foreach(itheta = 1:nthetas, .combine = rbind) %dorng% {
    qloglikelihood(thetas[itheta,], Inf, param_algo$M)
  }, nrow = nthetas)
  distances <- as.numeric(qvalues[,2:ncol(qvalues)])
  #
  threshold <- as.numeric(quantile(distances, probs = minimum_diversity))
  threshold_history <- threshold
  #
  thetas_history <- list()
  thetas_history[[1]] <- thetas
  qs_history <- list()
  qs_history[[1]] <- qvalues[,1]
  #
  for (istep in 1:nsteps){
    res <- smc_step(thetas, qvalues, threshold, qloglikelihood, target, param_algo)
    param_algo <- res$param_algo
    cat("step", istep, "/", nsteps, ", ")
    cat(" acceptance rates:", 100 * res$acceptrates, "%, threshold =", param_algo$threshold,
        ", min. dist. =", min(distances), "\n")
    #
    thetas <- res$thetas
    distances <- res$qvalues
    #
    nunique <- length(unique(thetas[,1])) # number of unique thetas
    #
    if (nunique/nthetas > minimum_diversity){
      # if diversity is large enough, we decrease the threshold so that the
      # expected diversity after resampling is still above the threshold
      # that expected diversity is
      g <- function(epsilon){
        weights <- as.numeric(distances <= epsilon)
        return((length(unique(thetas[,1] * weights))-1)/nrow(thetas))
      }
      opt <- optimize(f = function(e) (g(e) - minimum_diversity)^2, interval = c(0, param_algo$threshold))
      threshold <- opt$minimum
    }
    thetas_history[[istep+1]] <- thetas
    qs_history[[istep+1]] <- qs[,1]
    threshold_history <- c(threshold_history, threshold)
    #
    if (!is.null(savefile)){
      results <- list(thetas_history = thetas_history, distances_history = distances_history,
                      threshold_history = threshold_history,
                      param_algo = param_algo)
      save(results, file = savefile)
    }
  }
  return(list(thetas_history = thetas_history, distances_history = distances_history,
              threshold_history = threshold_history,
              param_algo = param_algo))
}



# ############################################################################################################################
# ############################################################################################################################
# filename <- paste0("~/Dropbox/ABCD/Results/univariatenormal/normaldata.n", nobservations, ".wsmc.RData")
# # wsmcresults <- smcsampler(compute_d, target, param_algo, savefile = filename)
# # save(wsmcresults, nobservations, file = filename)
# load(filename)
# wsmcresults <- results
#
#
# ## check with rejection sampler
# # m <- wsmcresults$param_algo$proposal$param_prop$mean
# # cov <- wsmcresults$param_algo$proposal$param_prop$cov
# #
# # epsilon <- wsmcresults$threshold_history[length(wsmcresults$threshold_history)]
# # epsilon <- 0.005
# #
# # theta <- fast_rmvnorm(1, m, cov)
# # compute_dm <- function(theta, metric = metricL2){
# #   m <- 10
# #   ds <- rep(0, m)
# #   for (i in 1:m){
# #     fake_rand <- target$generate_randomness(nobservations)
# #     fake_obs <- target$robservation(nobservations, theta, target$parameters, fake_rand)
# #     fake_obs_sorted <- sort(fake_obs)
# #     ds[i] <- metric(obs_sorted, fake_obs_sorted)
# #   }
# #   return(mean(ds))
# # }
# # d <- compute_dm(theta[1,])
# # d
#
#
# gthreshold <- qplot(x = 1:(length(wsmcresults$threshold_history)), y = wsmcresults$threshold_history, geom = "line")
# gthreshold <- gthreshold + xlab("step") + ylab("threshold")
# gthreshold
# tail(wsmcresults$threshold_history)
#
# wsmc_to_dataframe <- function(results, parameter_names){
#   th <- results$thetas_history
#   nsteps <- length(th)
#   df <- data.frame()
#   for (i in 1:(nsteps)){
#     df_ <- data.frame(cbind(th[[i]], rep(i, nrow(th[[i]]))))
#     names(df_) <- c(parameter_names, "step")
#     df <- rbind(df, df_)
#   }
#   names(df) <- c(parameter_names, "step")
#   return(df)
# }
#
# df <- wsmc_to_dataframe(wsmcresults, target$parameter_names)
# nsteps <- max(df$step)
# head(df)
#
# g <- ggplot(df, aes(x = mu, y = sigma, colour = step, group = step))
# g <- g + geom_point(alpha = 0.5)
# g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
# g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
# g
#
# g <- ggplot(df, aes(x = mu, colour = factor(step), group = step)) + geom_density(aes(y = ..density..))
# g <- g + geom_vline(xintercept = true_theta[1])
# g
#
#
# filename <- paste0("~/Dropbox/ABCD/Results/univariatenormal/normaldata.n", nobservations, ".metropolis.RData")
# load(filename)
#
# chainlist_to_dataframe <- function(chains_list){
#   nchains <- length(chains_list)
#   niterations <- nrow(chains_list[[1]])
#   chaindf <- foreach (i = 1:nchains, .combine = rbind) %do% {
#     data.frame(ichain = rep(i, niterations), iteration = 1:niterations, X = chains_list[[i]])
#   }
#   return(chaindf)
# }
# chaindf <- chainlist_to_dataframe(mh$chains)
#
# qqplot((chaindf %>% filter(iteration > 5000))$X.1, (df %>% filter(step == nsteps))$mu)
# mumcmc <- (chaindf %>% filter(iteration > 5000))$X.1
# sigmamcmc<- (chaindf %>% filter(iteration > 5000))$X.2
# # hist(mumcmc)
#
# g <- ggplot(df %>% filter(step %in% c(8, 9, 10, nsteps)), aes(x = mu, fill = factor(step), group = step))
# g <- g + geom_density(data=chaindf %>% filter(iteration > 5000), aes(x = X.1, y = ..density.., fill = NULL, group = NULL), fill = "black")
# g <- g + geom_density(aes(y = ..density..), alpha = 0.5)
# g <- g + geom_density(data=chaindf %>% filter(iteration > 5000), aes(x = X.1, y = ..density.., fill = NULL, group = NULL), colour = "black", size = 2)
# g <- g + geom_vline(xintercept = true_theta[1])
# g
#
# g <- ggplot(df %>% filter(step %in% c(8, 9, 10, 15, nsteps)), aes(x = sigma, fill = factor(step), group = step))
# g <- g + geom_density(data=chaindf %>% filter(iteration > 5000), aes(x = X.2, y = ..density.., fill = NULL, group = NULL), fill = "black")
# g <- g + geom_density(aes(y = ..density..), alpha = 0.5)
# g <- g + geom_density(data=chaindf %>% filter(iteration > 5000), aes(x = X.2, y = ..density.., fill = NULL, group = NULL), colour = "black", size = 2)
# g <- g + geom_vline(xintercept = true_theta[2])
# g
#
# ## produce more samples from the latest distributions
# thetas <- wsmcresults$thetas_history[[length(wsmcresults$thetas_history)]]
# distances <- wsmcresults$distances_history[[length(wsmcresults$distances_history)]]
# hist(thetas[,1], nclass = 20, prob = TRUE, col = rgb(1,0,0, alpha = 0.5))
# hist(mumcmc, nclass = 50, prob = TRUE, add = TRUE, col = rgb(0,0,0, alpha = 0.5))
# hist(thetas[,2], nclass = 50, prob = TRUE, col = rgb(1,0,0, alpha = 0.5))
# hist(sigmamcmc, nclass = 50, prob = TRUE, add = TRUE, col = rgb(0,0,0, alpha = 0.5))
# abline(v = sd(obs))
#
# hist(thetas[,1], nclass = 100)
#
# param_algo <- wsmcresults$param_algo
# thetas_list <- list()
# distances_list <- list()
# thetas_list[[1]] <- thetas
# distances_list[[1]] <- distances
# nmoves <- 10
# for (i in 1:nmoves){
#   step_results <- mh_rhit_step(thetas, distances, compute_d, target, param_algo)
#   # step_results <- mh_step(thetas, distances, compute_d, target, param_algo)
#   print(step_results$acceptrate)
#   thetas <- step_results$thetas
#   distances <- step_results$distances
#   thetas_list[[i+1]] <- thetas
#   distances_list[[i+1]] <- distances
# }
# filename <- paste0("~/Dropbox/ABCD/Results/univariatenormal/normaldata.n", nobservations, ".wsmc.rejuvenated.RData")
# save(thetas_list, distances_list, file = filename)
# load(filename)
# nmoves <- length(thetas_list)
# #
# hist(thetas_list[[nmoves]][,1], nclass = 50, prob = TRUE)
# ra <- range(thetas_list[[nmoves]][,1])
# lines(density(mumcmc[mumcmc > ra[1] && mumcmc < ra[2]]))
#
#
# hist(thetas_list[[nmoves]][,2], nclass = 50, prob = TRUE)
# ra <- range(thetas_list[[nmoves]][,2])
# lines(density(sigmamcmc[sigmamcmc > ra[1] && sigmamcmc < ra[2]]))
