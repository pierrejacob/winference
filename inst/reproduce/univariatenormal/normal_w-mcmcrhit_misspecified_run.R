#+ presets, echo = FALSE, warning = FALSE, message = FALSE
library(winference)
library(ggthemes)
library(doRNG)
registerDoMC(cores = 10)
rm(list = ls())
setmytheme()

set.seed(11)

target <- get_normal()

# number of observations
nobservations <- 10
load(file = "~/Dropbox/ABCD/Results/data/gammadata.RData")
obs <- obs[1:nobservations]

obs_sorted <- sort(obs)
# function to compute distance between observed data and data generated given theta
compute_d <- function(theta){
  fake_rand <- target$generate_randomness(nobservations)
  fake_obs <- target$robservation(nobservations, theta, target$parameters, fake_rand)
  fake_obs_sorted <- sort(fake_obs)
  return(metricL2(obs_sorted, fake_obs_sorted))
}

proposal <- randomwalk_proposal()
proposal$param_prop$cov <- diag(c(1, 1))
nmoves <- 10000
nchains <- 4
# adapts the proposal every 500 iterations
# does not attempt to run the chains in parallel
param_algo <- list(nmoves = nmoves, nchains = nchains, proposal = proposal,
                   R = 5, adaptation = 500, maxtrials = 1000)


#
# ini1 <- initialization(threshold, qloglikelihood, target, param_algo)
# res1 <- quasimh_steps(ini1$thetas, ini1$priorvalues, ini1$qvalues, ini1$distances, threshold, qloglikelihood,  target, param_algo)
# multiple steps of MH r-hit kernel, for multiple chains (not parallelized)

#'@export
rhitmh <- function(threshold, compute_d, target, param_algo){
  nmoves <- param_algo$nmoves
  nchains <- param_algo$nchains
  param_algo$threshold <- threshold
  qloglikelihood <- function(theta, threshold, M){
    d <- compute_d(theta)
    return(list(qvalue = log(d < threshold), distances = d))
  }
  param_algo$M <- 1
  ini <- initialization(threshold, qloglikelihood, target, param_algo, maxattempts = 1000)
  thetas <- ini$thetas
  distances <- ini$distances[,1]

  if (is.null(param_algo$adaptation)) param_algo$adaptation <- 0
  # store histories
  chains_ <- rep(list(matrix(nrow = nmoves, ncol = target$thetadim)), nchains)
  distances_ <- matrix(0, nrow = nmoves, ncol = nchains)
  # priorchains_ <- matrix(nrow = nmoves, ncol = nchains)
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
    # res <- quasimh_step(thetas, priorvalues, qvalues, distances, threshold, qloglikelihood, target, param_algo)
    for (ichain in 1:nchains){
      mh_step_results <- mhrhit_step_onetheta(thetas[ichain,], distances[ichain], compute_d, target, param_algo)
      thetas[ichain,] <- mh_step_results$theta
      distances[ichain] <- mh_step_results$distance
      chains_[[ichain]][imove,] <- thetas[ichain,]
      distances_[imove,ichain] <- distances[ichain]
      accepts <- accepts + mh_step_results$accepted / nchains
    }
  }
  chains_df <- foreach (i = 1:nchains, .combine = rbind) %do% {
    data.frame(ichain = rep(i, nmoves), iteration = 1:nmoves, X = chains_[[i]])
  }
  return(list(chains_list = chains_, chains_df = chains_df, distancechains = distances_, param_algo = param_algo))
}

threshold <- 0.5
filename <- paste0("~/Dropbox/ABCD/Results/univariatenormal/gammadata.n",
                   nobservations, ".wrhit.threshold", threshold, ".RData")
rhitres1 <- rhitmh(threshold, compute_d, target, param_algo)
save(rhitres1, file = filename)

threshold <- 0.25
param_algo <- rhitres1$param_algo
filename <- paste0("~/Dropbox/ABCD/Results/univariatenormal/gammadata.n",
                   nobservations, ".wrhit.threshold", threshold, ".RData")
rhitres2 <- rhitmh(threshold, compute_d, target, param_algo)
save(rhitres2, file = filename)

threshold <- 0.1
param_algo <- rhitres2$param_algo
filename <- paste0("~/Dropbox/ABCD/Results/univariatenormal/gammadata.n",
                   nobservations, ".wrhit.threshold", threshold, ".RData")
rhitres3 <- rhitmh(threshold, compute_d, target, param_algo)
save(rhitres3, file = filename)

# ggplot(res1$chains_df, aes(x = iteration, y = X.1, group = ichain)) + geom_line()
# ggplot(res1$chains_df, aes(x = iteration, y = X.2, group = ichain)) + geom_line()
#
# res1$param_algo

# # lower threshold
#
# filename <- paste0("~/Dropbox/ABCD/Results/univariatenormal/normaldata.n", nobservations, ".metropolis.misspecified.RData")
# load(filename)
# chainlist_to_dataframe <- function(chains_list){
#   nchains <- length(chains_list)
#   niterations <- nrow(chains_list[[1]])
#   chaindf <- foreach (i = 1:nchains, .combine = rbind) %do% {
#     data.frame(ichain = rep(i, niterations), iteration = 1:niterations, X = chains_list[[i]])
#   }
#   return(chaindf)
# }
# mcmc.df <- chainlist_to_dataframe(mh$chains)
#
# mumcmc <- (mcmc.df %>% filter(iteration > 5000))$X.1
# sigmamcmc<- (mcmc.df %>% filter(iteration > 5000))$X.2
#
#
#
# hist(muwmcmc, nclass = 50, prob = TRUE, col = rgb(1,0,0, alpha = 0.5))
# hist(mumcmc, nclass = 50, prob = TRUE, add = TRUE, col = rgb(0,0,0, alpha = 0.5))
#
# hist(sigmawcmc, nclass = 50, prob = TRUE, col = rgb(1,0,0, alpha = 0.5))
# hist(sigmamcmc, nclass = 50, prob = TRUE, add = TRUE, col = rgb(0,0,0, alpha = 0.5))
#
