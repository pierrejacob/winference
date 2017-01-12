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

# target parameters
target$parameters <- list(mu_0 = 0, nu = 1, alpha = 2, beta = 1)

# number of observations
nobservations <- 10000
load(file = "~/Dropbox/ABCD/Results/data/normaldata.RData")
obs <- obs[1:nobservations]

# Gaussian log-likelihood for one observation, and a matrix of particles
loglikelihood <- function(thetaparticles, observation){
  logdensities <- dnorm(observation, mean = thetaparticles[,1],
                        sd = thetaparticles[,2], log = TRUE)
  return(logdensities)
}

# Full log-likelihood
full_loglikelihood <- function(thetaparticles, observations){
  ll <- 0
  for (i in 1:length(observations)){
    ll <- ll + loglikelihood(thetaparticles, observations[i])
  }
  return(ll)
}

# Posterior density function (log)
posterior <- function(thetaparticles){
  logdens <- target$dprior(thetaparticles, target$parameters)
  which.ok <- which(is.finite(logdens))
  theta.ok <- thetaparticles[which.ok,,drop=FALSE]
  logdens[which.ok] <- logdens[which.ok] + full_loglikelihood(theta.ok, obs)
  return(logdens)
}

# Function to perform Metropolis-Hastings
metropolishastings <- function(posterior, tuning_parameters){
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
    if (iteration %% 1000 == 1){
      cat("iteration ", iteration, "/", niterations, "\n")
      cat("average acceptance:", naccepts / (niterations*nchains) * 100, "%\n")
    }
    if (iteration > 50 && tuning_parameters$adaptation > 0  && (iteration %% tuning_parameters$adaptation) == 0){
      # adapt the proposal covariance matrix based on the last < 50,000 samples of all chains
      mcmc_samples <- foreach(ichain = 1:nchains, .combine = rbind) %do% {
        matrix(chains[[ichain]][max(1, iteration - 50000):(iteration-1),], ncol = p)
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
  }
  cat("average acceptance:", naccepts / (niterations*nchains) * 100, "%\n")
  return(list(chains = chains, naccepts = naccepts, cov_proposal = cov_proposal))
}

# tuning parameters
nchains <- 4
# initialize chain from prior
init_chains <- target$rprior(nchains, parameters = target$parameters)
tuning_parameters <- list(niterations = 10000, nchains = nchains, cov_proposal = diag(1e-1, nrow = target$thetadim, ncol = target$thetadim),
                          adaptation = 1000, init_chains = init_chains)

filename <- paste0("~/Dropbox/ABCD/Results/univariatenormal/normaldata.n", nobservations, ".metropolis.RData")
mh <- metropolishastings(posterior, tuning_parameters)
save(mh, nobservations, file = filename)
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
# # plot
# chaindf.melt <- melt(chaindf, id.vars = c("ichain", "iteration"))
#
# ggplot(chaindf.melt  %>% filter(iteration > 1000), aes(x = iteration, y = value, group = interaction(variable,ichain), colour = variable)) + geom_line()
#
# library(tidyr)
# chain.bycomponent.df <- chaindf.melt %>% filter(iteration > 2000) %>% spread(variable, value)
# # ggplot(chain.bycomponent.df, aes(x = X.1)) + geom_density(aes(y = ..density.. ))
#
# g <- ggplot(chain.bycomponent.df, aes(x = X.1)) + geom_histogram(aes(y = ..density..), fill = "grey", binwidth = 0.01)
# g <- g + geom_vline(xintercept = true_theta[1])
# g
#
# g <- ggplot(chain.bycomponent.df, aes(x = X.2)) + geom_histogram(aes(y = ..density..), fill = "grey", binwidth = 0.01)
# g <- g + geom_vline(xintercept = true_theta[2])
# g
#
