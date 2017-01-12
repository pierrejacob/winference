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

# function to compute the quasi-loglikelihood, given a threshold
qloglikelihood <- function(theta, threshold = Inf, M = 1){
  ds <- rep(0, M)
  for (i in 1:M){
    ds[i] <- compute_d(theta)
  }
  return(list(qvalue = log(mean(ds < threshold)), distances = ds))
}

proposal <- randomwalk_proposal()
proposal$param_prop$cov <- diag(c(1, 1))
nmoves <- 10000
nchains <- 4
# adapts the proposal every 500 iterations
# does not attempt to run the chains in parallel
param_algo <- list(nmoves = nmoves, nchains = nchains, proposal = proposal,
                   M = 20, adaptation = 1000, parallel = FALSE)
#
# ini1 <- initialization(threshold, qloglikelihood, target, param_algo)
# res1 <- quasimh_steps(ini1$thetas, ini1$priorvalues, ini1$qvalues, ini1$distances, threshold, qloglikelihood,  target, param_algo)
threshold <- 0.5
filename <- paste0("~/Dropbox/ABCD/Results/univariatenormal/gammadata.n",
                   nobservations, ".wmcmc.threshold", threshold, ".RData")
res1 <- quasimh(threshold, qloglikelihood, target, param_algo)
save(res1, file = filename)

# lower threshold
threshold <- 0.25
param_algo <- res1$param_algo
#
filename <- paste0("~/Dropbox/ABCD/Results/univariatenormal/gammadata.n",
                   nobservations, ".wmcmc.threshold", threshold, ".RData")
res2 <- quasimh_continue(res1, threshold, qloglikelihood, target, param_algo)
save(res2, file = filename)

threshold <- 0.1
param_algo <- res2$param_algo
#
filename <- paste0("~/Dropbox/ABCD/Results/univariatenormal/gammadata.n",
                   nobservations, ".wmcmc.threshold", threshold, ".RData")
res3 <- quasimh_continue(res2, threshold, qloglikelihood, target, param_algo)
save(res3, file = filename)


# ggplot(res1$chains_df, aes(x = iteration, y = X.1, group = ichain)) + geom_line()
# ggplot(res1$chains_df, aes(x = iteration, y = X.2, group = ichain)) + geom_line()


# chains.df <- res2$chains_df %>% filter(iteration > param_algo$nmoves / 2)
# muwmcmc <- chains.df$X.1
# sigmawcmc <- chains.df$X.2
#
# threshold <- 0.035
# param_algo <- res2$param_algo
# param_algo$nmoves <- 1000
# param_algo$nchains <- 2
# param_algo$M <- 250
# #
# res3 <- quasimh_continue(res2, threshold, qloglikelihood, target, param_algo)
# # or we could start from scratch
# res4 <- quasimh(threshold, qloglikelihood, target, param_algo)
#
# res3$chains_df %>% tail
#
# chains.df <- res4$chains_df # %>% filter(iteration > res3$param_algo$nmoves / 2)
#
# muwmcmc <- chains.df$X.1
# sigmawcmc <- chains.df$X.2
#
# # head(chains.df)
# ggplot(chains.df, aes(x = iteration, y = X.1, group = ichain)) + geom_line()
# ggplot(chains.df, aes(x = iteration, y = X.2, group = ichain)) + geom_line()
# #
# # ggplot(chains.df, aes(x = X.1)) + geom_histogram(aes(y = ..density..), binwidth = 0.01) +
# #   stat_function(fun = dnorm)
# #
# # ggplot(chains.df, aes(x = X.2)) + geom_histogram(aes(y = ..density..), binwidth = 0.01) +
# #   stat_function(fun = function(x) dgamma(x, 2, 1))
#
#
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


# # rare event simulation approach
#
# compute_d_given_rand <- function(theta, fake_rand, metric = metricL2){
#   fake_obs <- target$robservation(nobservations, theta, target$parameters, fake_rand)
#   fake_obs_sorted <- sort(fake_obs)
#   return(metric(obs_sorted, fake_obs_sorted))
# }
#
# smcrareevent <- function(theta, N, Nacc, threshold_min, terminate = -Inf){
#   threshold_previous <- Inf
#   threshold_history <- c(threshold_previous)
#   x <- matrix(nrow = N, ncol = nobservations)
#   phis <- rep(0, N)
#   # step 1
#   for (i in 1:N){
#     x[i,] <- target$generate_randomness(nobservations)
#     phis[i] <- compute_d_given_rand(theta, x[i,])
#   }
#   logLhat <- 0
#   iter <- 0
#   while(threshold_previous > threshold_min && logLhat > terminate){
#     threshold_new <- max(phis[order(phis)[Nacc]], threshold_min)
#     It <- which(phis <= threshold_new)
#     logLhat <- logLhat + log(length(It)/N)
#     # ancestors <- sample(x = It, size = N, replace = TRUE)
#     weights <- (phis <= threshold_new) / length(It)
#     ancestors <- systematic_resampling(weights)
#     x <- x[ancestors,]
#     phis <- phis[ancestors]
#     for (i in 1:N){
#       res <- elliptical_update(x[i,], theta, 0, compute_d_given_rand, threshold_previous, function(distances, threshold){ return(log(distances <= threshold))})
#       x[i,] <- res$fprime
#       phis[i] <- res$dist
#     }
#     # cat("threshold = ", threshold_new, ", summary of phis = ", summary(phis), "\n")
#     threshold_previous <- threshold_new
#     threshold_history <- c(threshold_history, threshold_previous)
#     iter <- iter + 1
#     if (iter > 1e3){
#       print("already 1000 iterations in the SMC calculations... aborting")
#       return(NULL)
#     }
#   }
#   return(list(threshold_history = threshold_history, logLhat = logLhat))
# }
#
# mh_rareevent <- function(nmoves, theta, threshold, target, param_algo){
#   p <- length(theta)
#   # store whole chains
#   chain <- matrix(nrow = nmoves, ncol = p)
#   # current states of the chain
#   current_chain <- theta
#   # initialization of the chain
#   chain[1,] <- theta
#   # log target density values associated with the current states of the chains
#   current_dtarget <-  smcrareevent(current_chain, param_algo$N, param_algo$Nacc, threshold)$logLhat
#   logdens <- target$dprior(matrix(current_chain, nrow = 1), target$parameters)
#   current_dtarget <- current_dtarget + logdens
#   #
#   naccepts <- 0
#   # run the chains
#   for (iteration in 2:nmoves){
#     if (iteration %% 100 == 1){
#       cat("iteration ", iteration, "/", nmoves, "\n")
#       cat("average acceptance:", naccepts / nmoves * 100, "%\n")
#     }
#     if (iteration > 50 && param_algo$adaptation > 0  && (iteration %% param_algo$adaptation) == 0){
#       # adapt the proposal covariance matrix based on the last < 50,000 samples of all chains
#       mcmc_samples <- matrix(chain[max(1, iteration - 50000):(iteration-1),], ncol = p)
#       param_algo$proposal$param_prop$cov <- cov(mcmc_samples) / p
#     }
#     # proposals
#     proposal <- param_algo$proposal$r(matrix(current_chain, nrow = 1), param_algo$proposal$param_prop)
#     dproposal <- param_algo$proposal$d(proposal, param_algo$proposal$param_prop)
#     dcurrent <- param_algo$proposal$d(matrix(current_chain, nrow = 1), param_algo$proposal$param_prop)
#     # proposals' target density
#     proposal_dtarget <- target$dprior(matrix(proposal, nrow = 1), target$parameters)
#     uniform <- runif(n = 1)
#     terminate <- (log(uniform) + (current_dtarget - proposal_dtarget) + (dcurrent - dproposal))
#     if (is.finite(proposal_dtarget)){
#       proposal_dtarget <- proposal_dtarget + smcrareevent(proposal, param_algo$N, param_algo$Nacc, threshold, terminate)$logLhat
#     }
#     # log Metropolis Hastings ratio
#     acceptance_ratios <- (proposal_dtarget - current_dtarget) + (dproposal - dcurrent)
#     # acceptance decisions
#     accept <- (log(uniform) < acceptance_ratios)
#     naccepts <- naccepts + accept
#     # make the appropriate replacements
#     if (accept){
#       current_chain <- proposal
#       current_dtarget <- proposal_dtarget
#     }
#     # book keeping
#     chain[iteration,] <- current_chain
#   }
#   cat("average acceptance:", naccepts / nmoves * 100, "%\n")
#   return(list(chain = chain, naccepts = naccepts, param_algo = param_algo))
# }
#
# threshold
# param_algo_RE <- list(adaptation = 0, N = 256, Nacc = 128, proposal = param_algo$proposal)
# ini <- find_init_values(res2, threshold)
# theta <- ini$thetas[1,]
#
# nmoves <- 5000
# mh_RE <- mh_rareevent(nmoves, theta, threshold, target, param_algo_RE)
#
# matplot(mh_RE$chain, type = "l")
#
# hist(muwmcmc, nclass = 20, prob = TRUE, col = rgb(1,0,0, alpha = 0.5))
# hist(mumcmc, nclass = 50, prob = TRUE, add = TRUE, col = rgb(0,0,0, alpha = 0.5))
# hist(mh_RE$chain[,1], nclass = 50, prob = TRUE, add = TRUE, col = rgb(0,0,1, alpha = 0.5))
