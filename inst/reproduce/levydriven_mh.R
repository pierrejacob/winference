library(winference)
registerDoParallel(cores = 10)
rm(list = ls())
setmytheme()
setwd("~/Dropbox/")
set.seed(13)

target <- get_levydriven()
prefix <- ""
load(file = paste0(prefix, "levydrivendata.RData"))
nobservations <- 1000
obs <- obs[1:nobservations]
filename2 <- paste0(prefix, "levydriven.n", nobservations, ".lag1.wsmc.hilbert.summary.RData")
load(filename2)

thetas <- results$thetas_history[[results$thetas_history %>% length]]

post_cov <- cov(thetas)
post_mean <- colMeans(thetas)

variancescaling <- 1

rproposal <- function(n){
  return(fast_rmvnorm(n, post_mean, post_cov))
}

dproposal <- function(thetas){
  fast_dmvnorm(thetas, post_mean, post_cov)
}

# What if we did MCMC on this distribution?
nmcmc <- 5000
rinit <- function(){
  current <- rproposal(1)
  target_pdf <- dproposal(current)
  return(list(current = current, target_pdf = target_pdf))
}
kern <- function(current, target_pdf){
  proposal <- current + fast_rmvnorm(1, rep(0, length(current)), variancescaling * post_cov)
  proposal_pdf <- dproposal(proposal)
  if (log(runif(1)) < (proposal_pdf - target_pdf)){
    return(list(current = proposal, target_pdf = proposal_pdf, accept = 1))
  } else {
    return(list(current = current, target_pdf = target_pdf, accept = 0))
  }
}

nmcmc <- 1000
res_rinit <- rinit()
current <- res_rinit$current
target_pdf <- res_rinit$target_pdf
chain <- matrix(nrow = nmcmc, ncol = length(current))
naccepts <- 0
for (imcmc in 1:nmcmc){
  res <- kern(current, target_pdf)
  current <- res$current
  target_pdf <- res$target_pdf
  naccepts <- naccepts + res$accept
  chain[imcmc,] <- current
}
print(naccepts/nmcmc*100)
matplot(chain[,1:2], type = "l")
matplot(chain[,3:4], type = "l")
matplot(chain[,5], type = "l")
# it would work like a charm
# so now we envision PMCMC

model <- list()
model$rinit <- function(nparticles, theta){
  x <- matrix(nrow = nparticles, ncol = 2)
  for (i in 1:nparticles){
    x[i,] <- rgamma(2, shape = theta[3] * theta[3]/theta[4], scale = theta[4]/theta[3])
  }
  return(x)
}

model$rtransition <- function(xparticles, theta){
  rtransition_r <- levydriven_rtransition_rand(dim(xparticles)[1], theta)
  new_z <- exp(-theta[5]) * xparticles[,2] + rtransition_r$sum_weighted_e
  new_v <- (1/theta[5]) * (xparticles[,2] - new_z + rtransition_r$sum_e)
  xparticles[,1] <- new_v
  xparticles[,2] <- new_z
  return(xparticles)
}

model$dobs <- function(observations, time, xparticles, theta){
  return(dnorm(observations[time], mean = theta[1] + theta[2] * xparticles[,1], sd = sqrt(xparticles[,1]), log = TRUE))
}

### particle filter estimates of the log-likelihood
particle_filter <- function(nparticles, model, theta, observations){
  datalength <- nobservations
  # initialization
  xparticles <- model$rinit(nparticles, theta)
  logw <- rep(0, nparticles)
  # logw <- model$dobs(observations, 1, xparticles, theta)
  if (all(is.infinite(logw))){
    return(NA)
  }
  maxlw <- max(logw)
  w <- exp(logw - maxlw)
  # update log likelihood estimate
  ll <- maxlw + log(mean(w))
  normweights <- w / sum(w)
  #
  # step t > 1
  for (time in 1:datalength){
    if (time > 1){
      ancestors <- systematic_resampling_given_u(normweights, runif(1))
      xparticles <- xparticles[ancestors,,drop=F]
    }
    xparticles <- model$rtransition(xparticles, theta)
    logw <- model$dobs(observations, time, xparticles, theta)
    if (all(is.infinite(logw))){
      return(NA)
    }
    maxlw <- max(logw)
    w <- exp(logw - maxlw)
    # update log likelihood estimate
    ll <- ll + maxlw + log(mean(w))
    normweights <- w / sum(w)
    #
  }
  return(ll)
}

# How many particles do we need?
nparticles <- 8192
# nparticles <- 1024
# lls <- foreach(irep = 1:100, .combine = c) %dorng% {
#   particle_filter(nparticles, model, post_mean, obs)
# }
# summary(lls)
# sd(lls)

rinit_pmmh <- function(){
  valid <- FALSE
  while(!valid){
    current <- rproposal(1)
    prior_pdf <- target$dprior(current)
    valid <- is.finite(prior_pdf)
  }
  target_pdf <- prior_pdf + particle_filter(nparticles, model, current, obs)
  return(list(current = current, target_pdf = target_pdf))
}

kern_pmmh <- function(current, target_pdf){
  proposal <- current + fast_rmvnorm(1, rep(0, length(current)), variancescaling * post_cov)
  proposal_pdf <- target$dprior(proposal)
  if (is.finite(proposal_pdf)){
    proposal_pdf <- proposal_pdf + particle_filter(nparticles, model, proposal, obs)
    if (log(runif(1)) < (proposal_pdf - target_pdf)){
      return(list(current = proposal, target_pdf = proposal_pdf, accept = 1))
    } else {
      return(list(current = current, target_pdf = target_pdf, accept = 0))
    }
  } else {
    return(list(current = current, target_pdf = target_pdf, accept = 0))
  }
}
#
# nmcmc <- 10000
# res_rinit <- rinit_pmmh()
# current <- res_rinit$current
# target_pdf <- res_rinit$target_pdf
# chain <- matrix(nrow = nmcmc, ncol = length(current))
# chain_pdf <- rep(0, nmcmc)
# naccepts <- 0
# for (imcmc in 1:nmcmc){
#   res <- kern_pmmh(current, target_pdf)
#   current <- res$current
#   target_pdf <- res$target_pdf
#   naccepts <- naccepts + res$accept
#   chain[imcmc,] <- current
#   chain_pdf[imcmc] <- target_pdf
#   if (imcmc %% 10 == 1){
#     cat("iteration", imcmc, "\n")
#     cat("acceptance rate", naccepts / imcmc * 100, "%\n")
#     save(naccepts, imcmc, chain, chain_pdf, file = paste0("~/Dropbox/levydriven.n1000.N", nparticles, ".pmcmc.RData"))
#   }
# }
# save(naccepts, nmcmc, nparticles, chain, chain_pdf, file = paste0("~/Dropbox/levydriven.n1000.N", nparticles, ".pmcmc.RData"))
# print(naccepts/nmcmc*100)

load(paste0("~/Dropbox/levydriven.n1000.N", nparticles, ".pmcmc.RData"))
print(naccepts/nmcmc*100)
plot(chain_pdf[1:nmcmc], type = "l")

matplot(chain[,1:2], type = "l")
matplot(chain[,3:4], type = "l")
matplot(chain[,5], type = "l")

hist(thetas[,1], prob = TRUE, ylim = c(0,28), nclass = 50)
hist(chain[,1], prob = TRUE, add = TRUE, col = rgb(1,0,0,.5), nclass = 50)

hist(thetas[,2], prob = TRUE, ylim = c(0,8), nclass = 50)
hist(chain[,2], prob = TRUE, add = TRUE, col = rgb(1,0,0,.5), nclass = 50)

hist(thetas[,3], prob = TRUE, ylim = c(0,10), nclass = 50, xlim = c(0.1, 0.8))
hist(chain[,3], prob = TRUE, add = TRUE, col = rgb(1,0,0,.5), nclass = 50)

hist(thetas[,4], prob = TRUE, ylim = c(0,100), nclass = 50)
hist(chain[,4], prob = TRUE, add = TRUE, col = rgb(1,0,0,.5), nclass = 50)

hist(thetas[,5], prob = TRUE, ylim = c(0,500), nclass = 50)
hist(chain[,5], prob = TRUE, add = TRUE, col = rgb(1,0,0,.5), nclass = 50)

