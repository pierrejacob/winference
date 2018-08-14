library(winference)
registerDoParallel(cores = 10)
rm(list = ls())
setmytheme()
setwd("~/Dropbox/ABCDresults/levydrivenstochvol//")
set.seed(13)
target <- get_levydriven()
prefix <- ""
load(file = paste0(prefix, "levydrivendata.RData"))
nobservations <- 10000
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

# nparticles <- 2^10
# res <- particle_filter(nparticles, model, post_mean, obs)


library(foreach)
library(doParallel)
registerDoParallel(cores = 10)
nparticles <- 4096
proposal <- rproposal(10)
dproposals <- dproposal(proposal)
prior_eval <- target$dprior(proposal)
nrep <- dim(proposal)[1]
nrep_per_theta <- 50
results.df <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  lls <- rep(0, nrep_per_theta)
  times <- rep(0, nrep_per_theta)
  if (is.finite(prior_eval[irep])){
    for (irep_per_theta in 1:nrep_per_theta){
      pct <- proc.time()
      lls[irep_per_theta] <- particle_filter(nparticles, model, proposal[irep,], obs)
      times[irep_per_theta] <- as.numeric((proc.time() - pct)[3])
    }
  }
  data.frame(irep = irep, irep_per_theta = 1:nrep_per_theta, lls = lls, times = times)
}

results.df %>% group_by(irep) %>% summarise(sdlls = sd(lls), meantime = mean(times))
