library(winference)
registerDoParallel(cores = 10)
rm(list = ls())
setmytheme()
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


# particle_filter(nparticles, model, mean_proposal, obs)


library(foreach)
library(doParallel)
registerDoParallel(cores = 10)
nparticles <- 8192

proposal <- rproposal(1000)
dproposals <- dproposal(proposal)
prior_eval <- target$dprior(proposal)
nrep <- dim(proposal)[1]
res <- foreach(irep = 1:nrep, .combine = c) %dorng% {
  ll <- 0
  if (is.finite(prior_eval[irep])){
    ll <- particle_filter(nparticles, model, proposal[irep,], obs)
  }
  ll
}

is_weight <- prior_eval + res - dproposals
maxlw <- max(is_weight)
w <- exp(is_weight - maxlw)
normweights <- w / sum(w)
1/(sum(normweights^2))

# qplot(x = proposal[,1], geom = "histogram")
df <- data.frame(proposal)
df$normweights <- normweights
head(df)

ggplot(df, aes(x = X1)) + geom_density(aes(y = ..density..)) + geom_density(aes(weight = normweights, y = ..density..), colour = "blue")
ggplot(df, aes(x = X2)) + geom_density(aes(y = ..density..)) + geom_density(aes(weight = normweights, y = ..density..), colour = "blue")
ggplot(df, aes(x = X3)) + geom_density(aes(y = ..density..)) + geom_density(aes(weight = normweights, y = ..density..), colour = "blue")
ggplot(df, aes(x = X4)) + geom_density(aes(y = ..density..)) + geom_density(aes(weight = normweights, y = ..density..), colour = "blue")
ggplot(df, aes(x = X5)) + geom_density(aes(y = ..density..)) + geom_density(aes(weight = normweights, y = ..density..), colour = "blue") + scale_x_log10()

## One could fit a Normal on the weighted samples, and try again...
post_mean2 <- wmean(proposal, normweights)
post_cov2 <- wcovariance(proposal, normweights, post_mean2)

post_mean
post_mean2
post_cov
post_cov2

proposal2 <- fast_rmvnorm(1000, post_mean2, post_cov2)
dproposals2 <- dproposal(proposal2)
prior_eval2 <- target$dprior(proposal2)
nrep <- dim(proposal2)[1]
res2 <- foreach(irep = 1:nrep, .combine = c) %dorng% {
  ll <- 0
  if (is.finite(prior_eval2[irep])){
    ll <- particle_filter(nparticles, model, proposal2[irep,], obs)
  }
  ll
}
is_weight2 <- prior_eval2 + res2 - dproposals2
maxlw2 <- max(is_weight2)
w2 <- exp(is_weight2 - maxlw2)
normweights2 <- w2 / sum(w2)
1/(sum(normweights2^2))
# put in data.frame
df2 <- data.frame(proposal2)
df2$normweights <- normweights

ggplot(df, aes(x = X1)) + geom_density(aes(y = ..density..)) + geom_density(aes(weight = normweights, y = ..density..), colour = "blue") +
  geom_density(data = df2, aes(y = ..density..), colour = "purple") + geom_density(data = df2, aes(weight = normweights, y = ..density..), colour = "orange")

ggplot(df, aes(x = X2)) + geom_density(aes(y = ..density..)) + geom_density(aes(weight = normweights, y = ..density..), colour = "blue") +
  geom_density(data = df2, aes(y = ..density..), colour = "purple") + geom_density(data = df2, aes(weight = normweights, y = ..density..), colour = "orange")

ggplot(df, aes(x = X3)) + geom_density(aes(y = ..density..)) + geom_density(aes(weight = normweights, y = ..density..), colour = "blue") +
  geom_density(data = df2, aes(y = ..density..), colour = "purple") + geom_density(data = df2, aes(weight = normweights, y = ..density..), colour = "orange")

ggplot(df, aes(x = X4)) + geom_density(aes(y = ..density..)) + geom_density(aes(weight = normweights, y = ..density..), colour = "blue") +
  geom_density(data = df2, aes(y = ..density..), colour = "purple") + geom_density(data = df2, aes(weight = normweights, y = ..density..), colour = "orange")

ggplot(df, aes(x = X5)) + geom_density(aes(y = ..density..)) + geom_density(aes(weight = normweights, y = ..density..), colour = "blue") +
  geom_density(data = df2, aes(y = ..density..), colour = "purple") + geom_density(data = df2, aes(weight = normweights, y = ..density..), colour = "orange")

# plot_bivariate(results, 1, 2)
plot(thetas[,1], thetas[,2])
points(proposal2[,1], proposal2[,2], col = "red")

plot(thetas[,3], thetas[,4])
points(proposal2[,3], proposal2[,4], col = "red")

plot(thetas[,4], thetas[,5])
points(proposal2[,4], proposal2[,5], col = "red")


# ggplot(new_df, aes(x = X1)) + geom_density(aes(y = ..density..)) + geom_density(aes(weight = normweights, y = ..density..), colour = "blue")
# ggplot(new_df, aes(x = X2)) + geom_density(aes(y = ..density..)) + geom_density(aes(weight = normweights, y = ..density..), colour = "blue")
# ggplot(new_df, aes(x = X3)) + geom_density(aes(y = ..density..)) + geom_density(aes(weight = normweights, y = ..density..), colour = "blue")
# ggplot(new_df, aes(x = X4)) + geom_density(aes(y = ..density..)) + geom_density(aes(weight = normweights, y = ..density..), colour = "blue")
# ggplot(new_df, aes(x = X5)) + geom_density(aes(y = ..density..)) + geom_density(aes(weight = normweights, y = ..density..), colour = "blue")
#
# summary(new_is_weight[is.finite(new_is_weight)])
#
# # ggplot(df, aes(x = X1, y = X2)) + geom_density_2d() + geom_density_2d(aes(weight = normweights), colour = "red")
# # ggplot(df, aes(x = X3, y = X4)) + geom_density_2d() + geom_density_2d(aes(weight = normweights), colour = "red")
# # ggplot(df, aes(x = X4, y = X5)) + geom_density_2d() + geom_density_2d(aes(weight = normweights), colour = "red")
#
# # ancestors <- systematic_resampling(normweights)
# # thetas_resampled <- df[ancestors,1:4]
# # unique(ancestors)
# # fit <- mixmodCluster(data = data.frame(thetas_resampled), nbCluster = 5, dataType = "quantitative")


#timings
nparticles <- 1000
library(microbenchmark)
microbenchmark(particle_filter(nparticles, model, true_theta, obs), times = 1e2)
