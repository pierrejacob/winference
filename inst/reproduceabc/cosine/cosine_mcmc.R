library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()
set.seed(11)
prefix <- ""
target <- get_cosine()

nobservations <- 100

load(paste0(prefix, "cosinedata.RData"))
obs <- matrix(obs[1:nobservations], nrow = 1)
loglikelihood <- function(thetas, ys, ...){
  evals <- rep(0, nrow(thetas))
  for (itheta in 1:nrow(thetas)){
    backbone <- exp(thetas[itheta,4]) * cos(2 * pi * thetas[itheta,1] * (1:nobservations) + thetas[itheta,2])
    evals[itheta] <- sum(dnorm(ys, mean = backbone, sd = exp(thetas[itheta,3]), log = TRUE))
  }
  return(evals)
}


target$loglikelihood <- loglikelihood
theta_init <- target$rprior(8, target$parameters)

tuning_parameters <- list(niterations = 100000, nchains = nrow(theta_init),
                          cov_proposal = diag(0.1, nrow = target$thetadim, ncol = target$thetadim),
                          adaptation = 10000, init_chains = theta_init)
mhfile <- paste0(prefix, "cosine.mcmc.n", nobservations, ".RData")
mh <- metropolishastings(obs, target, tuning_parameters)
save(mh, file = mhfile)
load(mhfile)

burnin <- 50000
chain.df <- mhchainlist_to_dataframe(mh$chains)
chain.df %>% head
g <- ggplot(chain.df %>% filter(iteration > burnin, iteration %% 100 == 1), aes(x = iteration, y = X.1, group = ichain, colour = factor(ichain))) + geom_line()
g +  geom_hline(yintercept = true_theta[1], col = "red")

g <- ggplot(chain.df %>% filter(iteration > burnin, iteration %% 100 == 1), aes(x = iteration, y = X.2, group = ichain, colour = factor(ichain))) + geom_line()
g +  geom_hline(yintercept = true_theta[2], col = "red")

g <- ggplot(chain.df %>% filter(iteration > burnin, iteration %% 100 == 1), aes(x = iteration, y = X.3, group = ichain, colour = factor(ichain))) + geom_line()
g +  geom_hline(yintercept = true_theta[3], col = "red")

g <- ggplot(chain.df %>% filter(iteration > burnin, iteration %% 100 == 1), aes(x = iteration, y = X.4, group = ichain, colour = factor(ichain))) + geom_line()
g +  geom_hline(yintercept = true_theta[4], col = "red")

