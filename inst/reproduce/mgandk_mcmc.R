library(winference)
registerDoParallel(cores = 6)
rm(list = ls())
setmytheme()
set.seed(11)
target <- get_mgandk()
prefix <- ""

nobservations <- 500
load(paste0(prefix, "mgandkdata.RData"))
obs <- obs[,1:nobservations]
# load ABC posterior
load(paste0(prefix, "mgandk.wsmc.n1000.hilbert.RData"))
thetas <- tail(results$thetas_history, 1)[[1]]
theta_init <- thetas[sample(x = 1:nrow(thetas), 8, replace = TRUE),]
colMeans(thetas)
cov <- cov(thetas)
# test log-likelihood
target$loglikelihood(theta_init, obs)
tuning_parameters <- list(niterations = 100000, nchains = nrow(theta_init),
                          cov_proposal = cov,
                          adaptation = 1000, init_chains = theta_init)
mhfile <- paste0(prefix, "mgandk.mcmc.n", nobservations, ".mh.initfromABC.RData")
mh <- metropolishastings(obs, target, tuning_parameters, savefile = mhfile)
save(mh, file = mhfile)
load(mhfile)
if (exists("mh_results")){
  mcmc.df <- mhchainlist_to_dataframe(mh_results$chains)
  mcmc.df <- mcmc.df %>% filter(iteration < mh_results$iteration)
} else {
  mcmc.df <- mhchainlist_to_dataframe(mh$chains)
}

# mcmc.df %>% head
# burnin <- 5000
# # burnin <- 0
# ggplot(mcmc.df %>% filter(iteration > burnin, iteration %% 100 == 1), aes(x = iteration, y = X.7, group = ichain, colour = factor(ichain))) + geom_line(alpha = 0.5)
# ggplot(mcmc.df %>% filter(iteration > burnin, iteration %% 100 == 1), aes(x = iteration, y = X.9, group = ichain, colour = factor(ichain))) + geom_line(alpha = 0.5)
# ggplot(mcmc.df %>% filter(iteration > burnin), aes(x = X.1, group = ichain, colour = factor(ichain))) + geom_density() + geom_vline(xintercept = true_theta[1])
# ggplot(mcmc.df %>% filter(iteration > burnin), aes(x = X.2, group = ichain, colour = factor(ichain))) + geom_density() + geom_vline(xintercept = true_theta[2])
# ggplot(mcmc.df %>% filter(iteration > burnin), aes(x = X.3, group = ichain, colour = factor(ichain))) + geom_density() + geom_vline(xintercept = true_theta[3])
# ggplot(mcmc.df %>% filter(iteration > burnin), aes(x = X.4, group = ichain, colour = factor(ichain))) + geom_density() + geom_vline(xintercept = true_theta[4])
# ggplot(mcmc.df %>% filter(iteration > burnin), aes(x = X.5, group = ichain, colour = factor(ichain))) + geom_density() + geom_vline(xintercept = true_theta[5])
# ggplot(mcmc.df %>% filter(iteration > burnin), aes(x = X.6, group = ichain, colour = factor(ichain))) + geom_density() + geom_vline(xintercept = true_theta[6])
# ggplot(mcmc.df %>% filter(iteration > burnin), aes(x = X.7, group = ichain, colour = factor(ichain))) + geom_density() + geom_vline(xintercept = true_theta[7])
# ggplot(mcmc.df %>% filter(iteration > burnin), aes(x = X.8, group = ichain, colour = factor(ichain))) + geom_density() + geom_vline(xintercept = true_theta[8])
# ggplot(mcmc.df %>% filter(iteration > burnin), aes(x = X.9, group = ichain, colour = factor(ichain))) + geom_density() + geom_vline(xintercept = true_theta[9])
#

