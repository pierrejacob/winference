#+ presets, echo = FALSE, warning = FALSE, message = FALSE
library(winference)
library(ggthemes)
library(doRNG)
library(tidyr)
registerDoMC(cores = 10)
rm(list = ls())
setmytheme()

set.seed(11)

target <- get_normal()

# number of observations
nobservations <- 100
load(file = "~/Dropbox/ABCD/Results/data/gammadata.RData")
obs <- obs[1:nobservations]

## load results from ABC-SMC
load(paste0("~/Dropbox/ABCD/Results/univariatenormal/gammadata.n", nobservations, ".wsmc_rhit.RData"))
## load results from standard MCMC using the exact likelihood
filename <- paste0("~/Dropbox/ABCD/Results/univariatenormal/gammadata.n", nobservations, ".metropolis.RData")
load(filename)
chainlist_to_dataframe <- function(chains_list){
  nchains <- length(chains_list)
  niterations <- nrow(chains_list[[1]])
  chaindf <- foreach (i = 1:nchains, .combine = rbind) %do% {
    data.frame(ichain = rep(i, niterations), iteration = 1:niterations, X = chains_list[[i]])
  }
  return(chaindf)
}
chaindf <- chainlist_to_dataframe(mh$chains)
# plot
chaindf.melt <- melt(chaindf, id.vars = c("ichain", "iteration"))
mcmc.df <- chaindf.melt %>% filter(iteration > 2000) %>% spread(variable, value)
wsmcresults$threshold_history
wsmc.df <- wsmc_to_dataframe(wsmcresults, target$parameter_names)
#############
g <- ggplot(wsmc.df %>% filter(step > 2, step < 25), aes(x = mu, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g +  theme(legend.position = "none")
g <- g + geom_density(data = mcmc.df, aes(x = X.1, group = NULL, colour = NULL), colour = "black", linetype = 2)
g <- g + xlab(expression(mu))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
g
ggsave(filename = "~/Dropbox/ABCD/draft5/quasibayes_mu_gamma10_5_n100.pdf", plot = g, width = 7, height = 5 )

g <- ggplot(wsmc.df %>% filter(step > 2, step < 25), aes(x = sigma, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none")
g <- g + geom_density(data = mcmc.df, aes(x = X.2, group = NULL, colour = NULL), colour = "black", linetype = 2)
g <- g + xlab(expression(sigma))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
g
ggsave(filename = "~/Dropbox/ABCD/draft5/quasibayes_sigma_gamma10_5_n100.pdf", plot = g, width = 7, height = 5 )

# sqrt(10/25)
# wsmcresults$threshold_history
