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
nobservations <- 10
load(file = "~/Dropbox/ABCD/Results/data/normaldata.RData")
obs <- obs[1:nobservations]

## load results from ABC-MCMC
load("~/Dropbox/ABCD/Results/univariatenormal/normaldata.n10.wmcmc.threshold0.5.RData")
load("~/Dropbox/ABCD/Results/univariatenormal/normaldata.n10.wmcmc.threshold0.25.RData")
load("~/Dropbox/ABCD/Results/univariatenormal/normaldata.n10.wmcmc.threshold0.1.RData")
## load results from standard MCMC using the exact likelihood
filename <- paste0("~/Dropbox/ABCD/Results/univariatenormal/normaldata.n", nobservations, ".metropolis.RData")
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
#############

res1$chains_df %>% head

g <- ggplot(res1$chains_df, aes(x = X.1)) + geom_density(aes(y = ..density..), colour = "darkgrey")
g <- g + xlab(expression(mu))
g <- g + geom_density(data = res2$chains_df, colour = "darkgrey")
g <- g + geom_density(data = res3$chains_df, colour = "darkgrey")
g <- g + geom_density(data = mcmc.df, colour = "black", linetype = 2)
g

g <- ggplot(res1$chains_df, aes(x = X.2)) + geom_density(aes(y = ..density..), colour = "darkgrey")
g <- g + xlab(expression(sigma))
g <- g + geom_density(data = res2$chains_df, colour = "darkgrey")
g <- g + geom_density(data = res3$chains_df, colour = "darkgrey")
g <- g + geom_density(data = mcmc.df, colour = "black", linetype = 2)
g

