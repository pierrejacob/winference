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
nobservations <- 10
load(file = "~/Dropbox/ABCD/Results/data/normaldata.RData")
obs <- obs[1:nobservations]

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

ggplot(chaindf.melt  %>% filter(iteration > 1000), aes(x = iteration, y = value, group = interaction(variable,ichain), colour = variable)) + geom_line()

library(tidyr)
chain.bycomponent.df <- chaindf.melt %>% filter(iteration > 2000) %>% spread(variable, value)
# ggplot(chain.bycomponent.df, aes(x = X.1)) + geom_density(aes(y = ..density.. ))

g <- ggplot(chain.bycomponent.df, aes(x = X.1)) + geom_histogram(aes(y = ..density..), fill = "grey", binwidth = 0.01)
g <- g + geom_vline(xintercept = true_theta[1])  + xlab(expression(mu))
g

g <- ggplot(chain.bycomponent.df, aes(x = X.2)) + geom_histogram(aes(y = ..density..), fill = "grey", binwidth = 0.01)
g <- g + geom_vline(xintercept = true_theta[2]) + xlab(expression(sigma))
g

