#+ presets, echo = FALSE, warning = FALSE, message = FALSE
library(winference)
library(doMC)
library(doRNG)
library(dplyr)
library(ggthemes)
library(tidyr)

registerDoMC(cores = 6)

rm(list = ls())
set.seed(11)
setmytheme()
#
# model
target <- get_autoregressive()
#
# number of observations
nobservations <- 1000
load(file = "~/Dropbox/ABCD/Results/data/ar1data.RData")

### load the results obtained by only looking at the empirical distribution of the data,
### without delay reconstruction
filename <- paste0("~/Dropbox/ABCD/Results/autoregressive/ar1data.n1000.wsmc.RData")
load(filename)
wsmc.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps <- max(wsmc.df$step)

# g <- ggplot(wsmc.df, aes(x = rho, group = step)) + geom_density(aes(y = ..density..), colour = "darkgrey")
# g <- g +  theme(legend.position = "none")
# g <- g + xlab(expression(rho))
# g
#
# g <- ggplot(wsmc.df, aes(x = logsigma, group = step)) + geom_density(aes(y = ..density..), colour = "darkgrey")
# g <- g +  theme(legend.position = "none")
# g <- g + xlab(expression(log(sigma)))
# g


### load the results obtained with delay reconstruction, with a lag of 1
### with hilbert sort
# filename <- paste0("~/Dropbox/ABCD/Results/autoregressive/ar1data.n1000.wsmc-hilbert.RData")
# load(filename)
# wsmc.hilbert.df <- wsmc_to_dataframe(results, target$parameter_names)
# nsteps.hilbert <- max(wsmc.hilbert.df$step)
#
# g <- ggplot(wsmc.hilbert.df, aes(x = rho, y = logsigma, colour = step, group = step))
# g <- g + geom_point(alpha = 0.5)
# g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
# g <- g + xlab(expression(rho)) + ylab(expression(log(sigma)))
# g

### load the results obtained with delay reconstruction, with a lag of 1
### with Wasserstein distance
filename <- paste0("~/Dropbox/ABCD/Results/autoregressive/ar1data.n1000.wsmc_rhit.Lag1.wassersteinL2.RData")
load(filename)
wsmc.wasserstein.df <- wsmc_to_dataframe(results, target$parameter_names)
nsteps.wasserstein <- max(wsmc.wasserstein.df$step)


## load MH samples approximating posterior distribution
load("~/Dropbox/ABCD/Results/autoregressive/ar1data.n1000.metropolis.RData")
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
# check mixing
ggplot(chaindf.melt  %>% filter(iteration > 5000), aes(x = iteration, y = value, group = interaction(variable,ichain), colour = variable)) + geom_line()
#
chain.bycomponent.df <- chaindf.melt %>% filter(iteration > 5000) %>% spread(variable, value)
# g <- ggplot(chain.bycomponent.df, aes(x = X.1)) + geom_histogram(aes(y = ..density..), fill = "grey", binwidth = 0.01)
# g <- g + geom_vline(xintercept = true_theta[1])  + xlab(expression(rho))
# g
# g <- ggplot(chain.bycomponent.df, aes(x = X.2)) + geom_histogram(aes(y = ..density..), fill = "grey", binwidth = 0.01)
# g <- g + geom_vline(xintercept = true_theta[2]) + xlab(expression(log(sigma)))
# g
# g <- ggplot(chain.bycomponent.df %>% filter(iteration %% 10 == 1), aes(x = X.1, y = X.2, colour = ichain, group = ichain))
# g <- g + geom_point(alpha = 0.5)
# g <- g + theme(legend.position = "none")
# g <- g + xlab(expression(rho)) + ylab(expression(log(sigma))) + xlim(-1,1) + ylim(-4,4)
# g
#
# g <- ggplot(chain.bycomponent.df, aes(x = X.1, y = X.2))
# g <- g + geom_density2d()
# g <- g + theme(legend.position = "none")
# g <- g + xlab(expression(rho)) + ylab(expression(log(sigma)))  + xlim(-1,1) + ylim(-3,3)
# g



### now plot figures
g <- ggplot(wsmc.df, aes(x = rho, y = logsigma, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(phi)) + ylab(expression(log(sigma))) +  xlim(-1,1) + ylim(-4,4)
g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
# g
# ggsave(filename = "~/Dropbox/ABCD/draft4/ar1_identification_withoutlag.pdf", plot = g, width = 7, height = 5 )
ggsave(filename = "~/Dropbox/ABCD/draft4/ar1_identification_withoutlag.png", plot = g, width = 7, height = 5, dpi = 150)

g <- ggplot(wsmc.wasserstein.df, aes(x = rho, y = logsigma, colour = step, group = step))
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = floor(nsteps/2)) + theme(legend.position = "none")
g <- g + xlab(expression(phi)) + ylab(expression(log(sigma))) +  xlim(-1,1) + ylim(-4,4)
g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
# g
# ggsave(filename = "~/Dropbox/ABCD/draft4/ar1_identification_withlag1.pdf", plot = g, width = 7, height = 5 )
ggsave(filename = "~/Dropbox/ABCD/draft4/ar1_identification_withlag1.png", plot = g, width = 7, height = 5, dpi = 150)

g <- ggplot(wsmc.wasserstein.df %>% filter(step > 2), aes(x = rho, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(phi))
g <- g + geom_density(data=chain.bycomponent.df, aes(x = X.1, y = ..density.., group = NULL, colour = NULL),
                 linetype = 2, colour = "black")
g
ggsave(filename = "~/Dropbox/ABCD/draft4/ar1_withlag1_phi.pdf", plot = g, width = 7, height = 5)
# ggsave(filename = "~/Dropbox/ABCD/draft4/ar1_withlag1_phi.png", plot = g, width = 7, height = 5, dpi = 150)

g <- ggplot(wsmc.wasserstein.df %>% filter(step > 2), aes(x = logsigma, group = step, colour = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(log(sigma)))
g <- g + geom_density(data=chain.bycomponent.df, aes(x = X.2, y = ..density.., group = NULL, colour = NULL),
                      linetype = 2, colour = "black")
g
ggsave(filename = "~/Dropbox/ABCD/draft4/ar1_withlag1_logsigma.pdf", plot = g, width = 7, height = 5)
# ggsave(filename = "~/Dropbox/ABCD/draft4/ar1_withlag1_logsigma.png", plot = g, width = 7, height = 5, dpi = 150)

# g + geom_vline(xintercept = true_theta[2])




