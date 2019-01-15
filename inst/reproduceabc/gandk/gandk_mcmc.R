#
#+ presets, echo = FALSE, warning = FALSE, message = FALSE
library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()
set.seed(11)
target <- get_gandk()

prefix <- ""

nobs <- 250
load(paste0(prefix, "gandkdata.RData"))
obs <- obs[1:nobs]


loglikelihood <- function(thetas, ys, ...){
  evals <- rep(0, nrow(thetas))
  for (itheta in 1:nrow(thetas)){
    faster_ll <- function(ys, h = 1e-5, tolerance = 1e-10){
      all_ys <- c(ys-h, ys+h)
      o <- order(all_ys)
      x <- rep(0, length(all_ys))
      x[o[1]] <- gandkcdf(y = all_ys[o[1]], theta = thetas[itheta,], tolerance = tolerance)
      for (i in 2:length(all_ys)){
        x[o[i]] <- gandkcdf(y = all_ys[o[i]], theta = thetas[itheta,], tolerance = tolerance, lower = x[o[i-1]])
      }
      return(sum(log((x[(nobs+1):(2*nobs)] - x[1:nobs])/(2*h))))
    }
    evals[itheta] <- faster_ll(ys)
  }
  return(evals)
}

target$loglikelihood <- loglikelihood

theta_init <- target$rprior(5, target$parameters)

cov <- structure(c(0.000962572396210165, 0.00179098940779821, -0.00191646746855593,
                   -0.00107821160254063, 0.00179098940779821, 0.00485894811881567,
                   -8.6591393768244e-05, -0.00266660529656441, -0.00191646746855593,
                   -8.6591393768244e-05, 0.0175762001060228, 0.00128474129479652,
                   -0.00107821160254063, -0.00266660529656441, 0.00128474129479652,
                   0.00216663223962476), .Dim = c(4L, 4L))
# tuning_parameters <- list(niterations = 7500, nchains = nrow(theta_init),
#                           cov_proposal = diag(0.1, nrow = target$thetadim, ncol = target$thetadim),
#                           adaptation = 2500, init_chains = theta_init)
tuning_parameters <- list(niterations = 75000, nchains = nrow(theta_init),
                          cov_proposal = cov,
                          adaptation = 0, init_chains = theta_init)

mhfile <- paste0(prefix, "gandkmcmc.n", nobs, "mh.RData")
mh <- metropolishastings(obs, target, tuning_parameters)
save(mh, file = mhfile)
load(mhfile)

burnin <- 50000
chain.df <- mhchainlist_to_dataframe(mh$chains)
chain.df %>% head
burnin <- 0
g <- ggplot(chain.df, aes(x = iteration, y = X.1, group = ichain)) + geom_line()
g
g <- ggplot(chain.df, aes(x = iteration, y = X.2, group = ichain)) + geom_line()
g
g <- ggplot(chain.df, aes(x = iteration, y = X.3, group = ichain)) + geom_line()
g
g <- ggplot(chain.df, aes(x = iteration, y = X.4, group = ichain)) + geom_line()
g

# burnin <- 5500
# g <- ggplot(chain.df %>% filter(iteration > burnin), aes(x = X.1, y = X.2)) + geom_point()
# g <- g + geom_vline(xintercept = theta[1])
# g <- g + geom_hline(yintercept = theta[2])
# g
#
#
# burnin <- 5500
# g <- ggplot(chain.df  %>% filter(iteration > burnin), aes(x = X.3, y = X.4)) + geom_point()
# g <- g + geom_vline(xintercept = theta[3])
# g <- g + geom_hline(yintercept = theta[4])
# g

# g <- ggplot(chain.df %>% filter(iteration > burnin, iteration %% 100 == 1), aes(x = iteration, y = X.1, group = ichain)) + geom_line()
# g

# g <- ggplot(chain.df %>% filter(iteration > burnin), aes(x = X.1)) + geom_histogram(aes(y = ..density..), fill = "grey", binwidth = 0.005)
# g <- g + xlab(expression(rho)) + geom_vline(xintercept = true_theta[1])
# g

# g <- ggplot(chain.df %>% filter(iteration > burnin), aes(x = X.2)) + geom_histogram(aes(y = ..density..), fill = "grey", binwidth = 0.005)
# g <- g + xlab(expression(log(sigma))) + geom_vline(xintercept = true_theta[2])
# g
#
# g <- ggplot(chain.df %>% filter(iteration > burnin, iteration %% 100 == 1), aes(x = X.1, y = X.2)) + geom_point()
# g <- g + xlab(expression(log(sigma))) + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
# g


