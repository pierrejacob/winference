library(winference)
registerDoParallel(cores = detectCores()-2)
rm(list = ls())
setmytheme()
set.seed(11)
target <- get_mgandk()

prefix = ""

nobservations <- 500
load(paste0(prefix, "mgandkdata.RData"))
obs <- obs[,1:nobservations]
target$simulate <- function(theta) target$robservation(nobservations, theta)

library(transport)
compute_d <- function(z){
  sink(file = paste0(prefix,"tmp"))
  wdistance <- exact_transport_distance(obs, z, 1, 2)
  sink(NULL)
  if (wdistance > 1e10){
    wdistance <- 1e10
  }
  return(wdistance)
}


# #test
# set.seed(11)
# y_sim = target$simulat(target$rprior(1))
# compute_d(y_sim)

# thetas_ <- target$rprior(100, target$parameters)
# y_sim <- list()
# for (itheta in 1:nrow(thetas_)){
#   y_sim[[itheta]] <- target$simulate(thetas_[itheta,])
# }
#
# ds <- sapply(y_sim, compute_d)
# ds
#
# yy <- target$simulate(thetas_[1,])
# plot(yy[1,], yy[2,])
#
# plot(obs[1,], obs[2,])
# thetas_[1,]
# true_theta

param_algo <- list(nthetas = 2048, nmoves = 1, proposal = mixture_rmixmod(),
                   minimum_diversity = 0.5, R = 2, maxtrials = 100000)

savefile <- paste0(prefix, "mgandk.wsmc.n", nobservations, ".wasserstein.RData")
maxsimulation <- 2e6
results <- wsmc(compute_d, target, param_algo, savefile = savefile, maxsimulation = maxsimulation)
load(savefile)
# results$param_algo$maxtrials <- 10000
# results <- wsmc_continue(results, savefile = savefile, maxsimulation = 1e6)
# results$ncomputed %>% sum
# from_step <- 0
# wsmc.df <- wsmc_to_dataframe(results)
# index <- 7
# g <- ggplot(wsmc.df %>% filter(step > from_step), aes_string(x = target$parameter_names[index], colour = "time", group = "factor(time)")) + geom_density(aes(y = ..density..))
# g <- g + xlab(target$parameter_names[index]) # +  theme(legend.position = "none")
# g <- g + scale_color_gradient(name = "time (s)", low = rgb(1,0.5,0.5), high = "darkblue")
# g + geom_vline(xintercept = true_theta[index])
#
# qplot(x = cumsum(results$ncomputed), y = results$threshold_history, geom = "line") + scale_x_log10() + scale_y_log10()
# qplot(x = tail(results$distances_history, 1)[[1]], geom = "histogram")

