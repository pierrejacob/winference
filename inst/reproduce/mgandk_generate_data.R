library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()
set.seed(11)
target <- get_mgandk()
prefix = ""

nobservations <- 10000
true_theta <- c(3, 1, 1, 0.5, 4, 0.5, 2, 0.4, 0.6)

obs <- target$robservation(nobservations, true_theta)

save(true_theta, obs, file = paste0(prefix, "mgandkdata.RData"))
