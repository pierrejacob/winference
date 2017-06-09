library(winference)
registerDoParallel(cores = 8)
rm(list = ls())
setmytheme()
set.seed(11)
target <- get_mgandk()
prefix <- ""

nobservations <- 500
load(paste0(prefix, "mgandkdata.RData"))
obs <- obs[,1:nobservations]
target$simulate <- function(theta) target$robservation(nobservations, theta)

data_sets <- list()

x <- proc.time()
target$loglikelihood(target$rprior(1000, target$parameters), obs)
newx <- proc.time()
elapsed <- (newx - x)[3]
cat("time to compute likelihood:", elapsed/nrep, "\n")

nrep <- 1000
x <- proc.time()
for (irep in 1:nrep){
  data_sets[[irep]] <- target$simulate(true_theta)
}
newx <- proc.time()
elapsed <- (newx - x)[3]
cat("time to simulate data sets:", elapsed/nrep, "\n")
x <- proc.time()
for (irep in 1:nrep){
  sinkhorn_distance(obs, data_sets[[irep]])
}
newx <- proc.time()
elapsed <- (newx - x)[3]
cat("time to compute Sinkhorn distance:", elapsed/nrep, "\n")

x <- proc.time()
for (irep in 1:nrep){
  exact_transport_distance(obs, data_sets[[irep]])
}
newx <- proc.time()
elapsed <- (newx - x)[3]
cat("time to compute Wasserstein distance:", elapsed/nrep, "\n")

x <- proc.time()
for (irep in 1:nrep){
  hilbert_distance(obs, data_sets[[irep]])
}
newx <- proc.time()
elapsed <- (newx - x)[3]
cat("time to compute Hilbert distance:", elapsed/nrep, "\n")


cd <- get_mmd_to_y(obs)
x <- proc.time()
for (irep in 1:nrep){
  cd(data_sets[[irep]])
}
newx <- proc.time()
elapsed <- (newx - x)[3]
cat("time to compute MMD distance:", elapsed/nrep, "\n")

x <- proc.time()
for (irep in 1:nrep){
  swap_distance(obs, data_sets[[irep]], tolerance = 1e-5)
}
newx <- proc.time()
elapsed <- (newx - x)[3]
cat("time to compute Swap distance:", elapsed/nrep, "\n")


