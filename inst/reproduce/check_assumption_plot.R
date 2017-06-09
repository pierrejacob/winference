library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()
set.seed(13)


dgpname <- "gamma"
modelname <- "normal"

# dgp
if (dgpname == "gamma"){
   robs <- function(n) rgamma(n, 10, 5)
} else {
  robs <- function(n) rnorm(n, mean = 2, sd = 1)
}
nobs <- 500000

# model
if (modelname == "gamma"){
  rmodel <- function(n, theta) rgamma(n, theta, 5)
} else {
  rmodel <- function(n, theta) rnorm(n, mean = theta, sd = 1)
}
prefix <- ""
savefile <- paste0(prefix, "check.as6.dgp", dgpname, ".model", modelname, ".RData")

nseq <- 100
if (modelname == "gamma"){
  theta_seq <- seq(from = 9.5, to = 10.5, length.out = nseq)
} else {
  theta_seq <- seq(from = 1.5, to = 2.5, length.out = nseq)
}

# Wasserstein distance with p = 2
y_obs <- robs(nobs)
compute_d2 <- get_hilbert_to_y(matrix(y_obs, nrow = 1), p = 2, ground_p = 2)
results <- foreach(i = 1:nseq) %dorng% {

  theta <- theta_seq[i]
  y_fake <- rmodel(nobs, theta)
  # hilbert_distance(matrix(y_obs, nrow = 1), matrix(y_fake, nrow = 1))
  compute_d2(matrix(y_fake, nrow = 1))
}
save(results, file = savefile)
load(file = savefile)

wp2 <- sapply(results, function(v) v)

theta_star <- theta_seq[which.min(wp2)]
x <- abs(theta_seq - 2)
ysqrt <- sqrt(abs(wp2 - wp2[which.min(wp2)]))
y <- abs(wp2 - wp2[which.min(wp2)])

figurefile <- paste0(prefix, "check.as6.dgp", dgpname, ".model", modelname, ".theta-vs-wasserstein.pdf")
g <- qplot(x = theta_seq, y = wp2, geom = "point") + xlab(expression(theta)) + ylab("Wasserstein distance")
print(g)
ggsave(filename = figurefile, plot = g, width = 7, height = 5)

if (modelname == dgpname){
  figurefile <- paste0(prefix, "check.as6.dgp", dgpname, ".model", modelname, ".distancetheta.pdf")
  g <- qplot(x = x, y = y, geom = "point") + xlab("distance to theta star") + ylab("Wasserstein distance") + geom_smooth(method = "lm", se = FALSE)
  print(g)
  ggsave(filename = figurefile, plot = g, width = 7, height = 5)
} else {
  figurefile <- paste0(prefix, "check.as6.dgp", dgpname, ".model", modelname, ".distancetheta.pdf")
  g <- qplot(x = x, y = ysqrt, geom = "point") + geom_smooth(method = "lm", se = FALSE) + xlab(expression(abs(theta-2))) + ylab("sqrt(Wasserstein distance)")
  print(g)
  ggsave(filename = figurefile, plot = g, width = 7, height = 5)
}


