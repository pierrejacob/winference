library(winference)
registerDoParallel(cores = 10)
rm(list = ls())
setmytheme()
set.seed(13)

# load g-and-k model
target <- get_gandk()
# theta^star
true_theta <- c(3, 1, 2, 0.5)
# number of observations; with a large value, we're approximating the asymptotic case
nobservations <- 500000
# generate y ~ mu_{theta^star}, so that we can approximate W(mu_theta, mu_theta^star) by W(z,y) where
# z ~ mu_theta
obs <- target$robservation(nobservations, true_theta,
                           target$parameters, target$generate_randomness(nobservations))
# hist(obs)
# get Wasserstein distance to y; in 1-dimension this is equivalent to sorting and is implemented in get_hilbert_to_y
compute_d2 <- get_hilbert_to_y(matrix(obs, nrow = 1), p = 2, ground_p = 2)


# takes additive noise to be added to true_theta (that is, theta^star)
# and compute, for each resulting theta,
# the distance to true_theta, as well as the associated wasserstein distance W(mu_theta, mu_theta^star)
f_ <- function(noise){
  nth <- dim(noise)[1]
  thetas <- t(apply(noise, 1, function(v) v + true_theta))
  thetadistances <- apply(thetas, 1, function(v) sqrt(sum((v - true_theta)^2)))
  wdistances <- foreach(itheta = 1:nth, .combine = c) %dorng% {
    theta <- thetas[itheta,]
    z <- target$robservation(nobservations, theta, target$parameters, target$generate_randomness(nobservations))
    compute_d2(matrix(z, nrow = 1))
  }
  return(list(thetadistances = thetadistances, wdistances = wdistances))
}


nthetas <- 100
# move only component 1
noise1 <- matrix(runif(4*nthetas, min = -0.3, max = .3), nrow = nthetas)
noise1[,c(2,3,4)] <- 0
res1 <- f_(noise1)
# move only component 2
noise2 <- matrix(runif(4*nthetas, min = -0.3, max = .3), nrow = nthetas)
noise2[,c(1,3,4)] <- 0
res2 <- f_(noise2)
# move only component 3
noise3 <- matrix(runif(4*nthetas, min = -0.3, max = .3), nrow = nthetas)
noise3[,c(1,2,4)] <- 0
res3 <- f_(noise3)
# move only component 4
noise4 <- matrix(runif(4*nthetas, min = -0.3, max = .3), nrow = nthetas)
noise4[,c(1,2,3)] <- 0
res4 <- f_(noise4)
#
# # move all components
nthetas <- 1000
noise_all <- matrix(runif(4*nthetas, min = -.3, max = .3), nrow = nthetas)
res_all <- f_(noise_all)
# #
save(res_all, res1, res2, res3, res4, file = "gandk.checkassumption.RData")

load(file = "gandk.checkassumption.RData")
# par(mfrow = c(1,1))
# plot(res_all$thetadistances, res_all$wdistances)
# plot(res_all$wdistances, res_all$thetadistances)
gall <- qplot(x = res_all$wdistances, y = res_all$thetadistances, geom = "point") + xlab("Wasserstein distance") + ylab("parameter distance")
gall
ggsave(filename = "gandk.check.thetaall.png", plot = gall, width = 5, height = 5)

# plot(abs(As - true_theta[1]), distances - min(distances))
par(mfrow = c(2,2))
# plot(res1$thetadistances, res1$wdistances, log = "y")
# plot(res2$thetadistances, res2$wdistances, log = "y")
# plot(res3$thetadistances, res3$wdistances, log = "y")
# plot(res4$thetadistances, res4$wdistances, log = "y")

g1 <- qplot(res1$wdistances, res1$thetadistances, geom = "point") + xlab("Wasserstein distance") + ylab("parameter distance")
g2 <- qplot(res2$wdistances, res2$thetadistances, geom = "point") + xlab("Wasserstein distance") + ylab("parameter distance")
g3 <- qplot(res3$wdistances, res3$thetadistances, geom = "point") + xlab("Wasserstein distance") + ylab("parameter distance")
g4 <- qplot(res4$wdistances, res4$thetadistances, geom = "point") + xlab("Wasserstein distance") + ylab("parameter distance")
library(gridExtra)
grid.arrange(g1, g2, g3, g4, nrow = 2)

ggsave(filename = "gandk.check.theta1.png", plot = g1, width = 5, height = 5)
ggsave(filename = "gandk.check.theta2.png", plot = g2, width = 5, height = 5)
ggsave(filename = "gandk.check.theta3.png", plot = g3, width = 5, height = 5)
ggsave(filename = "gandk.check.theta4.png", plot = g4, width = 5, height = 5)


#
# plot(res2$wdistances, res2$thetadistances)
# plot(res3$wdistances, res3$thetadistances)
# plot(res4$wdistances, res4$thetadistances)
#
