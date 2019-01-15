library(winference)
library(doParallel)
library(doRNG)
library(dplyr)
library(ggthemes)

registerDoParallel(cores = 10)

rm(list = ls())
set.seed(1)
setmytheme()
my_colors <- get_my_colors()
my_colors <- c(my_colors, c("Sinkhorn" = "cornflowerblue"))
my_colors['Swap'] <- "009E73"
dimension <- 5
icomponent <- 2
target <- list()
target$parameters <- list()
# function to compute a dataset for each theta value
target$robservation <- function(nobservations, theta, parameters){
  mean <- rep(theta[1], dimension)
  covariance <- diag(theta[2], dimension, dimension)
  return(t(fast_rmvnorm(nobservations, mean, covariance)))
}
target$thetadim <- 2
target$ydim <- dimension
target$parameter_names <- c("mean", "var")

nobservations <- 5e2
w1 <- rep(1/nobservations, nobservations)
w2 <- rep(1/nobservations, nobservations)
target$simulate <- function(theta){
  target$robservation(nobservations, theta, target$parameters)
}

theta_dgp <- c(0,4)
y <- target$simulate(theta_dgp)
var(y[1,])
h_to_y <- get_hilbert_to_y(y, p = 2, ground_p = 2)
theta <- c(0, 4)
z <- target$simulate(theta)
#
C <- cost_matrix_L2(y, z)
hilbert <- hilbert_distance(y, z, p = 2, ground_p = 2)
exact <- as.numeric(exact_transport_given_C(w1, w2, C, p = 2))
swap <- swap_distance(y, z, p = 2, ground_p = 2, tolerance = 1e-5)$distance
cat(exact, swap, hilbert, "\n")
#
# compare computational costs
library(microbenchmark)
timings <- microbenchmark(generate = {z <- target$simulate(c(0,4));},
               hilbert = {z <- target$simulate(c(0,4)); hilbert <- hilbert_distance(y, z, p = 2, ground_p = 2);},
               swap = {z <- target$simulate(c(0,4)); swap <- swap_distance(y, z, p = 2, ground_p = 2, tolerance = 1e-5)$distance;},
               exact = {z <- target$simulate(c(0,4)); C <- cost_matrix_L2(y, z); exact <- as.numeric(exact_transport_given_C(w1, w2, C, p = 2))},
               times = 50)
timings

nthetas <- 100
thetas_ <- rbind(seq(from = -1, to  = 1, length.out = nthetas),
                seq(from = 0.1, to  = 3^2, length.out = nthetas))
thetas <- thetas_[icomponent,]


distances <- foreach(itheta = 1:nthetas, .combine = rbind) %dorng% {
  theta <- theta_dgp
  theta[icomponent] <- thetas[itheta]
  z <- target$simulate(theta)
  C <- cost_matrix_L2(y, z)
  hilb.distance <- hilbert_distance(y, z, p = 2, ground_p = 2)
  exact.distance <- as.numeric(exact_transport_given_C(w1, w2, C, p = 2))
  # sinkh1.distances[itheta] <- sinkhorn_given_C(w1, w2, C, p = 2, eps = 0.05, niterations = 100)$corrected
  # sinkh2.distances[itheta] <- sinkhorn_given_C(w1, w2, C, p = 2, eps = 0.025, niterations = 1000)$corrected
  swap.distance <- swap_distance(y, z, p = 2, ground_p = 2, tolerance = 1e-5)$distance
  c(exact.distance, swap.distance, hilb.distance)
}

exact.distances <- distances[,1]
swap.distances <- distances[,2]
hilb.distances <- distances[,3]

distances.df <- data.frame(thetas = thetas, h = hilb.distances, e = exact.distances, sw = swap.distances)
g <- ggplot(distances.df, aes(x = thetas, y = e)) + geom_point(aes(colour = "Wasserstein"))
g <- g + geom_point(aes(y = sw, colour = "Swap"))
g <- g + geom_point(aes(y = h, colour = "Hilbert"))
g <- g + scale_color_manual(name = "", values = my_colors) + xlab(expression(theta)) + ylab("distance") +
  geom_vline(xintercept = theta_dgp[icomponent], linetype = 2)
g <- g +  guides(colour = guide_legend(override.aes = list(size=10))) + theme(legend.text=element_text(size=30))
# g <- g + scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8))
if (icomponent==2){
  g <- g + geom_vline(xintercept = var(y[1,]), linetype = 3)
}
g

save(timings, distances.df, theta_dgp, y, dimension, icomponent, file = paste0("transportdistances.n", nobservations, ".d", dimension, ".i", icomponent, ".RData"))

# g <- ggplot(df, aes(x = thetas))
# g <- g + geom_point(aes(y = (sw-e)/e, colour = "Swap"))
# g <- g + geom_point(aes(y = (h-e)/e, colour = "Hilbert"))
# g <- g + scale_color_manual(name = "", values = my_colors) + xlab(expression(theta)) + ylab("distance") + ylim(0,1)
# g


# hilb.distances <- rep(0, nthetas)
# exact.distances <- rep(0, nthetas)
# sinkh1.distances <- rep(0, nthetas)
# sinkh2.distances <- rep(0, nthetas)
# swap.distances <- rep(0, nthetas)
# for (itheta in 1:nthetas){
#   theta <- theta_dgp
#   theta[icomponent] <- thetas[itheta]
#   z <- target$simulate(theta)
#   C <- cost_matrix_L2(y, z)
#   hilb.distances[itheta] <- hilbert_distance(y, z, p = 2, ground_p = 2)
#   exact.distances[itheta] <- as.numeric(exact_transport_given_C(w1, w2, C, p = 2))
#   sinkh1.distances[itheta] <- sinkhorn_given_C(w1, w2, C, p = 2, eps = 0.05, niterations = 100)$corrected
#   sinkh2.distances[itheta] <- sinkhorn_given_C(w1, w2, C, p = 2, eps = 0.025, niterations = 1000)$corrected
#   swap.distances[itheta] <- swap_distance(y, z, p = 2, ground_p = 2, tolerance = 1e-5)$distance
# }
# df <- data.frame(thetas = thetas, h = hilb.distances, e = exact.distances, s1 = sinkh1.distances, s2 = sinkh2.distances, sw = swap.distances)




# ggsave(filename = paste0("mvnorm.profile.component", icomponent, ".dimension", dimension, ".pdf"), plot = g, width=10, height = 7)

# g <- ggplot(df, aes(x = e, y = hilb.distances)) + geom_point(aes(colour = "Hilbert")) + geom_abline(slope = 1, intercept = 0)
# g <- g + geom_point(aes(y = sw, colour = "Swap")) + scale_color_manual(name = "", values = my_colors) + xlab("Wasserstein distance")
# g <- g + geom_point(aes(y = s2, colour = "Sinkhorn"))
# g <- g + ylab("distance")
# g <- g +  guides(colour = guide_legend(override.aes = list(size=10))) + theme(legend.text=element_text(size=30))
# g
# ggsave(filename = paste0("mvnorm.comparisonexact.component", icomponent, ".dimension", dimension, ".pdf"), plot = g, width=10, height = 7)




###########
# g <- ggplot(df, aes(x = thetas, y = s1)) + geom_point(aes(colour = "Sinkhorn 1"))
# g <- g + geom_point(aes(y = s2, colour = "Sinkhorn 2"))
# g <- g + geom_point(aes(y = e, colour = "Wasserstein"))
# g <- g + scale_color_colorblind(name = "") + xlab(expression(sigma)) + ylab("distance") +
#   geom_vline(xintercept = theta_dgp[icomponent], linetype = 2)
# if (dimension == 1){
#   g <- g + ylim(0, 2)
# }
# if (dimension == 2){
#   g <- g + ylim(0, 3)
# }
# if (dimension == 10){
#   g <- g + ylim(0, 10)
# }
# g
#
# g <- g + geom_point(aes(y = s1, colour = "Sinkhorn 1"))
# g <- g + geom_point(aes(y = s2, colour = "Sinkhorn 2"))
# g
#
# plot(thetas, hilb.distances, ylim = c(0, 10), col = "orange")
# points(thetas, exact.distances)
# points(thetas, sinkh1.distances, col = "red")
# points(thetas, sinkh2.distances, col = "red")
# points(thetas, swap.distances, col = "blue")
# abline(v = theta_dgp[icomponent])
#
# g <- ggplot(df, aes(x = e, y = hilb.distances)) + geom_point(aes(colour = "Hilbert")) + geom_abline(slope = 1, intercept = 0)
# g <- g + geom_point(aes(y = sw, colour = "Swapping"))
# g <- g + geom_point(aes(y = s1, colour = "Sinkhorn 1"))
# g <- g + geom_point(aes(y = s2, colour = "Sinkhorn 2")) + scale_color_colorblind(name = "") + xlab("Wasserstein distance")
# g <- g + ylab("Transport distance")
# g
