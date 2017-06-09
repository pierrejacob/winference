library(winference)
registerDoParallel(cores = 6)
rm(list = ls())
setmytheme()

set.seed(13)
# load model
target <- get_levydriven()
# prefix for file names
prefix <- ""
# load data
load(file = paste0(prefix, "levydrivendata.RData"))
# number of observations
nobservations <- 10000
# lag value
lagvalue <- 1
# subset observations
obs <- obs[1:nobservations]
# plot observations
qplot(x = 1:nobservations, y = obs, geom = "line") + xlab("time") + ylab(expression(y))

filename <- paste0(prefix, "levydriven.n", nobservations, ".lag", lagvalue, ".wsmc.hilbert.RData")
load(file = filename)
results1 <- results
results1$ncomputed %>% sum

g <- qplot(x = cumsum(results1$ncomputed), y = results1$threshold_history, geom = "line") + scale_y_log10() + scale_x_log10()
g <- g + xlab("# distances calculated") + ylab("threshold")
g

nsteps1 <- length(results1$thetas_history)
#

wsmc.df <- wsmc_to_dataframe(results1)

## plots of posterior after using only Hilbert
g <- ggplot(wsmc.df, aes(x = mu, y = beta, colour = step, group = step)) + geom_point(alpha = 0.5)
g <- g + theme(legend.position = "none") + xlab(expression(mu)) + ylab(expression(beta))
g <- g + scale_colour_gradient2(midpoint = floor(nsteps1/2))
g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
# g
ggsave(filename = paste0(prefix, "levydriven_mubeta.png"), plot = g, width = 7, height = 5, dpi = 150)

g <- ggplot(wsmc.df, aes(x = xi, y = omega2, colour = step, group = step)) + geom_point(alpha = 0.5)
g <- g + theme(legend.position = "none") + xlab(expression(xi)) + ylab(expression(omega^2))
g <- g + scale_colour_gradient2(midpoint = floor(nsteps1/2))
g <- g + geom_vline(xintercept = true_theta[3]) + geom_hline(yintercept = true_theta[4])
# g
ggsave(filename = paste0(prefix, "levydriven_xiomega2.png"), plot = g, width = 7, height = 5, dpi = 150)

g <- ggplot(wsmc.df, aes(x = lambda, colour = step, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(lambda))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
g <- g + geom_vline(xintercept = true_theta[5])
# g
ggsave(filename = paste0(prefix, "levydriven_lambda.pdf"), plot = g, width = 7, height = 5)
#

# plot of distances
threshold1 <- tail(results1$threshold_history, 2)[1]
qplot(x = tail(results1$distances_history, 1)[[1]], geom = "blank") + geom_histogram() + geom_vline(xintercept = threshold1)


#
f <- function(dataset){
  return(sum(acf(dataset^2, plot = F, lag.max = 50)$acf[,,1][1:50]))
}
f_obs <- f(obs)

fs <- as.numeric(foreach(i = 1:results1$param_algo$nthetas, .combine = c) %dorng% {
  f(results1$latest_y[[i]])
})

g <- qplot(x = results1$thetas_history[[nsteps1]][,5], y = fs, geom = "point") + geom_hline(yintercept = f_obs, linetype = 2)
g <- g + xlab(expression(lambda)) + ylab("ACF summary")
g
ggsave(filename = paste0(prefix, "levydriven_acfsummary.png"), plot = g, width = 7, height = 5, dpi = 150)

# does not match at all
filename2 <- paste0(prefix, "levydriven.n", nobservations, ".lag", lagvalue, ".wsmc.hilbert.summary.RData")
load(file = filename2)
results2 <- results
nsteps2 <- length(results2$thetas_history)
#
wsmc2.df <- wsmc_to_dataframe(results2)

g <- ggplot(wsmc2.df %>% filter(step > nsteps1), aes(x = mu, y = beta, colour = step, group = step)) + geom_point(alpha = 0.5)
g <- g + theme(legend.position = "none") + xlab(expression(mu)) + ylab(expression(beta))
g <- g + scale_colour_gradient2(midpoint = floor(nsteps1 + (nsteps2-nsteps1)/2))
g <- g + geom_vline(xintercept = true_theta[1]) + geom_hline(yintercept = true_theta[2])
g

g <- ggplot(wsmc2.df %>% filter(step > nsteps1), aes(x = xi, y = omega2, colour = step, group = step)) + geom_point(alpha = 0.5)
g <- g + theme(legend.position = "none") + xlab(expression(xi)) + ylab(expression(omega^2))
g <- g + scale_colour_gradient2(midpoint = floor(nsteps1 + (nsteps2-nsteps1)/2))
g <- g + geom_vline(xintercept = true_theta[3]) + geom_hline(yintercept = true_theta[4])
g

g <- ggplot(wsmc2.df %>% filter(step > nsteps1), aes(x = omega2, y = lambda, colour = step, group = step)) + geom_point(alpha = 0.5)
g <- g + theme(legend.position = "none") + xlab(expression(omega^2)) + ylab(expression(lambda))
g <- g + scale_colour_gradient2(midpoint = floor(nsteps1 + (nsteps2-nsteps1)/2))
g <- g + geom_vline(xintercept = true_theta[4]) + geom_hline(yintercept = true_theta[5])
g
ggsave(filename = paste0(prefix, "levydriven_omega2lambda.png"), plot = g, width = 7, height = 5, dpi = 150)

g <- ggplot(wsmc2.df  %>% filter(step > nsteps1), aes(x = omega2, colour = step, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(omega^2)) + geom_vline(xintercept = true_theta[4])
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
g

g <- ggplot(wsmc2.df  %>% filter(step > nsteps1), aes(x = lambda, colour = step, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(lambda)) + geom_vline(xintercept = true_theta[5])
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue")
g

# final posterior
g <- ggplot(wsmc2.df  %>% filter(step == nsteps2), aes(x = mu, colour = step, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(mu))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue") + geom_vline(xintercept = true_theta[1])
g

g <- ggplot(wsmc2.df  %>% filter(step == nsteps2), aes(x = beta, colour = step, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(beta))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue") + geom_vline(xintercept = true_theta[2])
g

g <- ggplot(wsmc2.df  %>% filter(step == nsteps2), aes(x = xi, colour = step, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(xi))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue") + geom_vline(xintercept = true_theta[3])
g

g <- ggplot(wsmc2.df  %>% filter(step > nsteps1), aes(x = omega2, colour = step, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(omega^2))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue") + geom_vline(xintercept = true_theta[4])
g
ggsave(filename = paste0(prefix, "levydriven_omega2_summar.pdf"), plot = g, width = 7, height = 5)

#
g <- ggplot(wsmc2.df  %>% filter(step > nsteps1), aes(x = lambda, colour = step, group = step)) + geom_density(aes(y = ..density..))
g <- g + theme(legend.position = "none") + xlab(expression(lambda))
g <- g + scale_color_gradient(low = rgb(1,0.5,0.5), high = "darkblue") + geom_vline(xintercept = true_theta[5])
g <- g + scale_x_log10(breaks = c(0.001,0.01, 0.1, 1))
g
ggsave(filename = paste0(prefix, "levydriven_lambda_summar.pdf"), plot = g, width = 7, height = 5)
