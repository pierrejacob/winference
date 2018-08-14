library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()
set.seed(11)
target <- get_gandk()

fig.height <- 5
fig.width <- 5

my_colors <- get_my_colors()

prefix <- ""

nobservations = 250
load(paste0(prefix, "gandkdata.RData"))
obs <- obs[1:nobservations]

# Wasserstein SMC
filename = paste0(prefix, "gandkwsmc.n", nobservations, ".RData")
load(filename)
results_was = results
wsmc.df = wsmc_to_dataframe(results_was)
nsteps = tail(wsmc.df$step,n=1)
ncomp = results_was$ncomputed

df = data.frame(steps = 1:nsteps, ncomp = ncomp)
g = ggplot(df, aes(x = steps, y = ncomp)) + geom_line(aes(colour = "Wasserstein")) + geom_point(aes(colour = "Wasserstein"))
g = g + scale_color_manual(name = "", values = my_colors) + xlab("steps") + ylab("# model simulations")
g = g + scale_y_log10(breaks = c(1e4,1e5,1e6,1e7,1e8)) + scale_x_continuous(breaks = c(10,20,40,30,50))
g <- g + geom_label(data = data.frame(x = c(34), y = c(6e7), method = c("Wasserstein")),
                    aes(x = x, y = y, colour = method, label = method), size = 7) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "gandk.n", nobservations, ".ncomp_vs_step.pdf"), plot = g, width = fig.width, height = fig.height)

