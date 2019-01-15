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

nobservations = 500
load(paste0(prefix, "mgandkdata.RData"))
obs <- obs[1:nobservations]

# Wasserstein SMC
filename = paste0(prefix, "mgandk.wsmc.n", nobservations, ".wasserstein.RData")
load(filename)
results_was = results
wsmc.df = wsmc_to_dataframe(results_was)
nsteps = tail(wsmc.df$step,n=1)
ncomp = results_was$ncomputed

# Swapping SMC
filename = paste0(prefix, "mgandk.wsmc.n", nobservations, ".swap.RData")
load(filename)
results_swap = results
swap.df = wsmc_to_dataframe(results_swap)
swap.nsteps = tail(swap.df$step,n=1)
swap.ncomp = results_swap$ncomputed

# Hilbert SMC
filename = paste0(prefix, "mgandk.wsmc.n", nobservations, ".hilbert.RData")
load(filename)
results_hilbert = results
hilbert.df = wsmc_to_dataframe(results_hilbert)
hilbert.nsteps = tail(hilbert.df$step,n=1)
hilbert.ncomp = results_hilbert$ncomputed

# MMD SMC
filename = paste0(prefix, "mgandk.wsmc.n", nobservations, ".mmd.RData")
load(filename)
results_mmd = results
mmd.df = wsmc_to_dataframe(results_mmd)
mmd.nsteps = tail(mmd.df$step,n=1)
mmd.ncomp = results_mmd$ncomputed


w.df = data.frame(steps = 1:nsteps, ncomp = ncomp)
s.df = data.frame(steps = 1:swap.nsteps, ncomp = swap.ncomp)
h.df = data.frame(steps = 1:hilbert.nsteps, ncomp = hilbert.ncomp)
m.df = data.frame(steps = 1:mmd.nsteps, ncomp = mmd.ncomp)

g = ggplot(w.df, aes(x = steps, y = ncomp)) + geom_line(aes(colour = "Wasserstein")) + geom_point(aes(colour = "Wasserstein"))
g = g + geom_line(data = s.df, aes(colour = "Swap")) + geom_point(data = s.df, aes(colour = "Swap"))
g = g + geom_line(data = h.df, aes(colour = "Hilbert")) + geom_point(data = h.df, aes(colour = "Hilbert"))
g = g + geom_line(data = m.df, aes(colour = "MMD")) + geom_point(data = m.df, aes(colour = "MMD"))
g = g + scale_color_manual(name = "", values = my_colors) + xlab("steps") + ylab("# model simulations")
g = g + scale_y_log10(breaks = c(0.5e4,1e4,2e4,4e4,1e5,2e5)) + scale_x_continuous(breaks = c(0,20,40,60,80))
g = g + geom_label(data = data.frame(x = c(60,70,60,17), y = c(7e3,1.5e4,1.7e5,1.9e4), method = c("Wasserstein","Swap","Hilbert","MMD")),
                    aes(x = x, y = y, colour = method, label = method), size = 7) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "mgandk.n", nobservations, ".ncomp_vs_step.pdf"), plot = g, width = fig.width, height = fig.height)
