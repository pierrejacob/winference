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
was.thresholds = results_was$threshold_history
was.ncomp = cumsum(results_was$ncomputed)

# Swapping SMC
filename = paste0(prefix, "mgandk.wsmc.n", nobservations, ".swap.RData")
load(filename)
results_swap = results
swap.thresholds = results_swap$threshold_history
swap.ncomp = cumsum(results_swap$ncomputed)

# Hilbert SMC
filename = paste0(prefix, "mgandk.wsmc.n", nobservations, ".hilbert.RData")
load(filename)
results_hilbert = results
hilbert.thresholds = results_hilbert$threshold_history
hilbert.ncomp = cumsum(results_hilbert$ncomputed)

# MMD SMC
filename = paste0(prefix, "mgandk.wsmc.n", nobservations, ".mmd.RData")
load(filename)
results_mmd = results
mmd.thresholds = results_mmd$threshold_history
mmd.ncomp = cumsum(results_mmd$ncomputed)


w.df = data.frame(threshold = was.thresholds, ncomp = was.ncomp)
s.df = data.frame(threshold = swap.thresholds, ncomp = swap.ncomp)
h.df = data.frame(threshold = hilbert.thresholds, ncomp = hilbert.ncomp)
m.df = data.frame(threshold = mmd.thresholds, ncomp = mmd.ncomp)

g = ggplot(w.df, aes(x = ncomp, y = threshold)) + geom_line(aes(colour = "Wasserstein")) + geom_point(aes(colour = "Wasserstein"))
g = g + geom_line(data = s.df, aes(colour = "Swap")) + geom_point(data = s.df, aes(colour = "Swap"))
g = g + geom_line(data = h.df, aes(colour = "Hilbert")) + geom_point(data = h.df, aes(colour = "Hilbert"))
g = g + geom_line(data = m.df, aes(colour = "MMD")) + geom_point(data = m.df, aes(colour = "MMD"))
g = g + scale_color_manual(name = "", values = my_colors) + xlab("# model simulations") + ylab("threshold")
g = g + scale_y_log10(breaks = c(1e-2,1,1e2,1e4,1e6)) + scale_x_log10()
g = g + geom_label(data = data.frame(x = c(3e5,5e4,7e5,1e4), y = c(5e2,4e4,7,0.07), method = c("Wasserstein","Swap","Hilbert","MMD")),
                   aes(x = x, y = y, colour = method, label = method), size = 7) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "mgandk.n", nobservations, ".threshold_vs_ncomp.pdf"), plot = g, width = fig.width, height = fig.height)

