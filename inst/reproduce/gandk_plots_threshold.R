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
#step = 20
#thresholds = results_was$threshold_history[-(1:(step-1))]
step = 1
thresholds = results_was$threshold_history
ncomp = cumsum(results_was$ncomputed)

# g = qplot(x = step:(step+length(thresholds)-1), y = thresholds, geom = "line")
# g = g + xlab("step") + ylab("threshold")
# g = g + scale_y_log10(breaks = c(1e4,100,10,1,0.1,0.01))
# g
# ggsave(filename = paste0(prefix, "gandk.n", nobservations, ".threshold.pdf"), plot = g, width = fig.width, height = fig.height)

df = data.frame(thresholds = thresholds, ncomp = ncomp)
g = ggplot(df, aes(x = ncomp, y = thresholds)) + geom_line(aes(colour = "Wasserstein")) + geom_point(aes(colour = "Wasserstein"))
g = g + scale_color_manual(name = "", values = my_colors) + xlab("# model simulations") + ylab("threshold")
g = g + scale_y_log10(breaks = c(1e4,1e3,100,10,1,0.1,0.01)) + scale_x_log10(breaks = c(1e4,1e6,1e8))
g <- g + geom_label(data = data.frame(x = c(1e7), y = c(0.4), method = c("Wasserstein")),
                    aes(x = x, y = y, colour = method, label = method), size = 7) + theme(legend.position = "none")
g
ggsave(filename = paste0(prefix, "gandk.n", nobservations, ".threshold_vs_ncomp.pdf"), plot = g, width = fig.width, height = fig.height)

