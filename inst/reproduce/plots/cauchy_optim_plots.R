#+ presets, echo = FALSE, warning = FALSE, message = FALSE
library(winference)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(foreach)
library(doMC)
library(doRNG)
library(reshape2)
library(gridExtra)
registerDoMC(cores = 4)
rm(list = ls())
setmytheme()
set.seed(11)

M = 1000
N = 20
m = 10^4
n = c(50,100,250,500,1000,5000,10000)


######## Plots used to generate draft4, 
#Load results
mewe1 = readRDS(paste0("~/Dropbox/ABCD/Results/cauchy_optim/mewe1_normal_ls_cauchy_M1000_N20_m10000.RDS"))

df_mewe1 = lapply(1:(length(n)), function(k) t(cbind(mewe1[[k]],rep(k,M))))
df_mewe1 = matrix(unlist(df_mewe1), ncol = 3, byrow = TRUE)
df_mewe1 = data.frame(df_mewe1)
names(df_mewe1) = c("mu","sigma","n")

df_mewe1_scaled = lapply(1:(length(n)), function(k) t(cbind(sqrt(n[k])*mewe1[[k]],rep(k,M))))
df_mewe1_scaled = matrix(unlist(df_mewe1_scaled), ncol = 3, byrow = TRUE)
df_mewe1_scaled = data.frame(df_mewe1_scaled)
names(df_mewe1_scaled) = c("mu","sigma","n")


#width of plots
w = 7
#height of plots
h = 5

#limits for axes
yliml = 0
ylimu = 1.5
xliml = 1
xlimu = 3

#legend position
leg.pos = "none"

#fontsizes
title.size = 14
a.title.x = 12
a.text.x = 12
a.title.y = 12
a.text.y = 12

#Scatterplot
g <- ggplot(data = df_mewe1, aes(x = mu, y = sigma, colour = n, group = n)) #+ ylim(0, 1.5) + xlim(1,3)
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = 4) + theme(legend.position = leg.pos, plot.title = element_text(size=title.size))
g <- g + xlab(expression(mu)) + ylab(expression(sigma)) 
g
ggsave(filename = "~/Dropbox/ABCD/draft4/cauchydata_normalfit_scatter_M1000_k20_m10000.png", 
       plot = g ,width = w, height = h, dpi = 150)


#Histograms of mu
g <- ggplot(data = df_mewe1_scaled, aes(x = mu, group = n)) #+ ylim(0, 1.5) + xlim(1,3)
g <- g + geom_density(aes(y=..density..))
g <- g + xlab(expression(mu))
g
ggsave(filename = "~/Dropbox/ABCD/draft4/cauchydata_normalfit_density_M1000_k20_m10000.pdf", 
       plot = g ,width = w, height = h)

