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

mle = readRDS("~/Dropbox/ABCD/Results/gamma_concentration/mle_normal_ls_gamma_10_5_M1000.RDS")
mewe1 = readRDS("~/Dropbox/ABCD/Results/gamma_concentration/mewe_normal_ls_gamma_10_5_M1000_N20_m10000.RDS")

df_mewe1 = lapply(1:(length(n)), function(k) t(cbind(mewe1[[k]],rep(k,M))))
df_mewe1 = matrix(unlist(df_mewe1), ncol = 3, byrow = TRUE)
df_mewe1 = data.frame(df_mewe1)
names(df_mewe1) = c("mu","sigma","n")


df_mle = lapply(1:(length(n)), function(k) t(cbind(mle[[k]],rep(k,M))))
df_mle = matrix(unlist(df_mle), ncol = 3, byrow = TRUE)
df_mle = data.frame(df_mle)
names(df_mle) = c("mu","sigma","n")


#width of plots
w = 5
#height of plots
h = 5
#legend position
leg.pos = "none"

#Scatterplot
g <- ggplot(data = df_mewe1, aes(x = mu, y = sigma, colour = n, group = n)) +ylim(0.4, 0.95) + xlim(1.65,2.35)
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = 4) + theme(legend.position = leg.pos, plot.title = element_text(size=title.size))
g <- g + xlab(expression(mu)) + ylab(expression(sigma))
g
ggsave(filename = "~/Dropbox/ABCD/draft4/gammadata_normalfit_scatter_M1000_k20_m10000.png", 
       plot = g ,width = w, height = h, dpi = 150)


#Scatterplot
g <- ggplot(data = df_mle, aes(x = mu, y = sigma, colour = n, group = n)) + ylim(0.4, 0.95) + xlim(1.65,2.35)
g <- g + geom_point(alpha = 0.5)
g <- g + scale_colour_gradient2(midpoint = 4) + theme(legend.position = leg.pos, plot.title = element_text(size=title.size))
g <- g + xlab(expression(mu)) + ylab(expression(sigma))
g
ggsave(filename = "~/Dropbox/ABCD/draft4/gammadata_normalfit_scatterMLE_M1000.png", 
       plot = g ,width = w, height = h, dpi = 150)

#Histograms of mu
g <- ggplot(data = df_mewe1 %>% filter(n==7), aes(x = mu, group = n)) #+ ylim(0, 1.5) + xlim(1,3)
g <- g + geom_density(data = df_mle %>% filter(n==7), linetype = 2)
g <- g + geom_density(aes(y=..density..))
g <- g + xlab(expression(mu))
g
ggsave(filename = "~/Dropbox/ABCD/draft4/gammadata_normalfit_muhist_M1000_k20_m10000.pdf", 
       plot = g ,width = w, height = h)

#Histograms of sigma
g <- ggplot(data = df_mewe1 %>% filter(n==7), aes(x = sigma, group = n)) #+ ylim(0, 1.5) + xlim(1,3)
g <- g + geom_density(data = df_mle %>% filter(n==7), linetype = 2)
g <- g + geom_density(aes(y=..density..))
g <- g + xlab(expression(sigma))
g
ggsave(filename = "~/Dropbox/ABCD/draft4/gammadata_normalfit_sigmahist_M1000_k20_m10000.pdf", 
       plot = g ,width = w, height = h)

