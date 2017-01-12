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

# number of observations
n = 100

#Other parameters
M = 500
N = c(1,5,10,20,50,100,500,1000)   #This is what we refer to as k in the paper
m = c(10,20,50,100,300,1000,5000,10000)

#Load results
load(paste0("~/Dropbox/ABCD/Results/gamma_optim/gamma_optim_data_n",n,"_interplay_k_m.RData"))

# ###Plotting parameters
# #width of plots
# w = 10
# #height of plots
# h = 5
# 
# #limits for axes
# yliml = 0
# ylimu = 1.5
# xliml = 1
# xlimu = 3
# 
# #legend position
# leg.pos = "none"
# 
# #fontsizes
# title.size = 10
# a.title.x = 10
# a.text.x = 8
# a.title.y = 10
# a.text.y = 8
# 
# 
# #Format the results
# df_mwek = data.frame(mu = mwe[1], sigma = mwe[2], k=1)
# df_mwem = data.frame(mu = mwe[1], sigma = mwe[2], m=1)
# 
# mewe = lapply(1:length(m),
#               function(i) lapply(1:length(N),
#                                  function(j) t(sapply(1:M,
#                                                       function(k) c(mewe_k_m[[k]][[1]][j,i],mewe_k_m[[k]][[2]][j,i])))))
# 
# df_mewe_list = list()
# for(i in 1:length(m)){
#   dummy = lapply(1:length(N), function(k) t(cbind(mewe[[i]][[k]],rep(k,M))))
#   dummy = matrix(unlist(dummy), ncol = 3, byrow = TRUE)
#   dummy = data.frame(dummy)
#   names(dummy) = c("mu","sigma","k")
#   df_mewe_list[[i]] = assign(paste("df_mewe_m",m[i],sep=""),dummy)
# }
# 
# ## Plot all levels of k for different m
# #On the same axes
# g = list()
# for(i in 1:length(m)){
#   g[[i]] <- ggplot(data = df_mewe_list[[i]], aes(x = mu, y = sigma, colour = k, group = k)) + ylim(0, 1.5) + xlim(1,3)
#   g[[i]] <- g[[i]] + geom_point(alpha = 0.5)
#   g[[i]] <- g[[i]] + scale_colour_gradient2(midpoint = 5) + theme(legend.position = leg.pos, plot.title = element_text(size=title.size))
#   g[[i]] <- g[[i]] + theme(axis.title.x=element_text(size=a.title.x), axis.text.x=element_text(size=a.text.x),
#                            axis.title.y=element_text(size=a.title.y), axis.text.y=element_text(size=a.text.y))
#   g[[i]] <- g[[i]] + theme(axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
#                            axis.title.y = element_text(margin = unit(c(0, 2, 0, 0), "mm")),
#                            plot.title = element_text(margin = unit(c(0, 2, 0, 0), "mm")))
#   g[[i]] <- g[[i]] + xlab(expression(mu)) + ylab(expression(sigma)) + ggtitle(paste("m = ",m[i],sep=""))
#   g[[i]] <- g[[i]] + geom_point(data = df_mwek, aes(x = mu, y = sigma), color="black")
# }
# ggsave(filename = "~/Dropbox/ABCD/draft4/gammadata_normalfit_n100_levelsk_diffm.pdf", plot = do.call(grid.arrange, c(g, ncol=4)),width = w, height = h)
# ggsave(filename = "~/Dropbox/ABCD/draft4/gammadata_normalfit_n100_levelsk_diffm.png", 
#        plot = do.call(grid.arrange, c(g, ncol=4)),width = w, height = h, dpi = 150)
# 
# #On scaled axes
# g = list()
# for(i in 1:length(m)){
#   g[[i]] <- ggplot(data = df_mewe_list[[i]] %>% filter(sigma> 0.2), aes(x = mu, y = sigma, colour = k, group = k)) #+ ylim(0, 1.5) + xlim(1,3)
#   g[[i]] <- g[[i]] + geom_point(alpha = 0.5)
#   g[[i]] <- g[[i]] + scale_colour_gradient2(midpoint = 5) + theme(legend.position = leg.pos, plot.title = element_text(size=title.size))
#   g[[i]] <- g[[i]] + theme(axis.title.x=element_text(size=a.title.x), axis.text.x=element_text(size=a.text.x),
#                            axis.title.y=element_text(size=a.title.y), axis.text.y=element_text(size=a.text.y))
#   g[[i]] <- g[[i]] + theme(axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
#                            axis.title.y = element_text(margin = unit(c(0, 2, 0, 0), "mm")),
#                            plot.title = element_text(margin = unit(c(0, 2, 0, 0), "mm")))
#   g[[i]] <- g[[i]] + xlab(expression(mu)) + ylab(expression(sigma)) + ggtitle(paste("m = ",m[i],sep=""))
#   g[[i]] <- g[[i]] + geom_point(data = df_mwek, aes(x = mu, y = sigma), color="black")
# }
# ggsave(filename = "~/Dropbox/ABCD/draft4/gammadata_normalfit_n100_scaledaxes_levelsk_diffm.pdf", plot = do.call(grid.arrange, c(g, ncol=4)),width = w, height = h)
# ggsave(filename = "~/Dropbox/ABCD/draft4/gammadata_normalfit_n100_scaledaxes_levelsk_diffm.png", plot = do.call(grid.arrange, c(g, ncol=4)),width = w, height = h, dpi = 150)
# 
# 
# 
# ### Stratify on k instead
# mewe_k = lapply(1:length(N),
#                 function(j) lapply(1:length(m),
#                                    function(i) t(sapply(1:M,
#                                                         function(k) c(mewe_k_m[[k]][[1]][j,i],mewe_k_m[[k]][[2]][j,i])))))
# 
# df_mewe_k_list = list()
# for(i in 1:length(N)){
#   dummy = lapply(1:length(m), function(k) t(cbind(mewe_k[[i]][[k]],rep(k,M))))
#   dummy = matrix(unlist(dummy), ncol = 3, byrow = TRUE)
#   dummy = data.frame(dummy)
#   names(dummy) = c("mu","sigma","m")
#   df_mewe_k_list[[i]] = assign(paste("df_mewe_k",N[i],sep=""),dummy)
# }
# 
# ### Plot all levels of m for different k
# #On the same axes
# g = list()
# for(i in 1:length(N)){
#   g[[i]] <- ggplot(data = df_mewe_k_list[[i]] %>% filter(sigma> 0.2), aes(x = mu, y = sigma, colour = m, group = m)) + ylim(yliml, ylimu) + xlim(xliml,xlimu)
#   g[[i]] <- g[[i]] + geom_point(alpha = 0.5)
#   g[[i]] <- g[[i]] + scale_colour_gradient2(midpoint = 5) + theme(legend.position = leg.pos, plot.title = element_text(size=title.size))
#   g[[i]] <- g[[i]] + theme(axis.title.x=element_text(size=a.title.x), axis.text.x=element_text(size=a.text.x),
#                            axis.title.y=element_text(size=a.title.y), axis.text.y=element_text(size=a.text.y))
#   g[[i]] <- g[[i]] + theme(axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
#                            axis.title.y = element_text(margin = unit(c(0, 2, 0, 0), "mm")),
#                            plot.title = element_text(margin = unit(c(0, 2, 0, 0), "mm")))
#   g[[i]] <- g[[i]] + xlab(expression(mu)) + ylab(expression(sigma)) + ggtitle(paste("k = ",N[i],sep=""))
#   g[[i]] <- g[[i]] + geom_point(data = df_mwem, aes(x = mu, y = sigma), color="black")
# }
# ggsave(filename = "~/Dropbox/ABCD/draft4/gammadata_normalfit_n100_levelsm_diffk.pdf", plot = do.call(grid.arrange, c(g, ncol=4)),width = w, height = h)
# ggsave(filename = "~/Dropbox/ABCD/draft4/gammadata_normalfit_n100_levelsm_diffk.png", plot = do.call(grid.arrange, c(g, ncol=4)),width = w, height = h, dpi = 150)
# 
# #On scaled axes
# g = list()
# for(i in 1:length(N)){
#   g[[i]] <- ggplot(data = df_mewe_k_list[[i]] %>% filter(sigma > 0.2), aes(x = mu, y = sigma, colour = m, group = m)) #+ ylim(0.2, 1.4) + xlim(1.2,2.6)
#   g[[i]] <- g[[i]] + geom_point(alpha = 0.5)
#   g[[i]] <- g[[i]] + scale_colour_gradient2(midpoint = 5) + theme(legend.position = leg.pos, plot.title = element_text(size=title.size))
#   g[[i]] <- g[[i]] + theme(axis.title.x=element_text(size=a.title.x), axis.text.x=element_text(size=a.text.x),
#                            axis.title.y=element_text(size=a.title.y), axis.text.y=element_text(size=a.text.y))
#   g[[i]] <- g[[i]] + theme(axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
#                            axis.title.y = element_text(margin = unit(c(0, 2, 0, 0), "mm")),
#                            plot.title = element_text(margin = unit(c(0, 2, 0, 0), "mm")))
#   g[[i]] <- g[[i]] + xlab(expression(mu)) + ylab(expression(sigma)) + ggtitle(paste("k = ",N[i],sep=""))
#   g[[i]] <- g[[i]] + geom_point(data = df_mwem, aes(x = mu, y = sigma), color="black")
# }
# #do.call(grid.arrange, c(g, ncol=4))
# ggsave(filename = "~/Dropbox/ABCD/draft4/gammadata_normalfit_n100_scaledaxes_levelsm_diffk.pdf", plot = do.call(grid.arrange, c(g, ncol=4)),width = w, height = h)
# ggsave(filename = "~/Dropbox/ABCD/draft4/gammadata_normalfit_n100_scaledaxes_levelsm_diffk.png", plot = do.call(grid.arrange, c(g, ncol=4)),width = w, height = h, dpi = 150)



## Plot subset of figures

###Plotting parameters
#width of plots
w = 5
#height of plots
h = 5
#legend position
leg.pos = "none"
#fontsizes
title.size = 10
a.title.x = 10
a.text.x = 8
a.title.y = 10
a.text.y = 8

#Format the results
df_mwek = data.frame(mu = mwe[1], sigma = mwe[2], k=1)
df_mwem = data.frame(mu = mwe[1], sigma = mwe[2], m=1)

mewe = lapply(1:length(m),
              function(i) lapply(1:length(N),
                                 function(j) t(sapply(1:M,
                                                      function(k) c(mewe_k_m[[k]][[1]][j,i],mewe_k_m[[k]][[2]][j,i])))))

df_mewe_list = list()
for(i in 1:length(m)){
  dummy = lapply(1:length(N), function(k) t(cbind(mewe[[i]][[k]],rep(k,M))))
  dummy = matrix(unlist(dummy), ncol = 3, byrow = TRUE)
  dummy = data.frame(dummy)
  names(dummy) = c("mu","sigma","k")
  df_mewe_list[[i]] = assign(paste("df_mewe_m",m[i],sep=""),dummy)
}


### Stratify on k instead
mewe_k = lapply(1:length(N),
                function(j) lapply(1:length(m),
                                   function(i) t(sapply(1:M,
                                                        function(k) c(mewe_k_m[[k]][[1]][j,i],mewe_k_m[[k]][[2]][j,i])))))

df_mewe_k_list = list()
for(i in 1:length(N)){
  dummy = lapply(1:length(m), function(k) t(cbind(mewe_k[[i]][[k]],rep(k,M))))
  dummy = matrix(unlist(dummy), ncol = 3, byrow = TRUE)
  dummy = data.frame(dummy)
  names(dummy) = c("mu","sigma","m")
  df_mewe_k_list[[i]] = assign(paste("df_mewe_k",N[i],sep=""),dummy)
}


N_sub = c(1,10,100,1000)
m_sub = c(10,100,1000,10000)

#Stratify on m
#On scaled axes
g = list()
for(i in 1:length(m_sub)){
  j = which(m==m_sub[i])
  g[[i]] <- ggplot(data = df_mewe_list[[j]] %>% filter(sigma> 0.2), aes(x = mu, y = sigma, colour = k, group = k)) #+ ylim(0, 1.5) + xlim(1,3)
  g[[i]] <- g[[i]] + geom_point(alpha = 0.5)
  g[[i]] <- g[[i]] + scale_colour_gradient2(midpoint = 5) + theme(legend.position = leg.pos, plot.title = element_text(size=title.size))
  g[[i]] <- g[[i]] + theme(axis.title.x=element_text(size=a.title.x), axis.text.x=element_text(size=a.text.x),
                           axis.title.y=element_text(size=a.title.y), axis.text.y=element_text(size=a.text.y))
  g[[i]] <- g[[i]] + theme(axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
                           axis.title.y = element_text(margin = unit(c(0, 2, 0, 0), "mm")),
                           plot.title = element_text(margin = unit(c(0, 2, 0, 0), "mm")))
  g[[i]] <- g[[i]] + xlab(expression(mu)) + ylab(expression(sigma)) + ggtitle(paste("m = ",m_sub[i],sep=""))
  g[[i]] <- g[[i]] + geom_point(data = df_mwek, aes(x = mu, y = sigma), color="black")
}
ggsave(filename = "~/Dropbox/ABCD/draft5/gammadata_normalfit_n100_scaledaxes_levelsk_diffm_subset.pdf", plot = do.call(grid.arrange, c(g, ncol=2)),width = w, height = h)
ggsave(filename = "~/Dropbox/ABCD/draft5/gammadata_normalfit_n100_scaledaxes_levelsk_diffm_subset.png", plot = do.call(grid.arrange, c(g, ncol=2)),width = w, height = h, dpi = 150)


#Stratify on k
g = list()
for(i in 1:length(N_sub)){
  j = which(N == N_sub[i])
  g[[i]] <- ggplot(data = df_mewe_k_list[[j]] %>% filter(sigma > 0.2), aes(x = mu, y = sigma, colour = m, group = m)) #+ ylim(0.2, 1.4) + xlim(1.2,2.6)
  g[[i]] <- g[[i]] + geom_point(alpha = 0.5)
  g[[i]] <- g[[i]] + scale_colour_gradient2(midpoint = 5) + theme(legend.position = leg.pos, plot.title = element_text(size=title.size))
  g[[i]] <- g[[i]] + theme(axis.title.x=element_text(size=a.title.x), axis.text.x=element_text(size=a.text.x),
                           axis.title.y=element_text(size=a.title.y), axis.text.y=element_text(size=a.text.y))
  g[[i]] <- g[[i]] + theme(axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
                           axis.title.y = element_text(margin = unit(c(0, 2, 0, 0), "mm")),
                           plot.title = element_text(margin = unit(c(0, 2, 0, 0), "mm")))
  g[[i]] <- g[[i]] + xlab(expression(mu)) + ylab(expression(sigma)) + ggtitle(paste("k = ",N_sub[i],sep=""))
  g[[i]] <- g[[i]] + geom_point(data = df_mwem, aes(x = mu, y = sigma), color="black")
}
ggsave(filename = "~/Dropbox/ABCD/draft5/gammadata_normalfit_n100_scaledaxes_levelsm_diffk_subset.pdf", plot = do.call(grid.arrange, c(g, ncol=2)),width = w, height = h)
ggsave(filename = "~/Dropbox/ABCD/draft5/gammadata_normalfit_n100_scaledaxes_levelsm_diffk_subset.png", plot = do.call(grid.arrange, c(g, ncol=2)),width = w, height = h, dpi = 150)

