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
# my_colors <- c(my_colors, c("Sinkhorn" = "cornflowerblue"))
# my_colors['Swap'] <- "009E73"
timings_all <- data.frame()
distances.df_all <- data.frame()

for (d in 2:5){
  load(paste0("transportdistances.n500.d", d, ".i2.RData"))
  timings_all <- rbind(timings_all, timings %>% group_by(expr) %>% summarise(m = median(time)/1e9) %>% mutate(dimension = d))
  distances.df$dimension <- d
  distances.df_all <- rbind(distances.df_all, distances.df)
}

timings_all
# timings %>% group_by(expr, dimension) %>% summarise(m = median(time)/1e9)
distances.df_all$dimension <- factor(distances.df_all$dimension, levels = 2:5, labels = paste0("d = ", 2:5))
head(distances.df_all)
g <- ggplot(distances.df_all, aes(x = thetas, y = e)) + geom_point(aes(colour = "Wasserstein"))
g <- g + geom_point(aes(y = sw, colour = "Swap"))
g <- g + geom_point(aes(y = h, colour = "Hilbert"))
g <- g + scale_color_manual(name = "", values = my_colors) + xlab(expression(sigma^2)) + ylab("distance") + facet_wrap(~ dimension, nrow = 1)
g <- g +  guides(colour = guide_legend(override.aes = list(size=5))) + theme(legend.text=element_text(size=20))
g <- g + scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8))
# if (icomponent==2){
#   g <- g + geom_vline(xintercept = var(y[1,]), linetype = 3) + geom_vline(xintercept = theta_dgp[icomponent], linetype = 2)
# }
g
ggsave(filename = "transportdistances.n500.i2.png", plot = g, width = 15, height = 5, dpi = 500)
