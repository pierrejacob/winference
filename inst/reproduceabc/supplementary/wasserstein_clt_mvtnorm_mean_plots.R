library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()
set.seed(11)

fig.height <- 5
fig.width <- 5

my_colors <- get_my_colors()

dset = c(2,5,10)
nset = c(10,100,1000)

for(dd in 1:3){
  for(nn in 1:3){

    d = dset[dd]

    #obs
    n = nset[nn]
    obs = fast_rmvnorm(n,rep(2,d),diag(1,d))
    obs_mean = apply(obs,2,mean)

    #fake
    m = nset[nn]

    compute_d = function(fake_obs){
      swap_distance(t(obs), t(fake_obs), tolerance = 1e-5)$distance
    }
    #test
    fake_obs = fast_rmvnorm(n,rep(1,d),diag(1,d))
    compute_d(fake_obs)



    #Grid
    grid = seq(obs_mean[1]-1.5,obs_mean[1]+1.5,length.out = 20)

    #Number of data sets per grid point
    N = 100

    #Variance of fake data sets
    sigma2 = 1

    # dists = foreach(k = 1:length(grid), .combine = rbind) %dorng% {
    #   store_distances = rep(0,N)
    #   store_meandiffs = rep(0,N)
    #   for(i in 1:N){
    #     #fake_obs = fast_rmvnorm(m,c(grid[k],obs_mean[2:d]),diag(sigma2,d))
    #     fake_obs = fast_rmvnorm(m,c(grid[k],rep(1,d-1)),diag(sigma2,d))
    #     fake_mean = apply(fake_obs,2,mean)
    #     store_distances[i] = compute_d(fake_obs)
    #     store_meandiffs[i] = sum((obs_mean-fake_mean)^2)^0.5
    #   }
    #   out = cbind(rep(grid[k],N),store_distances,store_meandiffs)
    #   return(out)
    # }
    #
    # dists.df = data.frame(theta = dists[,1], Wasserstein = dists[,2], xbardiff = dists[,3])
    # save(dists.df,file = paste0(prefix,"was_vs_xbar_mvtnorm_d",d,"_n",n,"_var",sigma2,".Rdata"))
    load(paste0(prefix,"was_vs_xbar_mvtnorm_d",d,"_n",n,"_var",sigma2,".Rdata"))

    wmax = max(dists.df$Wasserstein)
    wmin = min(dists.df$Wasserstein)
    wmax - wmin

    xmax = max(dists.df$xbardiff)
    xmin = min(dists.df$xbardiff)
    xmax - xmin

    ll = max(wmax-wmin,xmax-xmin)

    # g = ggplot(dists.df, aes(theta, Wasserstein, group = theta)) + geom_boxplot()
    # g = g + coord_cartesian(ylim = c(wmin-ll/20,wmin+ll))
    # g
    #
    # g = ggplot(dists.df, aes(theta, xbardiff, group = theta)) + geom_boxplot()
    # g = g + coord_cartesian(ylim = c(xmin-ll/20,xmin+ll))
    # g

    g = ggplot(data = dists.df, aes(theta, Wasserstein, group = theta)) + geom_boxplot(aes(fill = "Wasserstein"),alpha=0.7)
    g = g + geom_boxplot(data = dists.df, aes(theta, xbardiff, group = theta, fill = "Summary"),alpha=0.7)
    g = g + scale_fill_manual(name = "", values = my_colors) + xlab(expression(theta[1])) + ylab("distance")
    #g
    ggsave(filename = paste0(prefix,"was_vs_xbar_mvtnorm_d",d,"_n",n,"_var",sigma2,".pdf"), plot = g, width = fig.width, height = fig.height)
  }
}
