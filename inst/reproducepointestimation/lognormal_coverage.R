## post process output of lognormal.cluster.script.R
## to compute coverage of confidence intervals
rm(list = ls())

prefix = "output/"

true_theta = c(0,1)

M = 400
n = 1000  #Number of samples in each data set from dgp
N = 20    #Number of synthetic data sets for finding MEWE
m = 10^4  #Number of samples in the synthetic data sets
B = 1000  #Number of bootstrap samples per data set from dgp
alpha = 0.05 #Desired significance level, compute the alpha/2 and 1-alpha/2 quantiles
thetadim = 2 #Inference on all parameters

#Check if value is inside the interval
is.inside = function(value,interval){
  if(interval[1] <= value & value <= interval[2]){
    return(TRUE)
  } else{
    return(FALSE)
  }
}

#Store the number of intervals that cover
store.cover = matrix(NA,nrow = M, ncol = thetadim + 1)

for(i in 1:M){
  
  # Load file from run i
  filename = paste0(prefix,"lognormal_bootstrap_run",i,"_N",N,"_m",m,"_n",n,"_B",B,".Rdata")
  load(filename)
  
  # Compute the empirical quantiles of the bootstrap samples,
  # both separately for each parameter and with Bonferroni correction.
  q_mu = quantile(bs_store[,1],probs = c(alpha/4,alpha/2,1-alpha/2,1-alpha/4))
  q_sigma = quantile(bs_store[,2],probs = c(alpha/4,alpha/2,1-alpha/2,1-alpha/4))
  
  ci_mu = q_mu[2:3]
  ci_sigma = q_sigma[2:3]
  ci_all = rbind(q_mu[c(1,4)],q_sigma[c(1,4)])
  
  #Check if inside
  store.cover[i,1] = is.inside(true_theta[1],ci_mu)
  store.cover[i,2] = is.inside(true_theta[2],ci_sigma)
  
  if(is.inside(true_theta[1],ci_all[1,]) & is.inside(true_theta[2],ci_all[2,])){
    store.cover[i,3] = TRUE
  } else{
    store.cover[i,3] = FALSE
  }
}

coverages = apply(store.cover,2,mean)
