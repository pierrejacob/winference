rm(list = ls())

source("gandk_functions.R")
prefix = "output/"

true_theta = c(3,1,2,0.5)

M = 400
n = 1000  #Number of samples in each data set from dgp
N = 20    #Number of synthetic data sets for finding MEWE
m = 10^4  #Number of samples in the synthetic data sets
B = 1000  #Number of bootstrap samples per data set from dgp
alpha = 0.05 #Desired significance level, compute the alpha/2 and 1-alpha/2 quantiles
thetadim = 4 #Inference on all parameters

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
  filename = paste0(prefix,"gandk_bootstrap_run",i,"_N",N,"_m",m,"_n",n,"_B",B,".Rdata")
  load(filename)
  
  # Compute the empirical quantiles of the bootstrap samples,
  # both separately for each parameter and with Bonferroni correction.
  q_A = quantile(bs_store[,1],probs = c(alpha/8,alpha/2,1-alpha/2,1-alpha/8))
  q_B = quantile(bs_store[,2],probs = c(alpha/8,alpha/2,1-alpha/2,1-alpha/8))
  q_g = quantile(bs_store[,3],probs = c(alpha/8,alpha/2,1-alpha/2,1-alpha/8))
  q_k = quantile(bs_store[,4],probs = c(alpha/8,alpha/2,1-alpha/2,1-alpha/8))
  
  ci_A = q_A[2:3]
  ci_B = q_B[2:3]
  ci_g = q_g[2:3]
  ci_k = q_k[2:3]
  ci_all = rbind(q_A[c(1,4)],q_B[c(1,4)],q_g[c(1,4)],q_k[c(1,4)])
  
  #Check if inside
  store.cover[i,1] = is.inside(true_theta[1],ci_A)
  store.cover[i,2] = is.inside(true_theta[2],ci_B)
  store.cover[i,3] = is.inside(true_theta[3],ci_g)
  store.cover[i,4] = is.inside(true_theta[4],ci_k)
  
  if(is.inside(true_theta[1],ci_all[1,]) & is.inside(true_theta[2],ci_all[2,]) &
     is.inside(true_theta[3],ci_all[3,]) & is.inside(true_theta[4],ci_all[4,])){
    store.cover[i,5] = TRUE
  } else{
    store.cover[i,5] = FALSE
  }
}

coverages = apply(store.cover,2,mean)
