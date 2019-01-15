#Constructs the summary statistic in Bonassi
library(winference)
registerDoParallel(cores = detectCores())
rm(list = ls())
setmytheme()
set.seed(11)
prefix = ""

target <- get_toggleswitch()

#Parameters
p = 25000 #Number of data sets drawn from prior predictive
nobservations = 2000
R = 100

upper.bound = 2000
nbins = 50
breaks = seq(0, upper.bound, length.out = nbins)

target$simulate <- function(theta) matrix(target$robservation(nobservations, theta, target$parameters, target$generate_randomness(nobservations)), nrow = 1)

prior_draws = target$rprior(p, target$parameters)

filename = paste0(prefix,"toggle.binned.p",p,".nobs", nobservations, ".RData")

t = proc.time()
binned.data = foreach(i = 1:p) %dorng% {

  #simulate data set given draw from prior
  ysim = target$simulate(prior_draws[i,])

  #sets all values greater than R to R-0.00001
  ysim = sapply(ysim, function(x) min(x,upper.bound-0.00001))

  #bin the data
  binned = tabulate(findInterval(ysim, vec=breaks),nbins = nbins-1)

  return(list(binned,ysim))
}
t = proc.time() - t
save(binned.data,file = filename)

load(filename)

binned.data.matrix = t(sapply(binned.data, function(x) x[[1]]))

binned.data.svd = svd(binned.data.matrix)

first.col.A = binned.data.svd$u[,1]

#Choose references as the R percentiles of first column of A
reference.order = order(first.col.A)
reference.index = reference.order[c(1,(1:(R-1))*(p/R),p)]

#Set bandwith so as to mimic the histogram
bandwidth = upper.bound/nbins

#Continuous histograms
mids = breaks[-1] - (breaks[-1]-breaks[-nbins])/2

prob.mat = binned.data.matrix[reference.index,]/nobservations

#Check call
y = data.matrix[2,]
summary_full(y,bandwidth,prob.mat,mids)

data.matrix = t(sapply(binned.data, function(x) x[[2]]))

filename2 = paste0(prefix,"summary.full.p",p,".nobs", nobservations, ".RData")
t = proc.time()
summary.stat.full = foreach(i = 1:p) %dorng% {
  y = data.matrix[i,]
  out = summary_full(y,bandwidth,prob.mat,mids)
  return(out)
}
t = proc.time() - t
save(summary.stat.full,file = filename2)
load(filename2)

summary.full.matrix = t(sapply(summary.stat.full, function(x) x))
#replaces -Inf with -1000 (number on log scale)
summary.full.matrix = pmax(summary.full.matrix,-1000)

summary.full.svd = svd(summary.full.matrix)

#Number of principal components
r = 11
M = rbind(diag(1,r),matrix(0,R+1-r,r))
H = summary.full.svd$v%*%M
summary.reduced.matrix = summary.full.matrix%*%H
# filename3 = paste0("summary.H.p",p,".nobs", nobservations, ".RData")
# save(H,file = filename3)

summary.stat = function(y){
  S = summary_full(y,bandwidth,prob.mat,mids)
  S = pmax(S,-1000)
  S.red = S%*%H
  return(S.red)
}

# #Try
# t = proc.time()
# y = target$simulate(target$rprior(target$truetheta, target$parameters))
# t = proc.time()-t
# t
#
# t = proc.time()
# ss = summary.stat(y)
# t = proc.time()-t
# t
#
# t = proc.time()
# ss = sort(y)
# t = proc.time()-t
# t

