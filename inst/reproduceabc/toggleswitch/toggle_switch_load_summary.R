#Loads the summary statistic provided the
prefix = ""

#Parameters
p = 25000 #Number of data sets drawn from prior predictive
nobservations = 2000
R = 100

upper.bound = 2000
nbins = 50
breaks = seq(0, upper.bound, length.out = nbins)

#Set bandwith so as to mimic the histogram
bandwidth = upper.bound/nbins
#mids from histogram
mids = breaks[-1] - (breaks[-1]-breaks[-nbins])/2

filename = paste0(prefix,"toggle.binned.p",p,".nobs", nobservations, ".RData")

load(filename)

binned.data.matrix = t(sapply(binned.data, function(x) x[[1]]))

binned.data.svd = svd(binned.data.matrix)

first.col.A = binned.data.svd$u[,1]

#Choose references as the R percentiles of first column of A
reference.order = order(first.col.A)
reference.index = reference.order[c(1,(1:(R-1))*(p/R),p)]

prob.mat = binned.data.matrix[reference.index,]/nobservations

# #Check call
# summary_full(y,bandwidth,prob.mat,mids)

data.matrix = t(sapply(binned.data, function(x) x[[2]]))

#Load full summary stat for all R references
filename2 = paste0(prefix,"summary.full.p",p,".nobs", nobservations, ".RData")
load(filename2)

summary.full.matrix = t(sapply(summary.stat.full, function(x) x))
#replace -Inf with -1000
summary.full.matrix = pmax(summary.full.matrix,-1000)

summary.full.svd = svd(summary.full.matrix)

#Number of principal components
r = 11
M = rbind(diag(1,r),matrix(0,R+1-r,r))
H = summary.full.svd$v%*%M

summary.stat = function(y){
  S = summary_full(y,bandwidth,prob.mat,mids)
  S = pmax(S,-1000)
  S.red = S%*%H
  return(S.red)
}
