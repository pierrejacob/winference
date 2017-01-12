# R functions that call the cpp functions in wasserstein_semi_discrete.cpp.
# Computes W_p(p_hat,p_theta). Takes samples generated from p_theta as an argument.
# Uses c(x,y) = d(x,y)^p, where d is the Euclidean distance (L_2) and p can be specified. 
# Can change the cost function in the source file.

### Explantion of parameters 
# Give p_theta_samples as an n-by-d matrix, where n is the number of samples, d is the dimension of the space. 
# support denotes the support of the discrete distribution, give as m-by-d matrix.
# weights corresponding to the points in support, give as m-dim vector.
# choose p >= 1 for the p-Wasserstein distance.
# epsilon is the regularization parameter, set epsilon > 0.
# stepsize is a tuning parameter for the SGD, set stepsize > 0. 
# thresh is a threshold on the change in the vector v from the semi-dual formulation from Genevay et al. (2016).
# k determines how many samples are used for the naive estimation of the expectation in the semi-dual problem once --
# -- the vector v has been computed. uses the first k samples supplied in p_theta_samples

# sourceCpp("~/wasserstein_semi_discrete.cpp")

# Performs averaged SGD for a given number of iterations.
#'@export
wasserstein_p_semi_discrete = function(p_theta_samples, support, weights, p, epsilon, stepsize, k){
  v = dw_sgd_v_cpp(p_theta_samples,stepsize,weights,support,epsilon,p)
  distance = dw_est_cpp(p_theta_samples[1:k,],v,weights,support,epsilon,p)
  return(distance)
}

# Performes averaged SGD until a threshold on the change of v is reached.
#'@export
wasserstein_p_semi_discrete_thresh = function(p_theta_samples, support, weights, p, epsilon, stepsize, thresh, k){
  v = dw_sgd_v_thresh_cpp(p_theta_samples,stepsize,weights,support,epsilon,thresh,p)
  distance = dw_est_cpp(p_theta_samples[1:k,],v,weights,support,epsilon,p)
  return(distance)
}


# Approximate distance with a large number of samples and Sinkhorn's algorithm, uses the cost function specified in cost_matrix.cpp
#'@export
wasserstein_semi_discrete_sinkhorn = function(p_theta_samples, support, weights, epsilon, niterations){
  m = dim(p_theta_samples)[1]
  M = cost_matrix_L2(t(p_theta_samples), t(support))
  distance = wasserstein(rep(1/m,m),weights,M,epsilon,niterations)$distances
  return(distance)
}
