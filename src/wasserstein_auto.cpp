#include <RcppEigen.h>
#include <algorithm>
using namespace Rcpp;
using namespace std;
using namespace Eigen;

// [[Rcpp::export]]
List wasserstein_auto_(NumericVector p_, NumericVector q_, NumericMatrix cost_matrix_,
                  double epsilon, double desired_alpha){
  // compute distance between p and q
  // p corresponds to the weights of a N-sample
  // each q corresponds to the weights of a M-sample
  // Thus cost_matrix must be a N x M cost matrix
  // epsilon is a regularization parameter, equal to 1/lambda in some references
  int N = p_.size();
  int M = q_.size();
  
  Map<VectorXd> p(as<Map<VectorXd> >(p_));
  Map<VectorXd> q(as<Map<VectorXd> >(q_));
  Map<MatrixXd> cost_matrix(as<Map<MatrixXd> >(cost_matrix_));
  // avoid to take exp(k) when k is less than -500,
  // as K then contains zeros, and then the upcoming computations divide by zero
  MatrixXd K = (cost_matrix.array() * (-1./epsilon)); // K =  exp(- M / epsilon)
  for (int i = 0; i < N; i++){
    for (int j = 0; j < M; j++){
      if (K(i,j) < -500){
        K(i,j) = exp(-500);
      } else {
        K(i,j) = exp(K(i,j));
      }
    }
  }
  MatrixXd K_transpose = K.transpose();
  MatrixXd K_tilde = p.array().inverse().matrix().asDiagonal() * K; // diag(1/p) %*% K
  VectorXd u = VectorXd::Constant(N, 1./N);
  //
  VectorXd marginal1, marginal2;
  MatrixXd transportmatrix;
  VectorXd v;
  double alpha = 0;
  double beta = 0;
  int niterations_max = 1000;
  int iteration = 0;
  // for (int iteration = 0; iteration < niterations; iteration ++){
  while ((iteration < niterations_max) and (alpha < desired_alpha)){
    iteration ++;
    u = 1. / (K_tilde * (q.array() / (K_transpose * u).array()).matrix()).array();
    if (iteration % 10 == 1){
      // check if criterion is satisfied
      v = q.array() / (K_transpose * u).array();
      transportmatrix = u.col(0).asDiagonal() * K * v.col(0).asDiagonal();
      marginal1 = transportmatrix.rowwise().sum();
      marginal2 = transportmatrix.colwise().sum();
      alpha = 10;
      for (int i = 0; i < N; i++){
        beta = std::min(p(i) / marginal1(i), q(i) / marginal2(i));
        alpha = std::min(alpha, beta);
      }
      // cerr << "alpha = " << alpha << endl;
    }
  }
  v = q.array() / (K_transpose * u).array();
  // compute the optimal transport matrix between p and the first q
  transportmatrix = u.col(0).asDiagonal() * K * v.col(0).asDiagonal();
  // MatrixXd uXIv = u.array() * ((K.array() * cost_matrix.array()).matrix() * v).array();
  // NumericVector d = wrap(uXIv.colwise().sum());
  return List::create(Named("transportmatrix") = wrap(transportmatrix),
                      Named("u") = wrap(u),
                      Named("v") = wrap(v),
                      Named("iteration") = iteration);
}

