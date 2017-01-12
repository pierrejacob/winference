#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
using namespace Eigen;

// [[Rcpp::export]]
List wasserstein_(NumericVector p_, NumericVector q_, NumericMatrix cost_matrix_,
                  double epsilon, int niterations){
  // compute distance between p and q
  // p corresponds to the weights of a N-sample
  // each q corresponds to the weights of a M-sample
  // Thus cost_matrix must be a N x M cost matrix
  // epsilon is a regularization parameter, equal to 1/lambda in some references
  int N = p_.size();
  int M = q_.size();
  
  Map<MatrixXd> q(as<Map<MatrixXd> >(q_));
  Map<MatrixXd> cost_matrix(as<Map<MatrixXd> >(cost_matrix_));
  Map<VectorXd> p(as<Map<VectorXd> >(p_));
  // avoid to take exp(k) when k is less than -500,
  // as K then contains zeros, and then the upcoming computations divide by zero
  MatrixXd K = (cost_matrix.array() * (-1./epsilon)).exp(); // K =  exp(- M / epsilon)
  
  // MatrixXd K = (cost_matrix.array() * (-1./epsilon)); // K =  exp(- M / epsilon)
  // for (int i = 0; i < N; i++){
    // for (int j = 0; j < M; j++){
  //     if (K(i,j) < -500){
  //       K(i,j) = exp(-500);
  //     } else {
        // K(i,j) = exp(K(i,j));
  //     }
    // }
  // }
  
  MatrixXd K_transpose = K.transpose();
  MatrixXd K_tilde = p.array().inverse().matrix().asDiagonal() * K; // diag(1/p) %*% K
  MatrixXd u = VectorXd::Constant(N, 1./N);
  for (int iteration = 0; iteration < niterations; iteration ++){
    // cerr << iteration << endl;
    //  u is set to 1 / (K_tilde %*% (qs / (K_transpose %*% u)))
    u = 1. / (K_tilde * (q.array() / (K_transpose * u).array()).matrix()).array();
    // for (int i = 0; i < N; i ++) cerr << u(i,0) << endl;
  }
  MatrixXd v = q.array() / (K_transpose * u).array();
  // compute the optimal transport matrix between p and the first q
  MatrixXd transportmatrix = u.col(0).asDiagonal() * K * v.col(0).asDiagonal();
  MatrixXd uXIv = u.array() * ((K.array() * cost_matrix.array()).matrix() * v).array();
  NumericVector d = wrap(uXIv.colwise().sum());
  return List::create(Named("distances")=d,
                      Named("transportmatrix") = wrap(transportmatrix),
                      Named("u") = wrap(u),
                      Named("v") = wrap(v));
}

