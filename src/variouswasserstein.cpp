#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
using namespace Eigen;

// // [[Rcpp::export]]
// List wasserstein_(NumericVector p_, NumericMatrix qs_, NumericMatrix cost_matrix_,
//                   double epsilon = 1, int niterations = 100){
//   // compute distance between one p and each q
//   // p corresponds to the weights of a N-sample
//   // each q corresponds to the weights of a M-sample
//   // Thus cost_matrix must be a N x M cost matrix
//   // epsilon is a regularization parameter, equal to 1/lambda in some references
//   int N = p_.size();
//   int M = qs_.rows();
//   int k = qs_.cols();
//
//   Eigen::Map<Eigen::MatrixXd> qs(as<Eigen::Map<Eigen::MatrixXd> >(qs_));
//   Eigen::Map<Eigen::MatrixXd> cost_matrix(as<Eigen::Map<Eigen::MatrixXd> >(cost_matrix_));
//   Eigen::Map<Eigen::VectorXd> p(as<Eigen::Map<Eigen::VectorXd> >(p_));
//
//   Eigen::MatrixXd K = (cost_matrix.array() * (-1./epsilon)); // K =  exp(- M / epsilon)
//   for (int i = 0; i < N; i++){
//     for (int j = 0; j < M; j++){
//       if (K(i,j) < -1000){
//         K(i,j) = exp(-1000);
//       } else {
//         K(i,j) = exp(K(i,j));
//       }
//     }
//   }
//   Eigen::MatrixXd K_transpose = K.transpose();
//   Eigen::MatrixXd K_tilde = p.array().inverse().matrix().asDiagonal() * K; // diag(1/p) %*% K
//   Eigen::MatrixXd u = Eigen::MatrixXd::Constant(N, k, 1./N);
//   for (int iteration = 0; iteration < niterations; iteration ++){
//     u = 1. / (K_tilde * (qs.array() / (K_transpose * u).array()).matrix()).array();
//   }
//   Eigen::MatrixXd v = qs.array() / (K_transpose * u).array();
//   // compute the optimal transport matrix between p and the first q
//   Eigen::MatrixXd transportmatrix = u.col(0).asDiagonal() * K * v.col(0).asDiagonal();
//   Eigen::MatrixXd uXIv = u.array() * ((K.array() * cost_matrix.array()).matrix() * v).array();
//   NumericVector d = wrap(uXIv.colwise().sum());
//   return List::create(Named("distances")=d,
//                       Named("transportmatrix") = wrap(transportmatrix),
//                       Named("u") = wrap(u),
//                       Named("v") = wrap(v));
// }

// [[Rcpp::export]]
List wasserstein_SAG_(NumericVector mu, NumericVector nu, NumericMatrix cost, double epsilon, int niterations, double stepsize){
  int nparticles = mu.size();
  NumericVector v(nparticles);
  NumericVector d(nparticles);
  NumericMatrix G(nparticles, nparticles);
  NumericVector grad(nparticles);
  int index;
  NumericVector uniform;
  double denominator;
  std::fill(v.begin(), v.end(), 0.);
  std::fill(d.begin(), d.end(), 0.);
  std::fill(G.begin(), G.end(), 0.);
  for (int iteration = 0; iteration < niterations; iteration++){
    uniform = runif(1);
    index = floor(uniform(0)*nparticles);
    d = d - G(index,_);
    denominator = sum(exp((v - cost(index,_)) / epsilon) * nu);
    grad = nu - exp((v - cost(index,_)) / epsilon) * nu / denominator;
    G(index,_) = mu(index) * grad;
    d = d + G(index,_);
    v = v + stepsize * d;
  }
  // compute u
  NumericVector u(nparticles);
  for (int i = 0; i < nparticles; i++){
    u(i) = - epsilon * log(sum(nu * exp((v - cost(i,_))/epsilon)));
  }
  NumericMatrix transportmatrix(nparticles, nparticles);
  for (int i = 0; i < nparticles; i++){
    for (int j = 0; j < nparticles; j++){
      transportmatrix(i,j) = exp((u(i) + v(j) - cost(i,j)) / epsilon) * mu(i) * nu(j);
    }
  }
  return List::create(Named("gradient") = wrap(d),
                      Named("transportmatrix") = wrap(transportmatrix),
                      Named("u") = wrap(u),
                      Named("v") = wrap(v));
}

// [[Rcpp::export]]
List wasserstein_SAG_auto_(NumericVector mu, NumericVector nu, NumericMatrix cost,
                           double epsilon, double tolerance, double stepsize, int maxiterations){
  int nparticles = mu.size();
  NumericVector v(nparticles);
  NumericVector d(nparticles);
  NumericMatrix G(nparticles, nparticles);
  NumericVector grad(nparticles);
  int index;
  NumericVector uniform;
  double denominator;
  std::fill(v.begin(), v.end(), 0.);
  std::fill(d.begin(), d.end(), 0.);
  std::fill(G.begin(), G.end(), 0.);
  int iteration = 0;
  double error = 10;
  while (error > tolerance){
    iteration++;
    uniform = runif(1);
    index = floor(uniform(0)*nparticles);
    d = d - G(index,_);
    denominator = sum(exp((v - cost(index,_)) / epsilon) * nu);
    grad = nu - exp((v - cost(index,_)) / epsilon) * nu / denominator;
    G(index,_) = mu(index) * grad;
    d = d + G(index,_);
    v = v + stepsize * d;
    error = (double) sum(abs(d));
  }
  NumericVector u(nparticles);
  for (int i = 0; i < nparticles; i++){
    u(i) = - epsilon * log(sum(nu * exp((v - cost(i,_))/epsilon)));
  }
  NumericMatrix transportmatrix(nparticles, nparticles);
  for (int i = 0; i < nparticles; i++){
    for (int j = 0; j < nparticles; j++){
      transportmatrix(i,j) = exp((u(i) + v(j) - cost(i,j)) / epsilon) * mu(i) * nu(j);
    }
  }
  return List::create(Named("gradient") = wrap(d),
                      Named("transportmatrix") = wrap(transportmatrix),
                      Named("u") = wrap(u),
                      Named("v") = wrap(v),
                      Named("iteration") = wrap(iteration));
}
