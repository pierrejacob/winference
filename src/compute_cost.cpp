#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
using namespace Eigen;

// [[Rcpp::export]]
NumericMatrix compute_cost1_(const NumericVector & x, const NumericVector & y){
  NumericMatrix cost(x.size(), y.size());
  for (int ix = 0; ix < x.size(); ix ++){
    for (int iy = 0; iy < y.size(); iy ++){
      cost(ix,iy) = std::abs(x(ix)-y(iy));
    }
  }
  return cost;
}


// [[Rcpp::export]]
NumericMatrix compute_cost2_(const NumericVector & x, const NumericVector & y){
  NumericMatrix cost(x.size(), y.size());
  for (int ix = 0; ix < x.size(); ix ++){
    for (int iy = 0; iy < y.size(); iy ++){
      cost(ix,iy) = std::pow(x(ix)-y(iy), 2);
    }
  }
  return cost;
}

// // [[Rcpp::export]]
// NumericMatrix cost_matrix_(const NumericMatrix & x, const NumericMatrix & y){
//   NumericMatrix cost(x.rows(), y.rows());
//   for (int ix = 0; ix < x.rows(); ix ++){
//     for (int iy = 0; iy < y.rows(); iy ++){
//       cost(ix,iy) = 0;
//       // loop over components, stored in columns
//       for (int id = 0; id < x.cols(); id ++){
//         // cost(ix,iy) += std::abs(x(ix,id)-y(iy,id));
//         cost(ix,iy) += std::pow(x(ix,id)-y(iy,id),2);
//       }
//       cost(ix,iy) = std::sqrt(cost(ix,iy));
//     }
//   }
//   return cost;
// }

// Takes two matrices of size d times N1 and d times N2, returns an N1 times N2 matrix
// [[Rcpp::export]]
NumericMatrix cost_matrix_L1_(const NumericMatrix & x, const NumericMatrix & y){
  NumericMatrix cost(x.cols(), y.cols());
  for (int ix = 0; ix < x.cols(); ix ++){
    for (int iy = 0; iy < y.cols(); iy ++){
      cost(ix,iy) = 0;
      // loop over components, stored in columns
      for (int id = 0; id < x.rows(); id ++){
        // cost(ix,iy) += std::abs(x(ix,id)-y(iy,id));
        cost(ix,iy) += std::abs(x(id,ix)-y(id,iy));
      }
    }
  }
  return cost;
}
// Takes two matrices of size d times N1 and d times N2, returns an N1 times N2 matrix
// [[Rcpp::export]]
NumericMatrix cost_matrix_L2_(const NumericMatrix & x, const NumericMatrix & y){
  NumericMatrix cost(x.cols(), y.cols());
  for (int ix = 0; ix < x.cols(); ix ++){
    for (int iy = 0; iy < y.cols(); iy ++){
      cost(ix,iy) = 0;
      // loop over components, stored in columns
      for (int id = 0; id < x.rows(); id ++){
        // cost(ix,iy) += std::abs(x(ix,id)-y(iy,id));
        cost(ix,iy) += std::pow(x(id,ix)-y(id,iy),2);
      }
      cost(ix,iy) = std::sqrt(cost(ix,iy));
    }
  }
  return cost;
}

