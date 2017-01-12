#include <RcppEigen.h>
#include "weighted_averages.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector wmean_(const NumericMatrix & x, const NumericVector & unnormalized_w){
  int nrows = x.rows();
  int ncols = x.cols();
  double sumw = sum(unnormalized_w);
  NumericVector result(ncols);
  double cumsumxw;
  for (int icol = 0; icol < ncols; icol++){
    cumsumxw = 0.;
    for (int irow = 0; irow < nrows ; irow++){
      cumsumxw += unnormalized_w(irow) * x(irow, icol);
    }
    result(icol) = cumsumxw / sumw;
  }
  return result;
}

// [[Rcpp::export]]
NumericMatrix wcovariance_(const NumericMatrix & x, const NumericVector & unnormalized_w, const NumericVector & xbar){
  int nrows = x.rows();
  int ncols = x.cols();
  double sumw = sum(unnormalized_w);
  double sumsqw = sum(unnormalized_w*unnormalized_w);
  NumericMatrix result(ncols, ncols);
  std::fill(result.begin(), result.end(), 0);
  for (int i = 0; i < ncols; i++){
    for (int j = 0; j < ncols; j++){
      for (int irow = 0; irow < nrows ; irow++){
        result(i, j) += unnormalized_w(irow) * (x(irow, i) - xbar(i)) * (x(irow, j) - xbar(j));
      }
      result(i,j) /=  (sumw - sumsqw / sumw);
    }
  }
  return result ;
}

RCPP_MODULE(module_waverage) {
  function( "wmean", &wmean_ );
  function( "wcovariance", &wcovariance_ );
}



