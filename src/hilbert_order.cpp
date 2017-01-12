#include <RcppEigen.h>
#include "HilbertCode.h"
using namespace Rcpp;
using namespace std;
using namespace Eigen;

// [[Rcpp::export]]
NumericVector hilbert_order_(NumericMatrix x){
  int dx = x.nrow();
  int N = x.ncol();
  NumericVector order_indices(N);
  Hilbert_Sort_CGAL(REAL(x), dx, N, REAL(order_indices));
  return order_indices;
}
