#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double mmd_c(double first_term, double eps, const NumericMatrix & x, const NumericMatrix & y){
  int nobs = y.cols();
  int dimension = y.rows();
  double result = first_term;
  double cost_xx;
  double cost_xy;
  double second_term = 0;
  double third_term = 0;
  for (int i1 = 0; i1 < nobs; i1 ++){
    for (int i2 = 0; i2 < nobs; i2 ++){
      cost_xx = 0;
      cost_xy = 0;
      for (int id = 0; id < dimension; id ++){
        cost_xx += std::pow(x(id,i1)-x(id,i2),2);
        cost_xy += std::pow(x(id,i1)-y(id,i2),2);
      }
      second_term += exp(-cost_xx / (2*eps*eps));
      third_term += exp(-cost_xy / (2*eps*eps));
    }
  }
  result += second_term / (nobs * nobs) - 2 * third_term / (nobs*nobs);
  return result;
}
