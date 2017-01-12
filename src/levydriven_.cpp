#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
using namespace Eigen;


// [[Rcpp::export]]
List levydriven_rtransition_rand_cpp(int nparticles, NumericVector & theta){
  RNGScope scope;
  NumericVector sum_weighted_e(nparticles);
  NumericVector sum_e(nparticles);
  NumericVector k = rpois(nparticles, theta[4] * theta[2] * theta[2] / theta[3]);
  for (int i = 0; i < nparticles; i++){
    NumericVector c = runif(k(i));
    NumericVector e = rexp(k(i), theta[2]/theta[3]); // rcpp rexp is parametrized with rate
    sum_e(i) = sum(e);
    sum_weighted_e(i) = sum(exp(-theta[4] * c) * e);
  }
  return(List::create(Named("sum_weighted_e") = sum_weighted_e, Named("sum_e") = sum_e));
}
