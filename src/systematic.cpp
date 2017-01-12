#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
IntegerVector systematic_resampling_n_(const NumericVector & weights, int ndraws, double u){
  RNGScope scope;
  int nparticles = weights.size();
  IntegerVector ancestors(ndraws);
  u = u / ndraws;
  int j = 0;
  double csw = weights(0);
  for(int k = 0; k < ndraws; k++){
    while(csw < u){
      j++;
      csw += weights(j);
    }
    u = u + 1. / ndraws;
    ancestors(k) = j;
  }
  // int swap;
  // for (int i = 0; i < nparticles; i++){
  //   if (ancestors(i) != i && ancestors(ancestors(i)) != ancestors(i)){
  //     swap = ancestors(ancestors(i));
  //     ancestors(ancestors(i)) = ancestors(i);
  //     ancestors(i) = swap;
  //     i--;
  //   }
  // }
  return ancestors;
}
