#include <RcppEigen.h>
#include "resampling.h"
using namespace Rcpp;
using namespace std;


void permute(IntegerVector & a, const int & nparticles){
  int swap;
  for (int i = 0; i < nparticles; i++){
    if (a(i) != i && a(a(i)) != a(i)){
      swap = a(a(i));
      a(a(i)) = a(i);
      a(i) = swap;
      i--;
    }
  }
}

// weights have to sum to 1
void systematic(IntegerVector & ancestors, const NumericVector & weights, const int & nparticles){
  RNGScope scope;
  NumericVector u_vec = runif(1);
  double u = u_vec(0) / nparticles;
  int j = 0;
  double csw = weights(0);
  for(int k = 0; k < nparticles; k++){
    while(csw < u){
      j++;
      csw += weights(j);
    }
    ancestors(k) = j;
    u = u + 1./nparticles;
  }
}

// made for export to R
// [[Rcpp::export]]
IntegerVector systematic_resampling_(const NumericVector & weights){
  RNGScope scope;
  int nparticles = weights.size();
  IntegerVector ancestors(nparticles);
  NumericVector u_vec = runif(1);
  double u = u_vec(0) / nparticles;
  int j = 0;
  double csw = weights(0);
  for(int k = 0; k < nparticles; k++){
    while(csw < u){
      j++;
      csw += weights(j);
    }
    ancestors(k) = j;
    u = u + 1./nparticles;
  }
  return ancestors + 1;
}

inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// the weights in argument don't need to be normalised, and 'ancestors' don't need to be of the same size as 'weights'
void multinomial(IntegerVector & ancestors, const NumericVector & weights){
    RNGScope scope;
    int nparents = weights.size();
    NumericVector cumsumw = cumsum(weights);
    int nchildren = ancestors.size();
    NumericVector uniforms = runif(nchildren);
    double sumw = cumsumw(nparents - 1);
    double lnMax = 0;
    int j = nparents;
    for (int i = nchildren; i > 0; i--){
      lnMax += log(uniforms(i-1)) / i;
      uniforms(i-1) = sumw * exp(lnMax);
      while (uniforms(i-1) < cumsumw(j-1)){
        j --;
      }
      ancestors(i-1) = j;
    }
    std::random_shuffle(ancestors.begin(), ancestors.end(), randWrapper);
}
