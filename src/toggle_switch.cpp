#include <Rcpp.h>
using namespace Rcpp;


double normal_pdf(double x, double m, double s)
{
  static const double inv_sqrt_2pi = 0.3989422804014327;
  double a = (x - m) / s;
  
  return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}


// [[Rcpp::export]]
NumericVector cont_hist(NumericVector x, double h, NumericVector prob, NumericVector mids) {
  
  int n = x.size();
  int m = prob.size();
  
  double sum;
  NumericVector out(n);
  
  for(int i = 0; i<n; i++){
    sum = 0.0;
    for(int j = 0; j<m; j++){
      sum += prob[j]*normal_pdf(x[i],mids[j],h);
    }
    out[i] = sum;
  }
  
  return out;
}

// [[Rcpp::export]]
NumericVector summary_full(NumericVector y, double h, NumericMatrix prob_mat, NumericVector mids) {
  
  int n = prob_mat.nrow();
  int m = y.size();
  NumericVector out(n);
  NumericVector store(m);
  double sum;
  
  for(int i = 0; i<n; i++){
    //NumericMatrix::Row prob = prob_mat(i,_);
    store = cont_hist(y,h,prob_mat(i,_),mids);
    sum = 0.0;
    for(int j = 0; j<m; j++){
      sum += std::log(store[j])/m;
    }
    out[i] = sum;
  }
  return out;
}


