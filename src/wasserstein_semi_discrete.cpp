#include <Rcpp.h>
using namespace Rcpp;


// Wasserstein distance between one continuous and one discrete measure. See Genevay et al (2016) for details.

//The cost function is taken to be Euclidean distance to some power p \geq 1. Here we define this distance.
// [[Rcpp::export]]
double eucl_distn_cpp(NumericVector x, NumericVector y) {
  int n = x.size();
  double out = 0;
  
  for(int i = 0 ; i < n; ++i){
    out += pow(x[i]-y[i],2);
  }
  out = sqrt(out);
  return out;
}


/*Computes \bar{h}_\varepsilon(x, \bm{v}), which also depends on the given weights \bm{\nu} and corresponding
observations y of the discrete measure. Only takes epsilon > 0. p corresponds to the p-Wasserstein distance */
// [[Rcpp::export]]
double h_cpp(NumericVector x, NumericVector v, NumericVector nu, NumericMatrix y, double epsilon, double p) {
  int n = nu.size();
  double dummy1 = 0;  
  NumericVector dummy2(n); 
  double dummy3 = 0;  
  double dummy4 = 0;
  double out = 0;
  
  for(int i = 0; i<n; ++i){
    dummy1 += nu[i]*v[i];   //Computes first sum in expression of \bar{h}_\varepsilon(x, \bm{v})
    dummy2[i] = log(nu[i]) + (v[i]-pow(eucl_distn_cpp(x,y.row(i)),p))/epsilon;   //Computes log of terms in second sum.
  }
  
  dummy3 = *std::max_element(dummy2.begin(),dummy2.end());
  for(int i = 0; i<n; ++i){
    dummy4 += exp(dummy2[i] - dummy3);  //This normalization helped with machine precision problems when computing the sum.
  }
  
  out = dummy1 - epsilon*(log(dummy4)+dummy3)-epsilon;  //Equals \bar{h}_\varepsilon(x, \bm{v}), notice renormalization.
  return out;
}


/*Computes the gradient of \bar{h}_\varepsilon(x, \bm{v}) wrt \bm{\nu}, which also depends on the given weights \bm{\nu} and
corresponding observations y of the discrete measure. Only takes epsilon > 0. p corresponds to the p-Wasserstein distance*/
// [[Rcpp::export]]
NumericVector grad_h_cpp(NumericVector x, NumericVector v, NumericVector nu, NumericMatrix y, double epsilon, double p) {
  int n = nu.size();
  NumericVector sum_comp(n);
  NumericVector log_sum_comp(n);
  double sum = 0;
  double dummy = 0;
  NumericVector out(n);
  
  
  for(int i = 0; i<n; ++i){
    log_sum_comp[i] = log(nu[i]) + (v[i]-pow(eucl_distn_cpp(x,y.row(i)),p))/epsilon;  //log of terms in sum in \chi
  }
  dummy = *std::max_element(log_sum_comp.begin(),log_sum_comp.end());
  for(int i = 0; i<n; ++i){
    sum_comp[i] = exp(log_sum_comp[i] - dummy);  //normalization to resolve numerical issues when computing the sum.
    sum += sum_comp[i];  //normalized version of sum in \chi
  }
  
  for(int i = 0; i<n; ++i){
    out[i] = nu[i] - sum_comp[i]/sum;  //Components of the gradient. Notice that normalization cancels.
  }
  return out;
}


/*Performes Averaged SGD to compute the v attaining max in semi-dual formaulation, stops when threshold is reached */
// [[Rcpp::export]]
NumericVector dw_sgd_v_thresh_cpp(NumericMatrix samples, double c, NumericVector nu, NumericMatrix y,
                                  double epsilon, double thresh, double p){
  int n = samples.nrow();
  int m = nu.size();
  NumericVector v(m);
  NumericVector w(m);
  NumericVector dummy1(m);
  NumericVector dummy2(m);
  double dummy_test = 1;
  NumericVector out(m);
  
  for(int i = 0; i < n; ++i){
    dummy1 = clone(v);
    w = w + c/sqrt(i+1)*grad_h_cpp(samples.row(i),w,nu,y,epsilon,p);
    v = 1.0/((double)(i+1))*w + ((double) i)/((double)(i+1))*v;
    dummy_test = eucl_distn_cpp(v,dummy1)/eucl_distn_cpp(v,dummy2);  
    //Computes the distance between current v and its previous value, divided by the length of v. Use in tolerance check.
    
    if((dummy_test <  thresh) && (i > 100)){ //If we've done at least 100 iterations and dummy_test below threshold, stop. 
      break;
    }
  }
  out = v;
  return out;
}

/*Performes Averaged SGD to compute the v attaining max in semi-dual formaulation, stops after n iterations 
are performed. n is the number of samples simulted from the model. */
// [[Rcpp::export]]
NumericVector dw_sgd_v_cpp(NumericMatrix samples, double c, NumericVector nu, NumericMatrix y, double epsilon, double p){
  int n = samples.nrow();
  int m = nu.size();
  NumericVector v(m);
  NumericVector w(m);
  double dummy_test = 1;
  NumericVector out(m);
  
  for(int i = 0; i < n; ++i){
    w = w + c/sqrt(i+1)*grad_h_cpp(samples.row(i),w,nu,y,epsilon,p);
    v = 1.0/((double)(i+1))*w + ((double) i)/((double)(i+1))*v;
  }
  out = v;
  return out;
}



/*Uses a naive estimate of the Wasserstein distance by computing the sample mean of \bar{h}_\varepsilon(x, \bm{v}) 
at the \bm{v} output from dw_sgd_v_thres_cpp or dw_sgd_v_cpp, with samples x_j drawn from the model. nu, y, epsilon
and y here need to be the same as in the computation of v.*/
// [[Rcpp::export]]
double dw_est_cpp(NumericMatrix samples_for_est, NumericVector v, NumericVector nu,
                  NumericMatrix y, double epsilon, double p){
  
  int k = samples_for_est.nrow();
  double out = 0;
  
  for(int i = 0; i < k; ++i){
    out += h_cpp(samples_for_est.row(i),v,nu,y,epsilon,p)/k;
    out = pow(out,1/p);
  }
  return out;
}
