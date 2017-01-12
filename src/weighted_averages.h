#ifndef _INCL_WEIGHTED_AVRG_
#define _INCL_WEIGHTED_AVRG_
#include <RcppEigen.h>
using namespace Rcpp;


NumericVector wmean_(const NumericMatrix & x, const NumericVector & unnormalized_w);

NumericMatrix wcovariance_(const NumericMatrix & x, const NumericVector & unnormalized_w, const NumericVector & xbar);

#endif
