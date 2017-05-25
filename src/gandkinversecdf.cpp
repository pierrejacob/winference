#include <math.h>
#include <iostream>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
using namespace Eigen;

// [[Rcpp::export]]
NumericVector gandkinversecdf_(NumericVector x, NumericVector theta){
  NumericVector z = qnorm(x);
  double A = theta(0);
  double B = theta(1);
  double c = 0.8;
  double g = theta(2);
  double k = theta(3);
  return A + B * (1 + c * (1 - exp(- g * z)) / (1 + exp(- g * z))) * pow((1 + pow(z, 2.0)), k) * z;
}

// [[Rcpp::export]]
NumericVector gandkinversecdf_givennormals_(NumericVector normals, NumericVector theta){
  NumericVector z = normals;
  double A = theta(0);
  double B = theta(1);
  double c = 0.8;
  double g = theta(2);
  double k = theta(3);
  return A + B * (1 + c * (1 - exp(- g * z)) / (1 + exp(- g * z))) * pow((1 + pow(z, 2.0)), k) * z;
}

// [[Rcpp::export]]
double gandkcdf_(double y, NumericVector theta, int maxsteps = 1000, double tolerance = 1e-10,
                 double lower = 1e-20, double upper = 1-1e-20){
  double A = theta(0);
  double B = theta(1);
  double c = 0.8;
  double g = theta(2);
  double k = theta(3);

  int istep = 0;
  double current_try = (upper + lower) / 2;
  double current_size = (upper - lower) / 4;
  NumericVector dd(1);
  dd(0) = current_try;
  NumericVector z = qnorm(dd);
  double fattempt = A + B * (1 + c * (1 - exp(- g * z(0))) / (1 + exp(- g * z(0)))) * pow((1 + pow(z(0), 2.0)), k) * z(0);
  while (!(fattempt > y-tolerance && fattempt < y+tolerance) && (istep < maxsteps)){
    istep++;
    if (fattempt > y-tolerance){
      current_try = current_try - current_size;
      dd(0) = current_try;
      NumericVector z = qnorm(dd);
      fattempt = A + B * (1 + c * (1 - exp(- g * z(0))) / (1 + exp(- g * z(0)))) * pow((1 + pow(z(0), 2.0)), k) * z(0);
      current_size = current_size / 2;
    } else {
      current_try = current_try + current_size;
      dd(0) = current_try;
      NumericVector z = qnorm(dd);
      fattempt = A + B * (1 + c * (1 - exp(- g * z(0))) / (1 + exp(- g * z(0)))) * pow((1 + pow(z(0), 2.0)), k) * z(0);
      current_size = current_size / 2;
    }
  }
  return current_try;
}
