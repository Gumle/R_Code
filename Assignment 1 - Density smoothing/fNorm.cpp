#include <math.h>
using namespace Rcpp;
using namespace arma;
using namespace std;
//m [[Rcpp:export]]
double fNorm(NumericVector x, double h){
  int n = x.size();
  double par;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      par = par + exp(-0.5(x[i]-x[j])*(x[i]-x[j]));
    }
  }
  double sum = par/((h*h)*(n*n)*2*sqrt(M_PI*h*h));
}


//m [[Rcpp:export]]
double fNorm(NumericVector x, double h){
  int n = x.size();
  double par = 0.0;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      par = par + exp(-0.5*(x[i]-x[j])*(x[i]-x[j])/(2*h*h));
    }
  }
  double div = (n*n)*2*sqrt(M_PI*h*h);
  double sum = par/div;
  return sum;
}

//m [[Rcpp:export]]
double fNorm(NumericVector x, double h){
  int n = x.size();
  double par = 0.0;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      par = par + exp(-0.5*pow(x[i]-x[j],2)/(2*pow(h,2)));
    }
  }
  double div = pow(n,2)*2*sqrt(M_PI*pow(h,2));
  double sum = par/div;
  return sum;
}
