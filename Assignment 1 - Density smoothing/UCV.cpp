#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
using namespace std;
using namespace Rcpp;

NumericVector Kbar(NumericVector x){
  int n = x.size();
  NumericVector y(n);
  for(int i = 0; i < n; i++){
    if(abs(x[i]) <= 2) {
      y[i] = (3.0 / 16.0) * pow(x[i], 2)*pow(min(1,x[i]+1),3) - (9.0 / 16.0) *  pow(x[i], 2) * min(1,x[i]+1) - 
        (3.0 / 16.0) * pow(x[i],2) * pow(max(-1,x[i]-1),3) + (9.0 / 16.0) * pow(x[i],2) * max(-1,x[i]-1) - 
        (9.0 / 32.0) * x[i] * pow(min(1,x[i]+1),4) + (9.0 / 16.0) * x[i] * pow(min(1,x[i]+1), 2) +
        (9.0 / 32.0) * x[i] * pow(max(-1,x[i]-1), 4) - (9.0 / 16.0) * x[i] * pow(max(-1,x[i]-1), 2) + 
        (9.0 / 80.0) * pow(min(1,x[i]+1), 5) - (3.0 / 8.0) * pow(min(1,x[i]+1), 3) + (9.0 / 16.0) * min(1,x[i]+1) - 
        (9.0 / 80.0) * pow(max(-1,x[i]-1), 5) + (3.0 / 8.0) * pow(max(-1,x[i]-1), 3) - (9.0 / 16.0) * max(-1,x[i]-1) ;
      }
  }
  return y;
}

//m [[Rcpp:export]]
double fNorm(NumericVector x, double h){
  int n = x.size();
  double par = 0.0;
  dobule constant = 1/(pow(n,2)*h);
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      par = par + Kbar((x[i]-x[j])/h);
    }  
  }
  double sum = par/constant;
  return sum;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//