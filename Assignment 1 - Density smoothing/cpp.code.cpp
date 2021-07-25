// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

// Function computes weigths
// [[Rcpp::export]]
NumericVector kbin_c(NumericVector x, double lo, double hi, int m) {
  int n = x.size();
  NumericVector w(m);
  double delta = (hi - lo) / (m - 1);
  for (int i = 0; i < n; i++)
  {
    int ii = floor((x[i] - lo) / delta + 0.5);
    w[ii] = w[ii] + 1;
  }
  
  return w*1/n;
}

// Function constructs a grid of m equidistant points within a specified range
// [[Rcpp::export]]
NumericVector grid(double from, double to, int m){
  NumericVector grid = wrap(seq(0,m-1));
  return (grid/grid(m-1))*(to-from)+(from);
}


// Function evaluates the epanechnikov kernal in (x - xi)/h
// [[Rcpp::export]]
NumericVector kernC(NumericVector z) {
  int m = z.size();
  NumericVector out(m);
  for(int i = 0; i < m; ++i){
    if (std::abs(z[i]) <= 1) {
      out[i] = 0.75*(1-pow(z[i],2.0));
    }
    else {
      out[i] = 0;
    }
  }
  return out;
}

// Function evaluates kernal density estimate with the epanechinkov kernal in grid
// [[Rcpp::export]]
vec kernlC(NumericVector grid,  double h){
  NumericVector arg = (grid - grid[0])/h;
  vec y = vec(kernC(arg)/h);
  return y;
}

// Function that wraps evaluation points and values into a list
// [[Rcpp::export]]
List kdens_binc(NumericVector x, double h, int m){
  double minn = min(x)-3*h;
  double maxx = max(x)+3*h;
  NumericVector xx = grid(minn,maxx,m);
  vec weights = vec(kbin_c(x, xx[0], xx[xx.size() - 1], m));
  vec kerneval = kernlC(xx,h);
  mat kerndif = toeplitz(kerneval);
  mat y = kerndif.t()*weights;
  List fina(2);
  fina[0] = xx;
  fina[1] = y;
  return fina;
}


// [[Rcpp::export]]
NumericVector test(NumericVector x, double h, int m) {
  double minn = min(x)-3*h;
  cout << minn <<endl;
  double maxx = max(x)+3*h;
  cout << maxx << endl;
  NumericVector xx = grid(minn,maxx,m);
  NumericVector weights = kbin_c(x, xx[0], xx[(xx.size() - 1)], m);
  return weights;
}