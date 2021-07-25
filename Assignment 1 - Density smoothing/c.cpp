// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

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

//vec A = randu<vec>(5);
//mat X = toeplitz(A);

//double toeplitz_matrix(int x) {
  //toeplitz(randu<vec>(x))
  //vec A = randu<vec>(x);
  //return A[0];
//}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
library(ggplot2)
library(microbenchmark)
#x = rnorm(10000)
#x = seq(0,12,1)
#h = 0.8
#lo = x[1] -3*h
#hi = x[length(x)] + 3*h
##(NumericVector x, double lo, double hi, int m

#sum(kbin_c(x, lo, hi, 512) - kbin(x, lo, hi, 512))
#test(x, h, 512)
#lo = x[1] + -3*0.5
#hi = x[length(x)] + 3*0.5
#kbin_c(x, lo, hi, 512)
#grid(x[1] -3*0.5, x[length(x)] + 3*0.5, 512)


#l = kdens_binc(x,h, 512)
#l[[2]] = as.numeric(l[[2]])
#u = list(x = l[[1]], y = l[[2]])
#hist(x, breaks = 15, probability = T, col = '#2F4F4F', main = 'Epanechnikov kernal c++')
#lines(u, lwd = 2, col = 'red')

x = rnorm(10000)

  kern_bench = microbenchmark( times = 100,
                               kdens_bin(x, 0.2),
                               kdens_binc(x, 0.2, 512)
  )

summary(kern_bench)

autoplot(kern_bench) + geom_jitter(aes(color = expr), alpha = 0.4) + aes(fill = I("gray"))+
    scale_color_manual(values = c("#c91246", "#9400D3")) + theme_bw() +
    theme(legend.position = "none")

*/
