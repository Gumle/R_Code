library(ggplot2)
library(microbenchmark)
library(Rcpp)
library(RcppArmadillo)
setwd("~/Documents/Statistics - KU/Computational Statistics/Assignments/Assignment 1 - Density smoothing")

# Loop-based
kernDens_loop <- function(x, h, m = 512) {
  
  rg <- range(x) + c(-3 * h, 3 * h)
  xx <- seq(rg[1], rg[2], length.out = m) # evaluation grid around x : m points
  y <- numeric(m) # vector to store the function values - same length as xx
  const <- 3/(4* h*length(x))
  for (i in seq_along(xx)){
    condition <- xx[i]- x
    kde <- ifelse(abs(condition) <= h, 1 - (condition/h)^2, 0) 
    y[i] <- sum(kde)*const}
  list(x = xx, y = y) # return evauation grid and the function values
}

# Sapply

kernDens_apply <- function(x, h, m = 512){
  rg <- range(x) + c(-3 * h, 3 * h)
  xx <- seq(rg[1], rg[2], length.out = m) # evaluation grid around x : m points 
  y <- numeric(m) # vector to store the function values - same length as xx
  const <- (3/4) / (h * length(x))
  y <- sapply(xx, function(z) sum((1 - ((z - x) / h)^2 ) * (abs(x-z) <= h)))
  y <- const * y
  list(x = xx, y = y) # return evauation grid and the function values
}

# binning

bin <- function (x, lo, hi, m) {
  w <- numeric(m)
  delta <- (hi - lo)/(m - 1) # bin length
  for (i in seq_along(x)) {
    ii <- floor((x[i] - lo)/delta + 0.5) + 1 # determines which bin xi is in
    w[ii] <- w[ii] + 1 # counts the number of observations in each bin
  }
  w <- w/sum(w) #the weights are computed as the number of observatios in each bin devided by n
}

kernDens_bin <- function (x, h, m = 512) {
  rg <- range(x) + c(-3 * h, 3 * h)
  xx <- seq(rg[1], rg[2], length.out = m) # evaluation grid around x : m points 
  weights <- bin(x, rg[1], rg[2], m) # compute weights 
  kerneval <- (3/4) * (1 - ((xx - xx[1])/h) ^ 2) * (abs(xx - xx[1]) <= h) * 1/h
  kerndif <- toeplitz(kerneval)
  y <- colSums(weights * kerndif)
  list(x = xx, y = y) # return evaluation grid and function values 
}

# Maybe make this more generalaized:

# Cross-validation

cv <- function (x, h){
  n <- length(x)
  y <- 0
  for (i in 1:n){
    y[i] <- sum(3/4 * (1 - ((x[i] - x[-i])/h)^2) * (abs(x[i] - x[-i]) <= h)) * 1/(h * (n-1))
  }
  value <- sum(log(y[y!=0])) 
  value
}


cv.factory <- function(x) {
  f <- function(h) cv(x, h)
  #Vectorize(f)
}

x <- rnorm(1000)
x <- F12
f <- cv.factory(x)
h.cv <- optimize(f, interval = c(0.001, 2), maximum = TRUE)$maximum

par(mfrow = c(1,1))
hist(x, probability = T, col = 'dodgerblue4', main = 'Cross validation h',  ylim=c(0,0.6))
lines(kernDens_bin(x,h.cv), lwd = 2, col = 'chocolate1')


############################ UCV - PLOT THE HISTOGAM ##########################3
# Testing

setwd("~/Documents/Statistics - KU/Computational Statistics/Assignments/Assignment 1 - Density smoothing")
dat <- read.table("infrared.txt", header = TRUE, sep = "")
x <- rnorm(10)
h <-0.2
F12 <- log(dat$F12)
x <- F12

# UCV

Kbar <- function(x) {
  #browser()
  y <- numeric(length(x))
  idx <- which(abs(x) <= 2)
  l <- sapply(x, function(z) max(-1, z - 1))
  u <- sapply(x, function(z) min(1, z + 1))
  y[idx] <- (3 / 16) * x[idx] ^ 2 * u[idx] ^ 3 - (9 / 16) *  x[idx] ^ 2 * u[idx] - 
    (3 / 16) * x[idx] ^ 2 * l[idx] ^ 3 + (9 / 16) * x[idx] ^ 2 * l[idx] - 
    (9 / 32) * x[idx] * u[idx] ^ 4 + (9 / 16) * x[idx] * u[idx] ^ 2 +
    (9 / 32) * x[idx] * l[idx] ^ 4 - (9 / 16) * x[idx] * l[idx] ^ 2 + 
    (9 / 80) * u[idx] ^ 5 - (3 / 8) * u[idx] ^ 3 + (9 / 16) * u[idx] - 
    (9 / 80) * l[idx] ^ 5 + (3 / 8) * l[idx] ^ 3 - (9 / 16) * l[idx]
  y
}

x <- seq(-3,3,0.1)
Kbar(x)
cppFunction(
"
#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]  /* to show that this function is to be exported to R */  
double fNorm(NumericVector x, double h, Function Kbar){
 int n = x.size();
 double par = 0.0;
 double constant = 1/(pow(n,2)*h);
 for (int i = 0; i < n; i++){
   for (int j = 0; j < n; j++){
     par = par + Kbar((x[i]-x[j])/h);
   }  
 }
 double sum = par/constant;
 return sum;
}")

cppFunction("
using namespace std;
double Kbar(NumericVector x){
  int n = x.size();
  NumericVector y (n);
  for(int i = 0; i < n; i++){
    if (abs(x[i]) <= 2) {
      y[i] = (3.0 / 16.0) * pow(x[i],2) * pow(min(1,x[i]+1),3) - (9.0 / 16.0) *  pow(x[i], 2) * min(1,x[i]+1) - 
        (3.0 / 16.0) * pow(x[i],2) * pow(max(-1,x[i]-1),3) + (9.0 / 16.0) * pow(x[i],2) * max(-1,x[i]-1) - 
        (9.0 / 32.0) * x[i] * pow(min(1,x[i]+1),4) + (9.0 / 16.0) * x[i] * pow(min(1,x[i]+1), 2) +
        (9.0 / 32.0) * x[i] * pow(max(-1,x[i]-1), 4) - (9.0 / 16.0) * x[i] * pow(max(-1,x[i]-1), 2) + 
        (9.0 / 80.0) * pow(min(1,x[i]+1), 5) - (3.0 / 8.0) * pow(min(1,x[i]+1), 3) + (9.0 / 16.0) * min(1,x[i]+1) - 
        (9.0 / 80.0) * pow(max(-1,x[i]-1), 5) + (3.0 / 8.0) * pow(max(-1,x[i]-1), 3) - (9.0 / 16.0) * max(-1,x[i]-1) ;
      }
  }
  return y;
}")
K_bar.tobias <- function(x){
  F_lowerbound <- 3 * (6*(max(-1,x-1))^5 - 15*x*(max(-1,x-1))^4 + (max(-1,x-1))^2*(10*x^2 - 20) + 30*x*(max(-1,x-1)) -30*x^2 + 30) / 160
  F_upperbound <- 3 * (6*(min(1,x+1))^5 - 15*x*(min(1,x+1))^4 + (min(1,x+1))^2*(10*x^2 - 20) + 30*x*(min(1,x+1)) -30*x^2 + 30) / 160
  F_upperbound- F_lowerbound
}

Kbar <- function(x) {
  y <- numeric(length(x))
  idx <- which(abs(x) <= 2)
  l <- sapply(x, function(z) max(-1, z - 1))
  u <- sapply(x, function(z) min(1, z + 1))
  y[idx] <- (3 / 16) * x[idx] ^ 2 * u[idx] ^ 3 - (9 / 16) *  x[idx] ^ 2 * u[idx] - 
    (3 / 16) * x[idx] ^ 2 * l[idx] ^ 3 + (9 / 16) * x[idx] ^ 2 * l[idx] - 
    (9 / 32) * x[idx] * u[idx] ^ 4 + (9 / 16) * x[idx] * u[idx] ^ 2 +
    (9 / 32) * x[idx] * l[idx] ^ 4 - (9 / 16) * x[idx] * l[idx] ^ 2 + 
    (9 / 80) * u[idx] ^ 5 - (3 / 8) * u[idx] ^ 3 + (9 / 16) * u[idx] - 
    (9 / 80) * l[idx] ^ 5 + (3 / 8) * l[idx] ^ 3 - (9 / 16) * l[idx]
  y
}

f2norm <- function(x,h){
  n <- length(x)
  constant <- 1/(n^2*h)
  y <- 0
  for(i in 1:n){
    for(j in 1:n){
    y <- y +  Kbar((x[i]-x[j])/h)
    }
  }
  constant*y
}

f2norm(x = F12, h = 0.5)
f2norm.factory <- function(x){
  function(h) f2norm(x,h)
}
f2norm.data <-f2norm.factory(F12)

fSum <- function (x, h){
  n <- length(x)
  y <- 0
  for (i in 1:n){
    y[i] <- sum(3/4 * (1 - ((x[i] - x[-i])/h)^2) * (abs(x[i] - x[-i]) <= h)) * 1/(h * (n-1))
  }
  value <- 2*sum(y)/n
  value
}
# Generalized
EKernel <- function(x){ ifelse(-1 <= x & x<= 1,x <- (3/4)*(1-x^2), 0) }
fSumGeneral <- function (x, h, kernel){
  n <- length(x)
  y <- 0
  for (i in 1:n){
    y[i] <- sum(kernel((x[i]/h)-kernel(x[-i]/h )) ) * 1/(h * (n-1))
  }
  value <- 2*sum(y)/n
  value
}

fSum(F12, 0.5)
EKernel()
fSumGeneral(F12, 0.5, EKernel(F12))
fSum.factory <- function(x){
  function(h) fSum(x,h)  
}
fSum.factory()
fSum.data <- fSum.factory(F12)

UCV <- function(h){
  f2norm.data(h) - fSum.data(h)
}

h.opt <- optimize(UCV, interval = c(0.0001,2))$minimum
h.opt


hist(x, probability = T, col = 'dodgerblue4', main = 'UCV h',  ylim=c(0,0.6))
lines(kernDens_bin(x,h.opt), lwd = 2, col = 'chocolate1')


######################### Testing our UCV ##############
integral(x = F12,h = 0.2)
EKernel <- function(x){
  ifelse(-1 <= x & x<= 1,x <- (3/4)*(1-x^2), 0) 
}

integral2 <- function(x){
  f <- integrate(function(u) EKernel(u)*EKernel(x-u), lower = x-1, upper = 1)$value
  f
}
x = 1
K_bar.tobias(1.2)
Kbar3(1.2)
integral2(1.2)


kern_epa_conv(3)
integral2(3)
