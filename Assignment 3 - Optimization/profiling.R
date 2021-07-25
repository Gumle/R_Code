library(splines) # splineDesign (B-splines)
library(LaplacesDemon) # logit and invlogit
library(Rcpp) # c++ 
library(profvis) # profiling
library(numDeriv) # numerical methods for grad and hess
library(ggplot2) 
library(ggpubr)
library(microbenchmark) # benchmarking
library(bench) 
library(gridExtra)
library(dplyr) # data manipulation
library(Matrix) # sparse matrix
setwd("~/Documents/Statistics - KU/Computational Statistics/Assignments/Assignment 3 - Optimization")
source("Debugging_and_tracing.R")
source("Assignment-3-functions.R")
library(simts)
n = 10000

# creating knots for B-spline basis
knots <- function(inner_knots) {
  sort(c(rep(range(inner_knots), 3), inner_knots))
}
# f(x|beta) = (phi1(x), phi2(x),...,phim(x)) * beta
f <- function(par, x ,inner_knots) {
  
  if(length(par) != length(inner_knots) + 2) {
    stop("none compaterble dimensions")
  }
  phi <- splineDesign(knots(inner_knots), x) # designmatrix
  phi %*% par
}
xx <- seq(0, 1000, length.out = n)
inner_knotsxx <- seq(range(xx)[1], range(xx)[2], length.out = 3)
par0 <- rnorm(5)
pvaerdier <- function(x)
{
  f <- f(par0, x, inner_knotsxx)  #0.1 + 1.2*(x-0.5)^2 + 0.9*x^3  + 0.3*x^4 + 0.2*x^5
  
  exp(f)/(1 + exp(f))
}
yy <- rbinom(n, 1, pvaerdier(xx))
xx <- sample(xx)

profvis(Newton(par0, H, grad_H, hessian_H, maxiter = 500, d = 0.1, c = 0.1,  gamma0 = 1, epsilon = 1e-5, stop = 'func', cb = NULL, xx, yy, lambda, inner_knotsxx))
profvis(GD(par0, H, grad_H, hessian_H, maxiter = 10000, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-5, stop = 'func', cb = NULL, xx, yy, lambda, inner_knotsxx))
profvis(CG(par0, H, grad_H, hessian_H, maxiter = 10000, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-5, stop = 'func', cb = NULL, xx, yy, lambda, inner_knotsxx))

