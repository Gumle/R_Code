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
# Test setup values:
df <- read.csv('Horses.csv', header = T)
df <- na.omit(df) 
x <- df$Temperature
y <- df$dead

Newton <- function(par, H, grad_H, hessian_H, maxiter = 50, d = 0.1, c = 0.1,  gamma0 = 1, epsilon = 1e-5, stop = 'grad', cb = NULL, ...){
  browser()
  n <- 1
  repeat {
    # Status of objective function and derivatives
    H_val <- H(par, ...)
    grad <- grad_H(par, ...)
    grad_norm <- sum(grad^2)
    hessian <- hessian_H(par, ...)
    # Descent 
    rho <- -drop(solve(hessian, grad)) # Descent direction
    gamma <- gamma0 # Step size
    par1 <- par + gamma * rho   # Proposed descent step
    h_prime <- t(grad) %*% rho # inner prodict of gradient and decent direction. 
    if(!is.null(cb)) cb() # Tracing 
    while(H(par1, ...) > H_val + c * gamma * h_prime) { # Backtracking
      gamma <- d * gamma
      par1 <- par + gamma * rho 
    }
    H_val1 <- H(par1, ...)
    
    # Condition choices - depends on default value of stop.
    if (n > maxiter) { 
      stopping <- 'maxiter' 
      break 
    }
    else if (stop == 'func') {
      # Small relative descent for objective function
      if (H_val - H_val1 < epsilon * (H_val1 + epsilon)) {
        stopping <- 'converged' 
        break 
      }
    } 
    else if (stop == 'grad') {
      # Small gradient 
      if (grad_norm <= epsilon) {
        stopping <- 'converged' 
        break
      }
    }
    else if (stop == 'par') {
      # Small relative change in parameters
      if (sum((par1 - par)^2) <= epsilon * ((sum(par1)^2) + epsilon)) {
        stopping <- 'converged' 
        break
      } 
    }
    par <- par1; n <- n + 1 
  } 
  
  list(par = par, H_val = H_val, grad_norm = grad_norm, stop = stopping, iter = n) 
}
Newton(par0, H, grad_H, hessian_H, maxiter = 500, d = 0.1, c = 0.1,gamma0 = 1, epsilon = 1e-5, stop = 'grad', cb = NULL, x, y, lambda, inner_knots)
