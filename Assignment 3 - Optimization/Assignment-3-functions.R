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
# Test setup values:
df <- read.csv('Horses.csv', header = T)
df <- na.omit(df) 
x <- df$Temperature
y <- df$dead


inner_knots <- seq(range(x)[1], range(x)[2], length.out = 3)
par0 <- rnorm(5)
lambda <- 0.4

# Helper functions

knots <- function(inner_knots) {
  sort(c(rep(range(inner_knots), 3), inner_knots))
}

pen_mat <- function(inner_knots) {
  knots <- knots(inner_knots)
  d <- diff(inner_knots)  ## the vector of knot differences; b - a 
  g_ab <- splineDesign(knots, inner_knots, derivs = 2) 
  knots_mid <- inner_knots[-length(inner_knots)] + d / 2
  g_ab_mid <- splineDesign(knots, knots_mid, derivs = 2)
  g_a <- g_ab[-nrow(g_ab), ]
  g_b <- g_ab[-1, ]
  (crossprod(d * g_a,  g_a) + 
      4 * crossprod(d * g_ab_mid, g_ab_mid) + 
      crossprod(d * g_b, g_b)) / 6 
}


# f(x|beta) = (phi1(x), phi2(x),...,phim(x)) * beta
f <- function(par, x ,inner_knots) {
  if(length(par) != length(inner_knots) + 2) {
    stop("none compaterble dimensions")
  }
  phi <- splineDesign(knots(inner_knots), x) # designmatrix
  phi %*% par
} 

H <- function(par, x, y, lambda, inner_knots) {
  f_val <- f(par, x, inner_knots)
  obj_func <- -sum(y * f_val - log(1 + exp(f_val)))
  penalty <- lambda * crossprod(par, pen_mat(inner_knots) %*% par)
  (obj_func + penalty) / length(x)
}
# gradiant of objective

grad_H <- function(par, x, y, lambda, inner_knots) {
  phi <- splineDesign(knots(inner_knots), x)
  f_val <- phi%*%par  
  p <- invlogit(f_val)
  omega <- pen_mat(inner_knots)
  -(crossprod(phi, y - p) - lambda * omega  %*% par)/length(x) 
}

# hessian of objective
hessian_H <- function(par, x, y, lambda, inner_knots) {
  phi <- splineDesign(knots(inner_knots), x)
  f_val <-  phi%*%par
  p <- invlogit(f_val)
  W <- diag(as.vector(p * (1 - p)))
  omega <- pen_mat(inner_knots)
  (crossprod(phi, W %*% phi) + lambda * omega)/length(x)
}



f.sparse <- function(par, x ,inner_knots) {
  if(length(par) != length(inner_knots) + 2) {
    stop("none compaterble dimensions")
  }
  phi <- Matrix(splineDesign(knots(inner_knots), x)) # designmatrix sparse
  phi %*% par
} 


#objective sparse
H.sparse <- function(par, x, y, lambda, inner_knots) {
  f_val <- f(par, x, inner_knots) 
  obj_func <- -sum(y * f_val - log(1 + exp(f_val)))
  omega <-  pen_mat(inner_knots)
  penalty <- lambda * crossprod(par, omega %*% par) # t(par) %*% (omega%*%par) 
  (obj_func + penalty)/length(x)
}


# gradiant of objective sparse
grad_H.sparse <- function(par, x, y, lambda, inner_knots) {
  phi <- Matrix(splineDesign(knots(inner_knots), x))
  f_val <- phi%*%par  
  p <- invlogit(f_val)
  omega <- pen_mat(inner_knots)
  -(crossprod(phi, y - p) - lambda * omega  %*% par)/length(x) 
}

# hessian of objective sparse
hessian_H.sparse <- function(par, x, y, lambda, inner_knots) {
  phi <- Matrix(splineDesign(knots(inner_knots), x))
  f_val <-  phi%*%par
  p <- invlogit(f_val)
  W <- diag(as.vector(p * (1 - p)))
  omega <- pen_mat(inner_knots)
  (crossprod(phi, W %*% phi) + lambda * omega)/length(x)
}

#Newton algorithm

Newton <- function(par, H, grad_H, hessian_H, maxiter = 50, d = 0.1, c = 0.1,  gamma0 = 1, epsilon = 1e-5, stop = 'grad', cb = NULL, ...){
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

# Gradient descent

GD <- function(par, H, grad_H, hessian_H, maxiter = 50, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-5, stop = 'grad', cb = NULL, ...){
  n <--1
  repeat {
    # Status of objective function and derivatives
    H_val <- H(par, ...)
    grad <- grad_H(par, ...)
    grad_norm <- sum(grad^2)
    
    # Descent 
    rho <- drop(-grad) # Descent direction
    #print(rho)
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

#Conjugate gradient

CG <- function(par, H, grad_H, hessian_H, maxiter = 50, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-6, stop = 'grad', cb = NULL, ...){
  n <- 1
  m <- 1
  rho0 <- numeric(length(par))
  
  repeat {
    m <- 1 + m
    
    #Status of objective function and derivatives
    H_val <- H(par, ...)
    grad <- grad_H(par, ...)
    grad_norm <- sum(grad^2)
    
    # Descent 
    rho <- drop(-grad + grad_norm*rho0)# Descent direction
    gamma <- gamma0 # Step size
    par1 <- par + gamma * rho   # Proposed descent step
    h_prime <- drop(t(grad) %*% rho) # inner prodict of gradient and decent direction. 
    if(!is.null(cb)) cb() # Tracing 
    
    if (m > length(par) || h_prime >= 0)
    {
      rho <- -grad
      m <- 1
      h_prime <- grad_norm
    }
    
    rho0 <- rho /grad_norm
    
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


#############optim

# fix  x, y, lamda and inner_knots and choose initial guess p0

H_optim <- function(par){
  optim(par, function(par)  H(par, x, y, lambda, inner_knots))
}


#Newton probabilites
pvaerdier_Newton <- function(x) {
  drop(invlogit(f(par_New, x, inner_knots)))
}

#Gradient probabilites
pvaerdier_GD <- function(x) {
  drop(invlogit(f(par_GD, x, inner_knots)))
}

#Conjugent gradient probabilites
pvaerdier_CG <- function(x) {
  drop(invlogit(f(par_CG, x, inner_knots)))
}

#optim probabilities
pvaerdier_optim <- function(x) {
  drop(invlogit(f(par_optim, x, inner_knots)))
}


