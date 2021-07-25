setwd("~/Documents/Statistics - KU/Computational Statistics/Assignments/Assignment 3 - Optimization")
library(splines)
library(smooth)
library(microbenchmark)
data <- read.csv("Horses.csv", sep = ",", header = TRUE)
data <- na.omit(data)

knots <- function(inner_knots) {c(rep(range(inner_knots), 3), inner_knots)}
inner_knots <- data$Temperature
par <- rnorm(534, mean = 0, sd = 1)

# f
f <- function(par, x, inner_knots) {
  if(length(par) != length(inner_knots) + 2) {stop()}
  X <- splineDesign(knots(inner_knots), x)
  X %*% par
}
f(par, data$Temperature, inner_knots)
# H

H <- function(par, x, y, lambda, inner_knots) {
  f_val <- f(par, x, inner_knots)
  obj_func <- -sum(y * f_val - log(1 + exp(f_val)))
  penalty <- lambda * crossprod(par, pen_mat(inner_knots) %*% par)
  (obj_func + penalty) / length(z)
}

# H'

grad_H <- function(par, x, y, lambda, inner_knots) {
  phi <- splineDesign(knots(inner_knots), x)
  f_val <- f(par, x, inner_knots)
  p <- logit_inv(f_val) -((crossprod(phi, y - p) - lambda * pen_mat(inner_knots) %*% par) / length(z))
}

#H''

hessian_H <- function(par, x, y, lambda, inner_knots) {
  phi <- splineDesign(knots(inner_knots), x)
  f_val <- f(par, x, inner_knots)
  p <- logit_inv(f_val)
  W <- diag(as.vector(p * (1 - p))) -((-crossprod(phi, W %*% phi) - lambda * pen_mat(inner_knots)) / length(x))
}

#Evaluate the design matrix for the B-splines defined by knots at the values in x.
inner_knots <- data$Temperature

X <- splineDesign(knots(inner_knots), inner_knots)
par <- as.matrix(par)
t <- X %*% par
dim(X)
dim(par)
sum(X != 0)
dim(X)[1] * dim(X)[2] - sum(X != 0)


