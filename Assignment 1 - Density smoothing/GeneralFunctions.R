library(MASS)
library(stats)
library(microbenchmark)
library(ggplot2)
library(kedd)
library(caret)
library(Rcpp)
library(RcppArmadillo)
setwd("~/Documents/Statistics - KU/Computational Statistics/Assignments/Assignment 1 - Density smoothing")
infrared <- read.table("infrared.txt", header = TRUE)
F12 <- infrared$F12
F12 <- log(F12)

EKernel <- function(x){
  ifelse(-1 <= x & x<= 1,(3/4)*(1-x^2), 0) 
}

t <- seq(-2,2,0.1)
sum(EKernel(t-1))

kernDens_loop <- function(x, h, m = 512) {
  
  rg <- range(x) + c(-3 * h, 3 * h)
  xx <- seq(rg[1], rg[2], length.out = m) # evaluation grid around x : m points
  y <- numeric(m) # vector to store the function values - same length as xx
  const <- 3/(4* h*length(x))
  for (i in seq_along(xx)){
    condition <- xx[i]- x
    kde <- ifelse(abs(condition) <= h, 1 - (condition/h)^2, 0) 
    y[i] <- sum(kde)*const
  }
  list(x = xx, y = y) # return evauation grid and the function values
}
t1 <- kernDens_loop(F12,h = 1)
round(sum(t1$y))
EKernel <- function(x){ifelse(-1 <= x & x<= 1,x <- (3/4)*(1-x^2), 0) }
EKernel(F12)


EKernel <- function(x){
  ifelse(-1 <= x & x<= 1,(3/4)*(1-x^2), 0) 
}

kernDens_loop2 <- function(x, h, m = 512, kernel) {
  rg <- range(x) + c(-3 * h, 3 * h)
  xx <- seq(rg[1], rg[2], length.out = m) # evaluation grid around x : m points
  y <- numeric(m) # vector to store the function values - same length as xx
  for (i in seq_along(xx)){
    kde <- kernel((xx[i]- x)/h)
    y[i] <- sum(kde)/(h*length(x))
  }
  list(x = xx, y = y) # return evauation grid and the function values
}
xinit <- seq(-4,4,0.5)
t2 <- kernDens_loop2(F12, h = 1, kernel = EKernel)
sum(t2$y)
################## GENERAL ##############

kernDens_apply <- function(x, h, m = 512){
  rg <- range(x) + c(-3 * h, 3 * h)
  xx <- seq(rg[1], rg[2], length.out = m) # evaluation grid around x : m points 
  y <- numeric(m) # vector to store the function values - same length as xx
  const <- (3/4) / (h * length(x))
  y <- sapply(xx, function(z) sum((1 - ((z - x) / h)^2 ) * (abs(x-z) <= h)))
  y <- const * y
  list(x = xx, y = y) # return evauation grid and the function values
}
old <- kernDens_apply(F12, h = 1)

#More generalized
kernDens_apply <- function(x, h, m = 512, kernel){
  rg <- range(x) + c(-3 * h, 3 * h)
  xx <- seq(rg[1], rg[2], length.out = m) # evaluation grid around x : m points 
  y <- numeric(m) # vector to store the function values - same length as xx
  y <- sapply(xx, function(z) sum(kernel((x-z)/h)))
  y <- y/(h * length(x))
  list(x = xx, y = y) # return evauation grid and the function values
}
new <- kernDens_apply(F12,h = 1, kernel = EKernel)
max(abs(old$y-new$y))


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

n1 <- kernDens_bin(F12,h = 1)

kernDens_bin2 <- function (x, h, m = 512, kernel) {
  rg <- range(x) + c(-3 * h, 3 * h)
  xx <- seq(rg[1], rg[2], length.out = m) # evaluation grid around x : m points 
  weights <- bin(x, rg[1], rg[2], m) # compute weights 
  kerneval <- kernel((xx - xx[1])/h)/h
  kerndif <- toeplitz(kerneval)
  y <- colSums(weights * kerndif)
  list(x = xx, y = y) # return evaluation grid and function values 
}
n2 <- kernDens_bin2(F12,h = 1, kernel = EKernel)
max(abs(n1$y-n2$y))
x1 <- 1:5
toeplitz (x1)


###### CV #######

cv <- function (x, h,kernel = EKernel){
  n <- length(x)
  y <- 0
  for (i in 1:n){
    y[i] <- sum(kernel((x[i] - x[-i])/h)) * 1/(h * (n-1))
  }
  value <- sum(log(y[y !=0])) 
  value
}
h <- seq(1, 5, 0.05)
h_cv <- sapply(h, cv(x,h))
h_opt <- h[which.min(CV)]
qplot(h, CV) + geom_line() + geom_vline(xintercept = h_opt, color = "red")
cv(F12, h=1)
kernel <- EKernel
cv.factory <- function(x) {
  f <- function(h) cv(x, h)
  #Vectorize(f)
}
f <- cv.factory(F12)
h.cv <- optimize(f, interval = c(0.001, 4), maximum = TRUE)$maximum


Kgaybar <- function(x){
  ifelse(0 <= x & x <= 2, 3*((2-x)^3)*(x^2+6*x+4)/160,0)
} 


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

f2norm.factory <- function(x){
  function(h) f2norm(x,h)
}
f2norm.data <-f2norm.factory(F12)
fSum <- function (x, h, kernel = EKernel){
  n <- length(x)
  y <- 0
  for (i in 1:n){
    y[i] <- sum(kernel((x[i] - x[-i])/h)) * 1/(h * (n-1))
  }
  value <- 2*sum(y)/n
  value
}
fSum.data <- fSum.factory(F12)
UCV <- function(h){
  f2norm.data(h) - fSum.data(h)
}

h.opt <- optimize(UCV, interval = c(0.0001,2))$minimum
h.opt
