library(ggplot2)
library(microbenchmark)
library(Rcpp)
library(RcppArmadillo)
library(kedd)

setwd("~/Documents/Statistics - KU/Computational Statistics/Assignments/Assignment 1 - Density smoothing")
dat <- read.table("infrared.txt", header = TRUE, sep = "")
F12 <- log(dat$F12)
source("Assignment-1-functions.R") # Get the functions

normcalculator <- norm2factory(F12)

normcalculator(0.5)

LOOCVsum(F12, 0.2)

sumofleaveout <- LOOCVfac(F12)

min <- optimize(ucv, interval = c(0.0001, 4))$minimum

density(F12, bw ="ucv", kernel = "gaussian")

sourceCpp("cpp.code.cpp")

# Simulated data:

x <- rnorm(10000)
h <- 0.8
loop <- kernDens_loop(x,h)
apply <- kernDens_loop(x,h)
bin <- kernDens_loop(x,h)
# Are the implementations the same?
max(kernDens_loop(x,h)$x-kernDens_apply(x,h)$x)
max(kernDens_loop(x,h)$y-kernDens_apply(x,h)$y)
max(kernDens_apply(x,h)$x-kernDens_bin(x,h)$x)
max(kernDens_apply(x,h)$y-kernDens_bin(x,h)$y)
max(kernDens_bin(x,h)$y-density(F12, bw =h, kernel = "epanechnikov")$y)
# Kernel density estimation:

par(mfrow = c(1,2))
hist(x, probability = T, col = 'dodgerblue4', main = 'binning')
lines(kernDens_bin(x,h), lwd = 2, col = 'chocolate1')

hist(x, probability = T, col = 'dodgerblue4', main = 'density')
lines(density(x, bw = h, kernel = "epanechnikov"), lwd = 2, col = 'chocolate1') # rescale h such that it fits
par(mfrow=c(1,1))
# Difference:
plot(kernDens_bin(x,h)$x,abs(kernDens_bin(x,h)$y-density(x, bw = h, kernel = "epanechnikov")$y), type = "l", ylab = "Absoulute difference between binning", xlab = "x")

