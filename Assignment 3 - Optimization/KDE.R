library(ggplot2)
library(microbenchmark)
library(Rcpp)
library(RcppArmadillo)
setwd("~/Documents/Statistics - KU/Computational Statistics/Assignments/Assignment 1 - Density smoothing")
source("Assignment-1-functions.R")

# simulated data.
x <- rnorm(10000)
h <- 1
par(mfrow = c(2,2))

hist(x, probability = T, col = 'lightblue', main = 'loop')
lines(kdens_loop(x,h), lwd = 2, col = '#FF7F50')

hist(x, probability = T, col = 'lightblue', main = 'sapply')
lines(kdens_sapply(x, h), lwd = 2, col = '#B22222')

hist(x, probability = T, col = 'lightblue', main = 'binning')
lines(kdens_bin(x,h), lwd = 2, col = '#0000CD')

hist(x, probability = T, col = 'lightblue', main = 'density')
lines(density(x, bw = h, kernel = "epanechnikov"), lwd = 2, col = 'green') # rescale h such that it fits

F12 <- log(read.table("infrared.txt", header = T)$F12)

h <- 1

par(mfrow = c(2,2))

hist(F12, probability = T, col = 'lightblue', main = 'loop')
lines(kdens_loop(F12,h), lwd = 2, col = '#FF7F50')

hist(F12, probability = T, col = 'lightblue', main = 'sapply')
lines(kdens_sapply(F12, h), lwd = 2, col = '#B22222')

hist(F12, probability = T, col = 'lightblue', main = 'binning')
lines(kdens_bin(F12,h), lwd = 2, col = '#0000CD')

hist(F12, probability = T, col = 'lightblue', main = 'density')
lines(density(F12, bw = h, kernel = "epanechnikov"), lwd = 2, col = 'green') # rescale h such that it fits

x <- rnorm(10000)

kern_bench = microbenchmark(times = 100,
                            kdens_loop(x, h),
                            kdens_sapply(x, h),
                            kdens_bin(x, h)
)

summary(kern_bench)

autoplot(kern_bench) + geom_jitter(aes(color = expr), alpha = 0.4) + aes(fill = I("gray"))+
  scale_color_manual(values = c("#354e75", "#884091", "#c91246", "black")) + theme_bw() +
  theme(legend.position = "none")
