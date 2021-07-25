# Intro and packages
library(purrr)
library(microbenchmark)
library(Rcpp)
library(ggplot2)
library("profvis")
library(MASS)
library(parallel)
library(numDeriv)
library(doParallel)
setwd("~/Documents/Statistics - KU/Computational Statistics/Assignments/Assignment 2 - Rejection sampling")
source("Assignment-2-functions.R")

# Data
pois <- read.csv("poisson.csv", header = T)
zx <-  pois$z*pois$x
x <- pois$x
grid <- seq(0, 0.5, by = 0.001)

#Benchmark
p <- microbenchmark(times = 200, target_sum(grid), target_sapply(grid), target_outer(grid))
autoplot(p) + geom_jitter(aes(color = expr), alpha = 0.4) + aes(fill = I("gray"))+
  scale_color_manual(values = c("#354e75", "#884091", "#c91246", "black")) + theme_bw() +
  theme(legend.position = "none", text = element_text(size=15))
p
b <- "dodgerblue4"
o <- "chocolate1"
# Finding a gaussian envelope
mean <- integrate(function(v) target_outer(v)*v, 0, Inf)$value
sd <- sqrt(integrate(function(v) target_outer(v)*v^2, 0, Inf)$value - mean^2)
curve(target_outer(x), 0, 0.5, lwd = 3, col ="dodgerblue4", main = " Target density and gaussian density", ylab = "f(y), g(y)", xlab = "y")
curve(dnorm(x, mean, sd), 0, 0.5, lwd = 3, col ="chocolate1", add =T)

# Target density and gaussian density
s <- function(y) dnorm(y, mean, sd)/target_outer(y)

alpha <- optimize(s, c(0,0.1))$objective # global minimum in dm(s)
alpha
#First alpha:
curve(s(x), from = 0, 0.5, lwd = 4, col = "dodgerblue4", main = "g(y)/f(y)", xlab = "y", ylab = "alpha" )
abline(h = alpha, col = "chocolate1", lwd = 2)

#First envelope:
curve(alpha*target_outer(x), 0, 0.5, lwd = 3, col ="dodgerblue4", main = "Gaussian envelope", ylab = "", xlab = "", ylim = c(0,8))
curve(dnorm(x, mean, sd), 0, 0.5, lwd = 3, col ="chocolate1", add =T)

#Improved envelope:
curve(target_outer(x), 0, 0.5, lwd = 3, col ="dodgerblue4", main = " Target density and gaussian density", ylab = "f(y), g(y)", xlab = "y")
curve(dnorm(x, mean, sd*1.09), 0, 0.5, lwd = 3, col ="chocolate1", add =T)
s <- function(y) dnorm(y, mean, sd*1.09)/target_outer(y)
alpha.improved <- optimize(s, c(0,0.4))$objective # global minimum in dm(s)
alpha.improved

#Improved alpha:

curve(s(x), from = 0, 0.5, lwd = 4, col = "dodgerblue4", main = "g(y)/f(y)", xlab = "y", ylab = "alpha" ,ylim = c(0,2))
abline(h = alpha.improved, col = "chocolate1", lwd = 3)

#Improved envelope:
curve(alpha.improved*target_outer(x), 0, 0.5, lwd = 3, col ="dodgerblue4", main = "Improved Gaussian envelope", ylab = "", xlab = "", ylim = c(0,8))
curve(dnorm(x, mean, sd*1.09), 0, 0.5, lwd = 3, col ="chocolate1", add =T)

########## - Rejection sampling algorithms #############

propn <- proposal_factory(target_outer, function(x) dnorm(x, mean, sd*1.09), function(x) rnorm(x, mean, sd*1.09), alpha)

samples <- target_simulator(10000, propn, 1)
samples_slow <- target_simulator_slow(10000, propn)
sum(abs(samples-samples_slow))

par(mfrow=c(1,2))

hist(samples, probability = T, breaks = 20, col = "dodgerblue4", main = "samples target_simulator_slow ")
curve(target_outer(x), add = T, col = "chocolate1", lwd = 2)
hist(samples, probability = T, breaks = 20, col = "dodgerblue4", main = "samples target_simulator")
curve(target_outer(x), add = T, col = "chocolate1", lwd = 2)

# Benchmarking
grid <- seq(0.01, 0.5, length.out =10000)
grid1 <- seq(0.01, 0.5, length.out = 4)
fastRS <- function(propn) target_simulator(10000, propn, 1)
slowRS <- function(propn) target_simulator_slow(10000, propn)
AdaptEnv <-function(grid) envelope_factory(grid, target_outer)
AdaptEnv(grid)
p <- microbenchmark(times = 200, fastRS(propn), AdaptEnv(grid))
library(knitr)
summary(p)
autoplot(p) + geom_jitter(aes(color = expr), alpha = 0.4) + aes(fill = I("gray"))+
  scale_color_manual(values = c("#354e75", "#884091", "black")) + theme_bw() +
  theme(legend.position = "none",text = element_text(size=10))


### Adaptive envelope

grid <- seq(0.01, 0.5, length.out = 20) # points where to place picewise functions
length(grid)
tmp <- envelope_factory(grid, target_outer)
str(tmp)
samples <- target_simulator(10000, tmp, 1)
par(mfrow=c(2,2))
hist(samples, breaks = 20, probability = T)
curve(tmp$pdens(x), 0.1, 0.5, col = "red", lwd = 2, add = T)
curve(target_outer, 0, 0.5,add = T)

grid3 <- seq(0.1, 0.5, length.out = 3)
tmp3 <- envelope_factory(grid3, target_outer)
tmp3$alpha
grid4 <- seq(0.1, 0.5, length.out = 4)
tmp4 <- envelope_factory(grid4, target_outer)
tmp4$alpha
grid5 <- seq(0.1, 0.5, length.out = 5)
tmp5 <- envelope_factory(grid5, target_outer)
tmp5$alpha
grid6 <- seq(0.1, 0.5, length.out = 6)
tmp6 <- envelope_factory(grid6, target_outer)
tmp6$alpha
grid7 <- seq(0.1, 0.5, length.out = 7)
tmp7 <- envelope_factory(grid7, target_outer)
tmp7$alpha
grid8 <- seq(0.1, 0.5, length.out = 8)
tmp8 <- envelope_factory(grid8, target_outer)
tmp3$alpha
tmp4$alpha
tmp5$alpha
tmp6$alpha
tmp7$alpha
tmp8$alpha

par(mfrow=c(2,3))
hist(samples, breaks = 20, probability = T, main="alpha = 0.5352192")
curve(tmp3$pdens(x), 0.1, 0.5, col = "red", lwd = 2, add = T)
curve(target_outer, 0, 0.5,add = T)
hist(samples, breaks = 20, probability = T,main="alpha = 0.7788359")
curve(tmp4$pdens(x), 0.1, 0.5, col = "red", lwd = 2, add = T)
curve(target_outer, 0, 0.5,add = T)
hist(samples, breaks = 20, probability = T,main="alpha = 0.8703055")
curve(tmp5$pdens(x), 0.1, 0.5, col = "red", lwd = 2, add = T)
curve(target_outer, 0, 0.5,add = T)
hist(samples, breaks = 20, probability = T,main="alpha = 0.9156727")
curve(tmp6$pdens(x), 0.1, 0.5, col = "red", lwd = 2, add = T)
curve(target_outer, 0, 0.5,add = T)
hist(samples, breaks = 20, probability = T,main="alpha = 0.9404132")
curve(tmp7$pdens(x), 0.1, 0.5, col = "red", lwd = 2, add = T)
curve(target_outer, 0, 0.5,add = T)
hist(samples, breaks = 20, probability = T,main="alpha = 0.9554017")
curve(tmp8$pdens(x), 0.1, 0.5, col = "red", lwd = 2, add = T)
curve(target_outer, 0, 0.5,add = T)

grid <- seq(0, 0.5, length.out = 10000)
tmp <- envelope_factory(grid, target_outer)
tmp$alpha
par(mfrow=c(1,2))
hist(samples, breaks = 40, probability = T,main="Gaussian envelope")
curve(target_outer, 0, 0.5, col = "dodgerblue4", lwd = 2, add = T)
hist(samples, breaks = 40, probability = T,main="Adaptive envelope")
curve(tmp$pdens(x), 0.1, 0.5, col = "firebrick4", lwd = 2, add = T)
##### COMPARISON ########



AdaptEnv <- function(grid) envelope_factory(grid, target_outer)
grid <- seq(0, 0.5, length.out = 20000)
ad <- AdaptEnv(grid)
grid0 <- seq(0, 0.5, length.out = 200)
grid3 <- seq(0, 0.5, length.out = 50000)
propn <- proposal_factory(target_outer, function(x) dnorm(x, mean, sd*1.09), function(x) rnorm(x, mean, sd*1.09), alpha)
samples <- target_simulator(10000, propn, 1)

samples_slow <- target_simulator_slow(10000, propn)
grid <- seq(0, 0.5, length.out = 200)
samples <- target_simulator(10000, envelope_factory(grid, target_outer), 1)

FastRS <- function(propn) target_simulator(10000,propn,1)

p <- microbenchmark(times = 200, FastRS(propn), AdaptEnv(grid))
summary(p)
autoplot(p) + geom_jitter(aes(color = expr), alpha = 0.4) + aes(fill = I("gray"))+
  scale_color_manual(values = c("#354e75", "#884091", "black")) + theme_bw() +
  theme(legend.position = "none",text = element_text(size=10))



# Object-Oriented programming
proposal_factory <- function(target_dens, proposal_dens, proposal_sim, a)
{
  structure(list( tdens = target_dens,
                  pdens = proposal_dens,
                  sim = proposal_sim,
                  alpha = a), class = "RS")
}
plot.RS <- function(object, support)
{
  gridx <- seq(support[1], support[2], length.out = 512)
  plot(gridx, object$alpha*object$tdens(gridx), col = "#286e99", xlab = "", ylab = "",
       main = paste("Alpha equals", propn$alpha), type = "l", lwd = 2)
  points(gridx, object$pdens(gridx), col = "#a89e82")
}
simulate.RS <- function(object, nsim = 1, scale = 1, plot = FALSE)
{
  simulations <- target_simulator(nsim, object, scale)
  if (plot)
  {
    hist(simulations, prob = TRUE, main = "", ylab = "", xlab = "", col = "#a89e82")
    curve(object$tdens(x), add = TRUE)
  }
  simulations
}
plot.RS(propn)
