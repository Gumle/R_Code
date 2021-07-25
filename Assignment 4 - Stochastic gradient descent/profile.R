#p is first, mu1 is second, mu2 is third, sigma1 is fourth, sigma2 is fifth, nu1 is sixth, nu2 is seventh 
library(MASS)
library(stats)
library(microbenchmark)
library(ggplot2)
library(kedd)
library(caret)
library(metRology)
library(numDeriv)
library(parallel)
library(Rcpp)
library("profvis")

setwd("~/Documents/Statistics - KU/Computational Statistics/Assignments/Assignment 4 - Stochastic gradient descent")
source("Debugging_and_tracing.R")
source("loglikegrad.R")
source("Assignment-4-functions.R")
source("EM-Algorithm.R")
par <-c(0.3, 10, 7, 3, 0.5, 2, 4)
data <- data.sample(n = 10000, par = par)
par.init <- c(0.2, 12, 6, 4, 0.7, 2.4, 4.9)

profvis(emalgoritm(minusloglike, Mstepfunction, derivMstep, par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9),  epsilon = 10 ^ (-6),
                   x = data, gamma = 0.00006))

SGDFast.scale <- function(par, H, grad_H, x, divlength, gamma, epsilon = 10^(-6), scale.p, scale.rest,  cb = NULL){
  samplesize <- length(x)
  j_end <- floor(samplesize/divlength) 
  repeat{
    value <- H(par, x)
    index <- sample(1: samplesize, replace = FALSE)
    y <- x[index]
    for (j in 1: j_end){
      grad.0 <- grad_H(par, y[(1 + (j - 1) * divlength ): (divlength * j )])
      par[2:7] <- par[2:7] - gamma * scale.rest * grad.0[2:7]
      par[1] <- par[1] - (gamma*scale.p) * grad.0[1]
    }
    if(!is.null(cb)) cb()
    criteria <- value - H(par, x) - epsilon * (H(par,x) + epsilon)
    if(criteria <= 0) {
      break
    }
  }
  par
}



#profvis(SDG(par = par, x = data, divlength = 50, maxit = 100, gamma = 0.0001))
profvis(SGDFast.scale(par = par.init, H = logLike, grad_H = analyticGrad,scale.p = 1/10, scale.rest = 3, x = data, divlength = 100,  gamma = 0.005))

profvis(GD(par = par.init, x = data, H = logLike, grad_H = analyticGrad, gamma0 = 0.00005))

