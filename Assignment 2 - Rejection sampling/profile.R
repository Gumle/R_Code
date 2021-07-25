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

grid <- seq(0, 0.5, by = 0.001)
y <- pois
propn <- proposal_factory(target_outer, function(x) dnorm(x, mean, sd*1.09), function(x) rnorm(x, mean, sd*1.09), alpha)
profvis(target_simulator_slow(10000, propn))
profvis(samples <- target_simulator(10000, propn, 1))
grid <- seq(0.1, 0.5, length.out = 1000)
# Data
pois <- read.csv("poisson.csv", header = T)
zx <-  pois$z*pois$x
x <- pois$x
profvis(envelope_factory(grid, logderiv = target_outer))
profvis(tmp <- envelope_factory(grid))

        