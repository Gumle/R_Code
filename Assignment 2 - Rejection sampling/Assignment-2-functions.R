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

# Density functions
target_sapply <- function(y){sapply(y, function(yy) prod(exp(yy*zx - exp(x*yy)))) / 1.05009e-41 }

target_sum <- function(y){sapply(y, function(yy) exp( sum(yy*zx - exp(x*yy)) )) / 1.05009e-41}

target_outer <- function(y){exp( rowSums(outer(y, zx) - exp(outer(y, x))))/ 1.05009e-41}

# Function factory for our proposal density
proposal_factory <- function(target_dens, proposal_dens, proposal_sim, a){
  list( tdens = target_dens, pdens = proposal_dens, sim = proposal_sim, alpha = a)
}

# Rejection sampling - fast

target_simulator <- function(n, proposal_object, scale ){
  
  simulated_values <- numeric(n)  #we need n samples
  
  num_accepted <- 0 #we keep track of how many we have
  
  while (num_accepted < n){
    samples <- ceiling( scale*(n - num_accepted)/proposal_object$a ) #sample expected amount
  
    prop_samples <- proposal_object$sim(samples)  #get all the samples
    
    uniform_samples <- runif(samples)
    
    accept <- uniform_samples <=  proposal_object$a*proposal_object$tdens(prop_samples) /proposal_object$pdens(prop_samples)
    
    new_additions <- min(n - num_accepted, sum(accept))
    
    simulated_values[(num_accepted + 1):(num_accepted + new_additions)] <- (prop_samples[accept])[1:new_additions]
    
    num_accepted <- num_accepted + new_additions
  }
  simulated_values
}


# Rejection sampling - slow

target_simulator_slow <- function(n, proposal_object) {
  #browser()
  y <- numeric(n)
  for(i in 1:n) {
    reject <- TRUE
    while(reject) {
      y0 <- proposal_object$sim(1)
      u <- runif(1)
      reject <- u > proposal_object$a*proposal_object$tdens(y0) / proposal_object$pdens(y0)
    }
    y[i] <- y0
  }
  y
}


envelope_factory <- function(x,
                             target_dens,
                             logderiv = NULL,
                             lower_support = -Inf,
                             upper_support = Inf){
  
  #check if the log-derivative was supplied, otherwise we get it ourselves:
  if (is.null(logderiv)){
    logderiv <- function(xx) grad(function(x) log(target_dens(x)), xx)}
  
  #calculate a-vector and stop if there are integrability issues
  a <- logderiv(x)
  
  #check if we need to stop
  continue <- (a[1] > 0 & a[length(a)] < 0) | (a[length(a)] < 0 & lower_support > -Inf)|
    a[1] > 0 & upper_support < Inf
  
  if (!continue) {stop("Envelope is not integrable. Re-submit new x")}
  
  #now calculate b, z, Fz, Q, and const (c) 
  
  b <- log(target_dens(x)) - a*x
  z <- c(lower_support, -diff(b)/diff(a), upper_support)
  Fz <- numeric(length(x))
  
  for (i in seq_along(Fz))
  {
    Fz[i] <- exp(b[i]) * ( exp( a[i] * z[i+1]) - exp(a[i]*z[i])) / a[i]
  }
  
  Q <- c(0, cumsum(Fz))
  const <- Q[length(Q)]
  
  #now define simulator and envelope -- must be vectorized
  
  proposal_density <- function(x) {
    index <- findInterval(x, z)  #Given x, find z-interval that x belongs to
    
    #findInterval is optimized and O(log(length(z)*length(x))
    
    #now just evaluate the function and return
    V <- a[index]*x + b[index]
    exp(V)/const
  }
  
  simulator <- function(n){
    
    #we need n uniform
    q <- runif(n)
    
    #find the index to which const*q belongs to.
    index <- findInterval(const*q, Q) 
    
    #solve the equation for x
    log( (const*q - Q[index])*a[index ]*exp(-b[index])+ exp(a[index]*z[index]))/a[index ]}
  
  #return the list-object that can be used directly in previous functions.
  
  list( tdens = target_outer, pdens = proposal_density, sim = simulator, alpha = 1/const)
}

grid <- seq(0, 0.5, by = 0.001)
pois <- read.csv("poisson.csv", header = T)
zx <-  pois$z*pois$x
x <- pois$x
grid <- seq(0, 0.5, by = 0.001)
propn <- proposal_factory(target_outer, function(x) dnorm(x, mean, sd*1.09), function(x) rnorm(x, mean, sd*1.09), alpha)
target_simulator_slow(10000, propn)
grid <- seq(0, 0.5, length.out = 200)
envelope_factory(grid, logderiv = target_outer)
microbenchmark(samples <- target_simulator(10000, propn, 1),
envelope_factory(grid, logderiv = target_outer),times = 100)
