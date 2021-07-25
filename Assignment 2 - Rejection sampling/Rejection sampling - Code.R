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
set.seed(10)

# Read in data and vectors:

grid <-seq(0, 0.5, 0.001)
pois <- read.table("Poisson.csv", sep = ",", header = TRUE)
zx <- pois$z*pois$x
x <- pois$x


# Implement the target density

target_sapply <- function(y){sapply(y, function(yy) prod(exp(yy*zx - exp(x*yy)))) / 1.05009e-41 }

target_sum <- function(y){sapply(y, function(yy) exp( sum(yy*zx - exp(x*yy)) )) / 1.05009e-41}

target_outer <- function(y){exp( rowSums(outer(y, zx) - exp(outer(y, x))))/ 1.05009e-41}

y <- grid
length(y)
length(z)
out <- exp( rowSums(outer(y, zx) - exp(outer(y, x))))
out.test <- outer(y,zx)
sap <- sapply(y, function(yy) exp( sum(yy*zx - exp(x*yy)) ))
# Compare the fastest target density implementation:

p <- microbenchmark(times = 100, target_sum(grid), target_sapply(grid), target_outer(grid))
autoplot(p) + geom_jitter(aes(color = expr), alpha = 0.4) + aes(fill = I("gray"))+
  scale_color_manual(values = c("#354e75", "#884091", "#c91246", "black")) + theme_bw() +
  theme(legend.position = "none")

# Target density:

plot(grid,target_outer(grid), type = "l", lwd = 4, col = "#325b96", ylab = "Density", xlab = "")
ggplot() + geom_point(data = NULL, aes(x = grid, y = target_outer(grid)), color = "#325b96") + ylab("Density") + xlab ("") 

# We find the gaussian density using the dnorm, using our mean and sd when integrating our target density
## Gaussian density

mean <- integrate(function(v) target_sapply(v)*v, 0, Inf)$value
sd <- sqrt(integrate(function(v) target_sapply(v)*v^2, 0, Inf)$value - mean^2)
mean
sd
targetDens <- target_outer(grid)
gausEnv <- dnorm(grid, mean, sd)

# Ratio plot: g/f (gausEnv/targetDens)

plot(grid, targetDens, type = "l", lwd = 4, col = "#325b96", ylab = "Target density", xlab = "")
lines(grid, gausEnv, col = "#ee961d", lwd =3)
ggplot() + geom_line(data = NULL, aes(x = grid, y = target_outer(grid)), color = "#325b96", lwd = 2)+ geom_line(data = NULL, aes(x = grid, y = gausEnv), color = "#ee961d", lwd = 2) + ylab("Density") + xlab ("") 

# We now want to find our alpha-parameter so we need to find the minimum of 
plot(grid, gausEnv/targetDens, lwd =1, col = "#325b96" ,type = "l")
ggplot() + geom_line(data = NULL, aes(x = grid, y = gausEnv/targetDens),color = "#325b96", lwd = 2)

# Minimizing the ratio so we get the smallest but closest to 1 as possible
# And by using the optimize function we get that our alpha parameter is 0.2378
a <- optimize(function(z) dnorm(z, mean, sd)/target_outer(z), c(0, 0.1))
a <- a$objective

# Our gaussian envelope and our envelope:

plot(grid, a*target_outer(grid), type = "l", lwd = 3,  col = "#325b96")
lines(grid, dnorm(grid, mean, sd),  col = "#ee961d", lwd =3)

# And clearly this needs some improvement.

# Improved envelope
gausEnv2 <- dnorm(grid, mean, sd*1.09)
plot(grid, targetDens, type = "l", lwd = 4, col = "#325b96", ylab = "Target density", xlab = "")
lines(grid, gausEnv2, col = "#ee961d", lwd =3)

a <- optimize(function(z) dnorm(z, mean, sd*1.09)/target_outer(z), c(0, 0.1))
a <- a$objective

# First implementation
target_dens <- target_outer(grid)
proposal_dens <- dnorm(grid, mean, sd*1.09)
proposal_sim <- dnorm(grid, mean, sd*1.09)
proposal_sim(1)

# Proposal factory
#generic implementation

proposal_factory <- function(target_dens, proposal_dens, proposal_sim, a)
{
  list( tdens = target_dens, pdens = proposal_dens, sim = proposal_sim, alpha = a)
}
propn <- proposal_factory(target_outer, function(x) dnorm(x, 0.2388885, 0.061235426), function(x) rnorm(x, 0.2388885, 0.06123542), 0.9163125)
#propt <- proposal_factory(target_outer,function(x) dt((x - 0.24)/(0.0554*1.055), 61.83)/(0.0554*1.055),function(x) rt(x, 61.83)*0.0554*1.055 + 0.24,0.95)

# Slow implementation
target_simulator_slow <- function(n, proposal_object) {
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
estimate1 <- target_simulator_slow(10000, propn)
hist(estimate1, breaks = 20, probability = TRUE)
lines(grid, targetDens, col = "#325b96", lwd =3)
target_simulator <- function(n, proposal_object, scale )
{
  simulated_values <- numeric(n)  #we need n samples
  
  num_accepted <- 0 #we keep track of how many we have
  
  while (num_accepted < n)
  {
    samples <- ceiling( scale*(n - num_accepted)/proposal_object$a ) #sample expected amount
    
    prop_samples <- proposal_object$sim(samples)  #get all the samples
    uniform_samples <- runif(samples)
    
    accept <- uniform_samples <=  proposal_object$a*proposal_object$tdens(prop_samples) /
      proposal_object$pdens(prop_samples)
    
    new_additions <- min(n - num_accepted, sum(accept))
    
    simulated_values[(num_accepted + 1):(num_accepted + new_additions)] <- (prop_samples[accept])[1:new_additions]
    
    num_accepted <- num_accepted + new_additions
  }
  simulated_values
}
estimate2 <- target_simulator(10000, propn, 1)

hist(estimate2, breaks = 20, probability = TRUE)
lines(grid, targetDens, col = "#325b96", lwd =3)
# This does not work.
library(doParallel)
foreach(i = rep(10000/8, 8),  .combine = 'c') %dopar% {
    propn <- proposal_factory(target_outer, function(x) dnorm(x, mean, sd*1.09), function(x) rnorm(x, mean, sd*1.09), alpha)
    target_simulator(i, propn, 1)    
}

kern_bench <- microbenchmark(
  target_simulator_slow(10000, propn),
  target_simulator(10000, propn, 1)
)
autoplot(kern_bench) + 
  geom_jitter(position = position_jitter(0.2, 0), 
              aes(color = expr), alpha = 0.4) + 
  aes(fill = I("gray")) + 
  theme(legend.position = "none")
