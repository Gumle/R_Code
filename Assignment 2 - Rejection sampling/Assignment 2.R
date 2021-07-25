library(purrr)
library(microbenchmark)
library(Rcpp)
library(ggplot2)
library("profvis")
library(MASS)
library(parallel)
library(numDeriv)
setwd("~/Documents/Statistics - KU/Computational Statistics/Assignments")
grid <- seq(0,0.5,0.001)
pois <- read.table("Poisson.csv", sep = ",", header = TRUE)

# How to evaluate the target density

z <- pois$z
x <- pois$x

# Unnormalized

target_sapply <- function(y){
  sapply(y, function(yy) prod(exp(yy*zx - exp(x*yy))))
}

target_sum <- function(y){
  sapply(y, function(yy) exp( sum(yy*z*x - exp(x*yy)) ))
}

target_outer <- function(y){
  exp( rowSums(outer(y, z*x) - exp(outer(y, x))))
}

t <- integrate(target_sapply, 0, Inf)
value_sapply <- integrate(target_sapply, 0, Inf)$value
value_sum <- integrate(target_sum, 0, Inf)$value
value_outer <- integrate(target_outer, 0, Inf)$value

# Normalized

target_sapply <- function(y)
{
  sapply(y, function(yy) prod(exp(yy*zx - exp(x*yy)))) / value_sapply
}
target_sum <- function(y)
{
  sapply(y, function(yy) exp( sum(yy*zx - exp(x*yy)) )) / value_sum
}
target_outer <- function(y)
{
  exp( rowSums(outer(y, zx) - exp(outer(y, x))))/ value_outer
}

integrate(target_outer, 0, Inf)
grid <-seq(0, 0.4, 0.001)
target_outer(grid)

# See which method is the best.
p <- microbenchmark(times = 200, target_sum(grid), target_sapply(grid), target_outer(grid))
autoplot(p) + geom_jitter(aes(color = expr), alpha = 0.4) + aes(fill = I("gray"))+
  scale_color_manual(values = c("#354e75", "#884091", "#c91246", "black")) + theme_bw() +
  theme(legend.position = "none")


####################

grid <-seq(0, 0.4, 0.0001)
mean <- integrate(function(v) target_sapply(v)*v, 0, Inf)$value
sd <- sqrt(integrate(function(v) target_sapply(v)*v^2, 0, Inf)$value - mean^2)

mean
sd
#now we can plot them
grid <- seq(0, 0.5, 0.001)

plot(grid,target_outer(grid), type = "l", lwd = 4, col = "#325b96", ylab = "Density", xlab = "")
# ggplot:
ggplot() + geom_point(data = NULL, aes(x = grid, y = target_outer(grid)), color = "#325b96") + ylab("Density") + xlab ("") 
#plot(grid, a*target_sapply(grid), type = "l", lwd = 4, col = "#325b96")
#lines(grid, dnorm(grid, mean, sd), col = "#ee961d", lwd =3)

plot(grid, a*target_outer(grid), type = "l", lwd = 4, col = "#325b96", ylab = "Target density", xlab = "")
lines(grid, dnorm(grid, mean, sd), col = "#ee961d", lwd =3)

#what we need to minimize is
plot(grid, dnorm(grid, mean, sd), lwd =1, col = "#325b96" )
plot(grid, dnorm(grid, mean, sd)/target_outer(grid), lwd =1, col = "#325b96" )

#not good, try increasing sd by 9 percent
plot(grid, dnorm(grid, mean, sd*1.09)/target_sapply(grid), lwd =1, col = "#325b96" ,type = "l")

#the minimum occurs at that local minima in the middle

a <- optimize(function(z) dnorm(z, mean, sd*1.09)/target_sapply(z), c(0, 0.1))
a <- a$objective
sd
a
mean
#with the chosen a, we can replot the envelope and see what we get

plot(grid, a*target_sapply(grid), type = "l", lwd = 3,  col = "#325b96")
lines(grid, dnorm(grid, mean, sd*1.09),  col = "#ee961d", lwd =3)



#t distribution

plot(grid, at*target_sapply(grid), type = "l", lwd = 3,  col = "#325b96")
lines(grid, dt((grid - 0.24)/(0.0554*1.055), 61.83)/(0.0554*1.055),  col = "#ee961d", lwd =3)
help(dt)

plot(grid, dt((grid - 0.24)/(0.0554*1.055), 61.83)/(0.0554*1.055) / target_sapply(grid) )

at <- optimize(function(z) dt((z - 0.24)/(0.0554*1.055), 61.83)/(0.0554*1.055) / target_sapply(z), c(0, 0.4))$objective
#################### WITH THE ENVELOPE CHOSEN we can continue

#generic implementation

proposal_factory <- function(target_dens, proposal_dens, proposal_sim, a)
{
  list( tdens = target_dens, pdens = proposal_dens, sim = proposal_sim, alpha = a)
}

propn <- proposal_factory(target_outer, function(x) dnorm(x, 0.2388885, 0.061235426), function(x) rnorm(x, 0.2388885, 0.06123542), 0.9163125)
propt <- proposal_factory(target_outer,function(x) dt((x - 0.24)/(0.0554*1.055), 61.83)/(0.0554*1.055),function(x) rt(x, 61.83)*0.0554*1.055 + 0.24,0.95)


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
sim <- target_simulator(100000, propt, 1)
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
sim2 <- target_simulator_slow(1000, propn)
hist(sim)
hist(sim2, breaks = 20)

