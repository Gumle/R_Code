library(MASS)
library(stats)
library(microbenchmark)
library(ggplot2)
library(kedd)
library(caret)
setwd("~/Documents/Statistics - KU/Computational Statistics")
infrared <- read.table("infrared.txt", header = TRUE)
F12 <- infrared$F12
F12 <- log(F12)

EKernel <- function(x){
  ifelse(-1 <= x & x<= 1,(3/4)*(1-x^2), 0) 
  }
Kernel <- EKernel(F12)

# LOOCV
n <- length(F12)
h <- 1

LOOCV <- function(x,n,h){
  y <- 0
  for(i in 1:n){
    y[i] <- sum(1-(x[i]-x[-i])/h)^2
  }
  value <- sum(log(y))/h*(n-1)
  value
}
CV <- LOOCV(Kernel,n, h)
CV
LOOCV2 <- function (x, n, h){
  y <- 0
  for (i in 1:n){
    y[i] <- sum(3/4 * (1 - ((x[i] - x[-i])/h)^2) * (x[-i] - h <= x[i]) * (x[i] <= x[-i] + h)) * 1/(h * (n-1))
  }
  value <- sum(log(y))
  value
}
CV <- LOOCV(F12,n, h)
CV
VLOOCV <- function(x, n) {
  f <- function(h) LOOCV(x, n, h)
  f
}

VLOOCV1 <- VLOOCV(x = F12, n = 628)

VLOOCV2 <- Vectorize(VLOOCV1)

LOO <- VLOOCV2(seq(0.01, 4, 0.01))

which.max(LOO)
285*0.01
VLOOCV1(2.85)



### Kernel estimation:

kd <- function(x, h, m=512) {
  rg <- range(x)
  xx <-  seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m)
  const <-   const <- 3/(4*h*length(x))
  condition <-  lapply(xx, function(z) z - x)
  values <-  lapply(condition, function(z) sapply(z, function(x) ifelse(abs(x) <= h, 1 - (x/h)^2, 0  ) ) )
  y <-  sapply(values, function(z) sum(z)*const)
  
  list(x = xx, y = y)
}

x <- F12
m <- 512
h <- 1
rg <- range(x)
xx <-  seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m)
const <-   const <- 3/(4*h*length(x))
condition <- outer(xx,x, function(zz,z) zz-z)
#condition <-  lapply(xx, function(z) z - x)
sapply(condition, function(x) ifelse(abs(x) <= h, 1 - (x/h)^2, 0  ) )
outer(condition, function(x) ifelse(abs(x) <= h, 1 - (x/h)^2, 0  ) )
values <-  lapply(condition, function(z) sapply(z, function(x) ifelse(abs(x) <= h, 1 - (x/h)^2, 0  ) ) )
values <- outer(condition, function (z) outer(z, function(x) ifelse(abs(x) <= h, 1-(x/h)^2,0) ) )
y <-  sapply(values, function(z) sum(z)*const)


kd3 <- function (x, h, m = 512){
    rg <- range(x)
    ## xx is equivalent to grid points in 'density'
    xx <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m)
    y <- numeric(m) ## The evaluations, initialized as a vector of zeroes
    ## The actual computation is done using nested for-loops. The outer loop
    ## is over the grid points, and the inner loop is over the data points.
    const <- 3/(4* h*length(x))
    for (i in seq_along(xx))
      for (j in seq_along(x))
        if(abs(xx[i] - x[j]) <= h){
          y[i] <- y[i] + (1 - ((xx[i] - x[j])/h)^2)
        } else {
          y[i] <- y[i]
        }
    y <- y * const 
    list(x = xx, y = y)
}
EKernel <- function(x){
  ifelse(-1 <= x & x<= 1,x <- (3/4)*(1-x^2), 0) 
}

  kd3 <- function (x, h, m = 512) {
    rg <- range(x)
    ## xx is equivalent to grid points in 'density'
    xx <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m)
    #y <- numeric(m) ## The evaluations, initialized as a vector of zeroes
    ## The actual computation is done using nested for-loops. The outer loop
    ## is over the grid points, and the inner loop is over the data points.
    const <- h*length(x)
    y <- outer(xx, x, function(zz,z) EKernel(zz-z)/h)
    y <- rowMeans(y) * const 
    list(x = xx, y = y)
  }
  
  f_hat3 <- kd3(F12, h = 1)
  hist(F12, probability = TRUE)
  lines(f_hat3, col = "Blue")
  lines(f_hat_dens, col = "Red")  
  EKernel <- function(x){
    ifelse(-1 <= x & x<= 1,x <- (3/4)*(1-x^2), 0) 
  }