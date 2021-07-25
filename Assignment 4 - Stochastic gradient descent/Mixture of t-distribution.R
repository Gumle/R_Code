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
setwd("~/Documents/Statistics - KU/Computational Statistics/Assignments/Assignment 4 - Stochastic gradient descent")
#detectCores()
# Sample data
set.seed(5000)
data.sample <- function(n, par){
  mixt <- numeric(n)
  logi <- sample(c(TRUE,FALSE), size = n, prob = c(par[1],1-par[1]), replace = TRUE)
  mixt[logi] <- rt.scaled(n = sum(logi), df = par[6], mean = par[2], sd = sqrt(par[4]))
  mixt[!logi] <- rt.scaled(n = n-sum(logi), df = par[7], mean = par[3], sd = sqrt(par[5]))
  mixt
}

par <- c(1,2,3,1,1,4,1)
data <- data.sample(n = 10000, par = par)


# LogLikelihood

logLike <- function(par, x){
  -sum(log(par[1] * (gamma((par[6] + 1) /2) / ( sqrt(pi * par[6] * par[4]) * gamma(par[6]/2)) *
                       (1 + (((x - par[2])^2)/(par[6] * par[4]))) ^ (-(par[6] + 1) / 2)) +
             (1 - par[1]) * (gamma((par[7] + 1) /2) / ( sqrt(pi * par[7] * par[5]) * gamma(par[7]/2)) *
                               (1 + ((x - par[3])^2)/(par[7] * par[5])) ^ (-(par[7] + 1) / 2))))
}

par <- c(1,2,3,1,1,4,1)
par <- c(0.5,0.5,0.5,0.5,0.5,0.5,0.5)

x <- 1

logLike(par,x)

# Function factory so we can use the grad function.
logLike.factory <- function(x){
  function(par) logLike(par,x)
}

#c++ implementation
cppFunction(
  "#include <math.h>
  double logLikeCpp(NumericVector par, NumericVector x){
  int n = x.size();
  double sum = 0.0;
  for(int i = 0; i < n; i++){
sum = sum + log(par[0]*(tgamma((par[5]+1)/2)/(sqrt(M_PI*par[5]*par[3])*tgamma(par[5]/2)))*
  (pow((1+pow(x[i]-par[1],2)/(par[5]*par[3])),-(par[5]+1)/2))
  +(1-par[0])*(tgamma((par[6]+1)/2)/(sqrt(M_PI*par[6]*par[4])*tgamma(par[6]/2)))*
  (pow((1+pow(x[i]-par[2],2)/(par[6]*par[4])),-(par[6]+1)/2)));                
  }
  return -sum;
}")

logLikeCpp.factory <- function(x){
  function(par) logLikeCpp(par,x)
}

par <- c(0.5,0.5,0.5,0.5,0.5,0.5,0.5)
x <- c(1)
logLikeCpp(par,x)

source("analyticGradient.R")
source("loglikegrad.R")

par <- c(0.3,0.2,0.1,0.7,0.3,0.3,0.8)
x <- 1
loglikegrad(par,x)

SDG <- function(par, x, divlength, maxit, gamma){
  samplesize <- length(x)
  for(i in 1:maxit){
    index <- sample(1:samplesize, replace = FALSE)
    y <- x[index]
    for(j in 1:floor(samplesize/divlength)){
      f <- logLike.factory(x = y[(1 + ( j - 1)*divlength  ) : ( divlength * j ) ] )
      par <- par - gamma*gradient(f,par)
    }
  }
  par
}
par <- c(0.3,0.2,0.1,0.7,0.3,0.3,0.8)
SDG(par = par, x = data, divlength = 50, maxit = 100, gamma = 0.000001)


SGDFast <- function(par, H, difH, x, divlength, gamma, epsilon = 10^(-6), cb = NULL){
  samplesize <- length(x)
  j_end <- floor(samplesize/divlength) 
  repeat{
    value <- H(par, x)
    index <- sample(1: samplesize, replace = FALSE)
    y <- x[index]
    for (j in 1: j_end){
      curgrad <- difH(par, y[(1 + (j - 1) * divlength ): (divlength * j )])
      par[2:7] <- par[2:7] - gamma * 3 * curgrad[2:7]
      par[1] <- par[1] - (gamma/10) * curgrad[1]
    }
    if(!is.null(cb)) cb()
    criteria <- value - H(par, x) - epsilon * (H(par,x) + epsilon)
    if(criteria <= 0) {
      break
    }
  }
  par
}

sgd_tracer <- tracer(c("criteria"), N = 1)

sgdFast(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), H = logLike, difH = loglikegrad, x = data, divlength = 100, 
     gamma = 0.000005)


par <- c(0.5,0.5,0.5,0.5,0.5,0.5,0.5)

value <- SDG2(par = par, x = data, divlength = 50, maxit = 100, gamma = 0.0001)

system.time(SDG(par = par, x = data, divlength = 50, maxit = 100, gamma = 0.000001))
system.time(SDG2(par = par, x = data, divlength = 50, maxit = 100, gamma = 0.0001))

GD <- function(par,x, H, d = 0.8,  c = 0.1,  gamma0 = 0.01,   epsilon = 1e-4,  cb = NULL) {
  H <- H(x)
  repeat {
    value <- H(par)
    curgrad <- grad(H,par)
    h_prime <- sum(curgrad^2)
    if(!is.null(cb)) cb()
    ## Convergence criterion based on gradient norm
    if(h_prime <= epsilon) break
    gamma <- gamma0
    ## First proposed descent step
    par1 <- par - gamma * curgrad
    ## Backtracking while descent is insufficient
    while(H(par1) > value - c * gamma * h_prime) {
      gamma <- d * gamma
      par1 <- par - gamma * curgrad
    }
    par <- par1
  }
  par
}

#GD(par = par,x = data, H = logLike.factory)
par <- c(0.3,2,3,1,1,4,1)
data <- data.sample(n = 1000, par = par)
# Compare the two descent methods.
gradDesc <- GD(par = par, x = data, H = logLike.factory, gamma0 = 0.000001)
stochastic <- SDG(par =par, x = data, divlength = 50, maxit = 100, gamma = 0.000001)
system.time(GD(par = par, x = data, H = logLike.factory, gamma0 = 0.000001))
system.time( SDG(par =par, x = data, divlength = 50, maxit = 100, gamma = 0.000001))
abs(gradDesc-stochastic)

# Paralell package --------------------#
set.seed(5000)
par <- c(0.5,0.5,0.5,0.5,0.5,0.5,0.5)
par0 <- c(1,2,3,1,1,4,1)
RNGkind("L'Ecuyer-CMRG")  
set.seed(10)

SGD_parhat <- mclapply(
  1:12,function(m) SDG(par = par, x = data.sample(n = 10000/12, par = par0), divlength = 50, maxit = 100, gamma = 0.0001),
  mc.cores = 12
)
SGD_parhat

GD_parhat <- mclapply(
  1:12,function(m)  GD(par = par, x = data.sample(n = 10000/12, par = par0), H = logLike.factory, gamma0 = 0.0001),
  mc.cores = 12
)
GD_parhat



#---- SGD & GD work, however, the GD is very slow because and so is the logLike function.












# EM - Algorithm
Estep <- function(par, x){
  g1 <- par[1]*dt.scaled(x = x, df = par[6], mean = par[2], sd = sqrt(par[4]))
  g2 <- (1-par[1])*dt.scaled(x = x, df = par[7], mean = par[3], sd = sqrt(par[5]))
  g1/(g1+g2)
}

par <- c(0.4, 1,2, 3, 0.5, 0.7, 2)
x <- data
delta0 <- E(par,x)
delta0

MStep.factory <- function(par, x, delta){
 d1 <- delta*log(par[1]*dt.scaled(x = x, df = par[6], mean = par[2], sd = sqrt(par[4])))
 d2 <- (1-delta)*log((1-par[1])*dt.scaled(x = x, df = par[7], mean = par[3], sd = sqrt(par[5])))
 -sum(d1+d2)
}
MStep.factory(c(0.5,1,1,1,1,1,1), x, delta0)

Mstep <- function(par, x, delta, d = 0.8,  c = 0.1,  gamma0 = 0.00001,   epsilon = 1e-4,  cb = NULL ){
  H <- function(par) MStep.factory(par,x,delta)
  value <- H(par)
  curgrad <- grad(H,par)
  h_prime <- sum(curgrad^2)
  if(!is.null(cb)) cb()
  ## Convergence criterion based on gradient norm
  if(h_prime <= epsilon) break
  gamma <- gamma0
  ## First proposed descent step
  par1 <- par - gamma * curgrad
  ## Backtracking while descent is insufficient
  while(H(par1) > value - c * gamma * h_prime) {
    gamma <- d * gamma
    par1 <- par - gamma * curgrad
  }
  par <- par1
  par
}

Mstep(par,data, delta = 1)

EM <- function(par, y, epsilon = 10^(-6), trace = NULL) {
  repeat{
    par0 <- par
    par <- Mstep(par = par , x = y, delta = Estep(par, y))
    if(!is.null(trace)) trace()
    if(sum((par - par0)^2) <= epsilon * (sum(par^2) + epsilon))
      break
  } 
  par  ## Remember to return the parameter estimate
}
source("Debugging_and_tracing.R")
EM_tracer <- tracer("par",N = 10)
EM(par = c(0.4, 1,2, 3, 0.5, 0.7, 2), y = data, trace = EM_tracer$trace)
EM_tracer[20]

