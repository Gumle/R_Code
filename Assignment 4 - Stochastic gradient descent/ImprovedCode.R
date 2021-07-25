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
setwd("~/Documents/Statistics - KU/Computational Statistics/Assignments/Assignment 4 - Stochastic gradient descent")
source("Debugging_and_tracing.R")


set.seed(5000)

datasampler <- function(n, par){
  mixt <- numeric(n)
  logi <- sample(c(TRUE,FALSE), size = n, prob = c(par[1], 1 - par[1]), replace = TRUE)
  mixt[logi] <- rt.scaled(n = sum(logi), df = par[6], mean = par[2], sd = sqrt(par[4]))  
  mixt[!logi] <- rt.scaled(n = n - sum(logi), df = par[7], mean = par[3], sd = sqrt(par[5]))
  mixt
}

library(microbenchmark)

data <- datasampler(n = 10000, par = c(0.3, 10, 7, 3, 0.5, 2, 4))

# Object oriented programming

ingredient <- function(data, par, divlength, gamma){
  structure(list("data" = data, "par" = par, "divlength" = divlength, "gamma" = gamma) , class = "sgdob")
}

objir <- ingredient(data = data, par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), divlength = 100, gamma = 0.0005)

mix2t <- function(obj, ...){
  UseMethod("mix2t")
}

mix2t.sgdob <- function(obj){
  values<- sgd2(obj$par, logLike, analyticGrad, obj$data, obj$divlength, 
                obj$gamma, epsilon = 10^(-6), cb = NULL)
  values
}

mix2t(objir)


logLike <- function(par, x){
  -sum(log(par[1] * (gamma((par[6] + 1) /2) / ( sqrt(pi * par[6] * par[4]) * gamma(par[6]/2)) *
                       (1 + (((x - par[2])^2)/(par[6] * par[4]))) ^ (-(par[6] + 1) / 2)) +
             (1 - par[1]) * (gamma((par[7] + 1) /2) / ( sqrt(pi * par[7] * par[5]) * gamma(par[7]/2)) *
                               (1 + ((x - par[3])^2)/(par[7] * par[5])) ^ (-(par[7] + 1) / 2))))
}

logLike.factory <- function(x){
  function(par) logLike(par, x)
}

source("loglikegrad.R")

# Need to test if it is correctly computed

library(numDeriv)

data <- datasampler(n = 10000, par = c(0.3, 10, 7, 3, 0.5, 2, 4))

lf <- logLike.factory(data)

grad(lf, c(0.3, 1, 2, 3, 4, 5, 6))

analyticGrad(par = c(0.3, 1, 2, 3, 4, 5, 6), data)


# now it is tested if it is faster.

library(microbenchmark)

p <- microbenchmark(times = 200, analyticGrad(par = c(0.5, 1, 1, 1, 2, 2, 2), data), 
               grad(logLike.factory(data), c(0.5, 1, 1, 1, 2, 2, 2)))

summary(p)

autoplot(p) + geom_jitter(aes(color = expr), alpha = 0.4) + aes(fill = I("gray"))+
  scale_color_manual(values = c("#354e75", "#884091", "black")) + theme_bw() +
  theme(legend.position = "none",text = element_text(size=10))

#Now stochastic gradient descent is implemented. 

sgd <- function(par, x, divlength, maxiter, gamma){
  samplesize <- length(x)
  for (i in 1: maxiter){
    index <- sample(1: samplesize, replace = FALSE)
    y <- x[index]
    for (j in 1: floor(samplesize/divlength)){
      f <- logLike.factory(x = y[(1 + (j - 1) * divlength ): (divlength * j )])
      par <- par - gamma * grad(f, par)
    }
  }
  par
}

sgd(par = par, x = data, divlength = 50, maxit = 100, gamma = 0.000001)


sgd2 <- function(par, H, difH, x, divlength, gamma, epsilon = 10^(-6), cb = NULL){
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

sgd2(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), H = logLike, difH = analyticGrad, x = data, divlength = 100, 
     gamma = 0.000005)

data2 <- datasampler(n = 1000000, par = c(0.3, 10, 7, 3, 0.5, 2, 4))


sgd2(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), H = logLike, difH = analyticGrad, x = data2, divlength = 100, 
     gamma = 0.00005, cb = sgd_tracer$trace)
summary(sgd_tracer)

#Continous problem with bad convergence of the nu parameters. 

sgd3 <- function(par, H, difH, x, divlength, gamma, epsilon = 10^(-6), cb = NULL){
  samplesize <- length(x)
  j_end <- floor(samplesize/divlength) 
  repeat{
    value <- H(par, x)
    index <- sample(1: samplesize, replace = FALSE)
    y <- x[index]
    for (j in 1: j_end){
      subtract <- difH(par, y[(1 + (j - 1) * divlength ): (divlength * j )])
      par[6:7] <- par[6:7] - gamma * 400 * subtract[6:7]
      par[2:5] <- par[2:5] - gamma * 3 * subtract[2:5]
      par[1] <- par[1] - (gamma/10) * subtract[1]
    }
    print(par)
    if(!is.null(cb)) cb()
    criteria <- value - H(par, x) - epsilon * (H(par,x) + epsilon)
    if(criteria <= 0) {
      break
    }
  }
  par
}

sgd3(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), H = logLike, difH = analyticGrad, x = data, divlength = 100, 
     gamma = 0.0005)


sgd4 <- function(par, H, difH, x, divlength, gamma0, gamma1, epsilon = 10^(-6), cb = NULL){
  samplesize <- length(x)
  j_end <- floor(samplesize/divlength) 
  repeat{
    value <- H(par, x)
    index <- sample(1: samplesize, replace = FALSE)
    y <- x[index]
    for (j in 1: j_end){
      subtract <- difH(par, y[(1 + (j - 1) * divlength ): (divlength * j )])
      par <- par - sign(subtract) * gamma0 * sum(abs(subtract))
      par <- par - gamma1 * subtract
    }
    print(par)
    if(!is.null(cb)) cb()
    criteria <- value - H(par, x) - epsilon * (H(par,x) + epsilon)
    if(criteria <= 0) {
      break
    }
  }
  par
}



#Trying to find the optimal learning rate

sgd_tracer <- tracer(c("criteria"), N = 5)

sgd2(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), H = logLike, difH = analyticGrad, x = data, divlength = 100, 
     gamma = 0.0005, cb = sgd_tracer$trace)

times_sg <- microbenchmark(sgd2(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), H = logLike, difH = analyticGrad, x = data[1:300], divlength = 100, 
                    gamma = 0.0005),
               sgd2(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), H = logLike, difH = analyticGrad, x = data[1:900], divlength = 100, 
                    gamma = 0.0005),
               sgd2(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), H = logLike, difH = analyticGrad, x = data[1:2700], divlength = 100, 
                    gamma = 0.0005),
               sgd2(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), H = logLike, difH = analyticGrad, x = data[1:8100], divlength = 100, 
                    gamma = 0.0005), times = 2)

summary(times_sg)
times_sg

sgd_tracer2 <- tracer(c("criteria"), N = 1)

data2 <- datasampler(n = 1000000, par = c(0.3, 10, 7, 3, 0.5, 2, 4))

sgd2(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), H = logLike, difH = analyticGrad, x = data2, divlength = 100, 
gamma = 0.00005)

summary(sgd_tracer2)

## Now for nagd. The good thing can be that it pushes in a specific direction. 

nagd <- function(par, H, difH, alpha, x, gamma, epsilon = 10 ^(-6), cb = NULL){
  par_x <- par
  pary_0 <- par
  repeat{
    value <- H(par_x, x)
    pary_1 <-par_x  - gamma * difH(par_x, x) 
    par_x <- pary_1 + alpha * (pary_1 - pary_0)
    pary_0 <- pary_1
    if(!is.null(cb)) cb()
    criteria <- value - H(par_x, x) - epsilon * (H(par_x,x) + epsilon) 
    if (criteria <= 0) {
      break
    }
  }
  par_x
}

nagd_tracer <- tracer(c("criteria"), N = 5)

#0.989 for learning rate. 
nagd(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), H = logLike, difH = analyticGrad, alpha = 0.9, x = data, 
      gamma = 0.000000989, cb = nagd_tracer$trace )

timesnagd <- microbenchmark(nagd(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), H = logLike, difH = analyticGrad, alpha = 0.9, x = data[1:300], 
                                 gamma = 0.000000989),
                            nagd(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), H = logLike, difH = analyticGrad, alpha = 0.9, x = data[1:900], 
                                 gamma = 0.000000989 ),
                            nagd(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), H = logLike, difH = analyticGrad, alpha = 0.9, x = data[1:2700], 
                                 gamma = 0.000000989 ),
                            nagd(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), H = logLike, difH = analyticGrad, alpha = 0.9, x = data[1:8100], 
                                 gamma = 0.000000989 ), times = 2)

timesnagd

#now for nasgd


nasgd <- function(par, H, difH, alpha, x, gamma, epsilon = 10^(-6), divlength = 50, cb = NULL){
  par_x <- par
  pary_0 <- par
  samplesize <- length(x)
  j_end <- floor(samplesize/divlength)
  repeat{
    value <- H(par_x, x)
    index <- sample(1: samplesize, replace = FALSE)
    y <- x[index]
    for (j in 1:j_end){
      pary_1 <-par_x  - gamma * difH(par_x, y[(1 + (j - 1) * divlength ): (divlength * j )]) 
      par_x <- pary_1 + alpha * (pary_1 - pary_0)
      pary_0 <- pary_1
    }
    if(!is.null(cb)) cb()
    criteria <- value - H(par_x, x) - epsilon * (H(par_x,x) + epsilon) 
   if (criteria <= 0) {
     break
   }
  }
  par_x
}

nasgd_tracer <- tracer(c("criteria"), N = 5)

nasgd(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), H = logLike, difH = analyticGrad, alpha = 0.9 
      , x = data, gamma = 0.0000025, cb = nasgd_tracer$trace)

timesnasgd <- microbenchmark(nasgd(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), H = logLike, difH = analyticGrad, alpha = 0.9 
               , x = data[1:300], gamma = 0.0000025), 
               nasgd(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), H = logLike, difH = analyticGrad, alpha = 0.9 
                                                           , x = data[1:900], gamma = 0.0000025),
               nasgd(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), H = logLike, difH = analyticGrad, alpha = 0.9 
                     , x = data[1:2700], gamma = 0.0000025),
               nasgd(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), H = logLike, difH = analyticGrad, alpha = 0.9 
                     , x = data[1:8100], gamma = 0.0000025), times = 2)

  
### Now for gd

GD <- function(par,
               x, 
               H,
               derivH, 
               d = 0.8, 
               c = 0.1, 
               gamma0 = 0.01, 
               epsilon = 1e-6, 
               cb = NULL) {
  repeat {
    value <- H(par, x)
    curgrad <- derivH(par, x)
    h_prime <- sum(curgrad^2)
    if(!is.null(cb)) cb()
    ## Convergence criterion based on gradient norm
    #if(h_prime <= epsilon) break
    gamma <- gamma0 
    ## First proposed descent step
    par1 <- par - gamma * curgrad
    ## Backtracking while descent is insufficient
    while(H(par1, x) > value - c * gamma * h_prime) {
      gamma <- d * gamma
      par1 <- par - gamma * curgrad
    }
    criteria <- value - H(par1, x) - epsilon * (H(par1,x) + epsilon)
    if (criteria <= 0) {
      break
    }
    par <- par1
  }
  par
}

gd_tracer <- tracer(c("criteria"), N = 5)

GD(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), x = data, H = logLike, derivH = analyticGrad, gamma0 = 0.000065, 
   cb = gd_tracer$trace)


timesgd <- microbenchmark(GD(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), x = data[1:300], H = logLike, derivH = analyticGrad,
                             gamma0 = 0.000065),
                          GD(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), x = data[1:900], H = logLike, derivH = analyticGrad,
                             gamma0 = 0.000065 ),
                          GD(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), x = data[1:2700], H = logLike, derivH = analyticGrad, 
                             gamma0 = 0.000065 ),
                          GD(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), x = data[1:8100], H = logLike, derivH = analyticGrad, 
                             gamma0 = 0.000065 ), times = 3)

timesgd
###############ggplot2 part#############################

library(ggplot2)

nasgd_tracer_data <- summary(nasgd_tracer)

nagd_tracer_data <- summary(nagd_tracer)

sgd_tracer_data <- summary(sgd_tracer)

gd_tracer_data <- summary(gd_tracer)


ggplot(data = nasgd_tracer_data, aes(.time, criteria)) + geom_point(col = "red") + 
  geom_point(data = nagd_tracer_data, aes(.time, criteria)) + geom_point(data = gd_tracer_data,
    aes(.time, criteria), col = "blue") + 
  geom_point(data = sgd_tracer_data, aes(.time, criteria, col = "yellow")) + scale_y_log10()



linesdatasgd <- data.frame("points" = c(300, 900, 2700, 8100), 
                          "seconds" = c(1.1487591, 0.6430827, 0.7299827,1.0825354))

linesdatanagd <- data.frame("points" = c(300, 900, 2700, 8100), 
                        "seconds" = c(48.26398, 46.42349, 57.49411, 79.74165))

linesdatanasgd <- data.frame("points" = c(300, 900, 2700, 8100),
                             "seconds" = c(29.02544, 24.75910, 29.62852, 37.92659))

linesdatagd <- data.frame("points" = c(300, 900, 2700, 8100), 
                          "seconds" =c(14.43648, 14.57062, 17.03292, 23.16934))

linesdataem <- data.frame("points" = c(300, 900, 2700, 8100),
                          "seconds" = 18.62617, 16.88839, 22.72270, 27.33615)

linesdataemgau <- data.frame("points" = c(300, 900, 2700, 8100),
                             "seconds" = 2.829419, 3.324136, 4.581812, 4.896264)

ggplot(data = linesdatasgd, aes(points, seconds)) + geom_line(col = "red") + 
  geom_line(data = linesdatanagd, aes(points, seconds), col = "blue") + 
  geom_line(data = linesdatanasgd, aes(points, seconds), col = "green") + 
  geom_line(data = linesdatagd, aes(points, seconds), col = "black")  +
  geom_line(data = linesdataem, aes(points, seconds), col = "yellow") + 
  geom_line(data= linesdataemgau, aes(points, seconds), col = "purple" ) + 
  scale_x_log10()

#############################

Estep <- function(par, x){
  g1 <- par[1] * (gamma((par[6] + 1) /2) / ( sqrt(pi * par[6] * par[4]) * gamma(par[6]/2)) *
    (1 + (((x - par[2])^2)/(par[6] * par[4]))) ^ (-(par[6] + 1) / 2))
  g2 <- (1 - par[1]) * (gamma((par[7] + 1) /2) / ( sqrt(pi * par[7] * par[5]) * gamma(par[7]/2)) *
                          (1 + ((x - par[3])^2)/(par[7] * par[5])) ^ (-(par[7] + 1) / 2))
  g1 / (g1 + g2)
}

kf <- Estep(par = c(0.3, 1, 2, 3, 4, 5, 6), x = c(1,2))

#Mstepfunction <- function(par, x, delta){
# -sum(delta * log(par[1] * dt.scaled(x, mean = par[2], df = par[6], sd = sqrt(par[4]))) + 
#        (1- delta) * log((1 - par[1]) * dt.scaled(x, mean = par[3], df = par[7], sd = sqrt(par[5]))))
#}

Mstepfunction <- function(par, x, delta){
  -sum(delta * log(par[1] * (gamma((par[6] + 1) /2) / ( sqrt(pi * par[6] * par[4]) * gamma(par[6]/2)) *
                               (1 + (((x - par[2])^2)/(par[6] * par[4]))) ^ (-(par[6] + 1) / 2)))) -
         sum((1- delta) * log((1 - par[1]) * (gamma((par[7] + 1) /2) / ( sqrt(pi * par[7] * par[5]) * gamma(par[7]/2)) *
                                            (1 + ((x - par[3])^2)/(par[7] * par[5])) ^ (-(par[7] + 1) / 2))))
}


Mstepfunction(par = c(0.5, 1, 1, 1, 1, 1, 1), x = c(1, 2, 3, 4), delta = kf)

library(numDeriv)

#Need to find an expression for the gradient below

nu1gradfunc <- function(x, par){
  (((1 + (((x - par[2])^2)/ (par[6] * par[4]))) ^ ( - (par[6] + 1) / 2)) * 
     (( digamma( (par[6] + 1) / 2)  * gamma((par[6] + 1)/2) * sqrt(pi * par[6] * par[4]) * gamma(par[6]/2) * 1/2) -
        gamma((par[6] + 1)/2) * (sqrt(pi * par[6] * par[4]) * gamma(par[6]/2) * digamma(par[6]/2) * 1/2 + 
                                   gamma(par[6]/2) * 1/2 * ((pi * par[6] * par[4]) ^ (-1/2)) * pi * par[4])) / (
                                     (pi * par[6] * par[4] * (gamma(par[6]/2) ^2)))) + 
    
    ((gamma((par[6] + 1) / 2) / (sqrt(pi * par[6] * par[4]) * gamma(par[6] / 2))) * 
       ((1 + (((x - par[2]) ^ 2) / (par[4] * par[6]))) ^ (- (par[6] + 1)/2)) * 
       
       ((-1/2 * log(1 + (((x - par[2])^2)  /(par[4] * par[6])))) + (   ((par[6] + 1) / 2) * 
                                                                         
                                                                         
      (  (((x - par[2])^2) / (par[4] * (par[6]^2))) /  (1 + (((x - par[2])^2) / (par[4] * par[6])) )))))
  
}

nu2gradfunc <- function(x, par){
  (((1 + (((x - par[3])^2)/ (par[7] * par[5]))) ^ ( - (par[7] + 1) / 2)) * 
     (( digamma( (par[7] + 1) / 2)  * gamma((par[7] + 1)/2) * sqrt(pi * par[7] * par[5]) * gamma(par[7]/2) * 1/2) -
        gamma((par[7] + 1)/2) * (sqrt(pi * par[7] * par[5]) * gamma(par[7]/2) * digamma(par[7]/2) * 1/2 + 
                                   gamma(par[7]/2) * 1/2 * ((pi * par[7] * par[5]) ^ (-1/2)) * pi * par[5])) / (
                                     (pi * par[7] * par[5] * (gamma(par[7]/2) ^2)))) + 
    
    ((gamma((par[7] + 1) / 2) / (sqrt(pi * par[7] * par[5]) * gamma(par[7] / 2))) * 
       ((1 + (((x - par[3]) ^ 2) / (par[5] * par[7]))) ^ (- (par[7] + 1)/2)) * 
       
       ((-1/2 * log(1 + (((x - par[3])^2)  /(par[5] * par[7])))) + (   ((par[7] + 1) / 2) * 
                                                                         
                                                                         
        (  (((x - par[3])^2) / (par[5] * (par[7]^2))) /  (1 + (((x - par[3])^2) / (par[5] * par[7])) )))))
}

derivMstep <- function(par, x, delta){
  n <- length(delta)
  
  fvalue1 <- (gamma((par[6] + 1)/2) / (sqrt(pi * par[6] * par[4]) * gamma(par[6]/2))) *
    (1 + ((x - par[2])^2)/ (par[6] * par[4])) ^ (- (par[6] + 1) /2)
  fvalue2 <- (gamma((par[7] + 1)/2) / (sqrt(pi * par[7] * par[5]) * gamma(par[7]/2))) *
    (1 + ((x - par[3])^2)/ (par[7] * par[5])) ^ (- (par[7] + 1) /2)
  
  pgrad <- -sum(delta / par[1]) +  sum((1 - delta)/(1 - par[1]))
  
  mu1grad <- sum((- delta * (gamma((par[6] + 1)/2) /(sqrt(pi * par[6] * par[4]) * gamma(par[6] / 2))) *
    (((par[6] + 1) * (((((x - par[2])^2)/(par[4] * par[6])) + 1) ^ ((- (par[6] + 3)) / 2)) * (x - par[2])) / 
       (par[4] * par[6])))/fvalue1)
  
  mu2grad <- - sum((1 - delta) * (gamma((par[7] + 1)/2) /(sqrt(pi * par[7] * par[5]) * gamma(par[7] / 2))) *
    (((par[7] + 1) * (((((x - par[3])^2)/(par[5] * par[7])) + 1) ^ ((- (par[7] + 3)) / 2)) * (x - par[3])) / 
       (par[5] * par[7])) /fvalue2)
  
  sigma1grad <- -sum((-delta * ((gamma((par[6]+1)/2) * (par[6]^2) * (par[4] - (x - par[2])^2)) / 
    (2 * sqrt(pi) * gamma(par[6]/2) * ((par[6] * par[4])^(3/2)) * (par[6] * par[4] + (x - par[2])^2) * 
       (1 + (((x - par[2])^2)/(par[4] * par[6] ))) ^((par[6] + 1)/2)))) /fvalue1)
  
  sigma2grad <- sum(- ((1 - delta) * (-(gamma((par[7]+1)/2) * (par[7]^2) * (par[5] - (x - par[3])^2)) / 
    (2 * sqrt(pi) * gamma(par[7]/2) * ((par[7] * par[5])^(3/2)) * (par[7] * par[5] + (x - par[3])^2) * 
       (1 + (((x - par[3])^2)/(par[5] * par[7] ))) ^((par[7] + 1)/2)))/ fvalue2))
  
  nu1grad <-  -sum((delta * nu1gradfunc(x, par))/fvalue1)
  
  nu2grad <- - sum((1-delta) * nu2gradfunc(x, par)/fvalue2)
  
  c(pgrad, mu1grad, mu2grad, sigma1grad, sigma2grad, nu1grad, nu2grad)
}

Mstep_tracer <- tracer(c("curgrad", "curgrad2"), N = 1)


library(numDeriv)

Mstep <- function(H, difH, par, x, delta, gamma0 = 0.00006, d = 0.8, c = 0.1, cb = NULL){
  if(!is.null(cb)) cb()
  value <- H(par, x, delta)
  curgrad <- difH(par, x, delta)
  h_prime <- sum(curgrad^2)
  ## Convergence criterion based on gradient norm
  gamma <- gamma0
  ## First proposed descent step
  par1 <- par - gamma * curgrad
  ## Backtracking while descent is insufficient
  while(H(par1, x, delta) > value - c * gamma * h_prime) {
    gamma <- d * gamma
    par1 <- par - gamma * curgrad
  }
  par <- par1
  par
}

Mstep(Mstepfunction, derivMstep, par = c(0.3, 1, 2, 3, 4, 5, 6), x = c(1, 2), delta = kf)

emalgoritm <- function(loglike, H, difH, par, x, gamma, epsilon = 10^(-4)){
  repeat{
   value <- loglike(par, x)
   delta <- Estep(par = par, x = x)
   par <- Mstep(H, difH, par = par, x = x, delta = delta, gamma0 = gamma)
   criteria <- value - loglike(par, x) - epsilon * (loglike(par,x) + epsilon)
   if (criteria <= 0) {
     break
   }
   par
   }
}
emalgoritm(logLike, Mstepfunction, derivMstep, par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9),  epsilon = 10 ^ (-6), x = data, gamma = 0.00006)

microbenchmark(emalgoritm(logLike, Mstepfunction, derivMstep, par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), x = data[1:300], gamma = 0.00006, epsilon = 10 ^ (-6)),
               emalgoritm(logLike, Mstepfunction, derivMstep, par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), x = data[1:900], gamma = 0.00006, epsilon = 10 ^ (-6)),
               emalgoritm(logLike, Mstepfunction, derivMstep, par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), x = data[1:2700], gamma = 0.00006, epsilon = 10 ^ (-6)),
               emalgoritm(logLike, Mstepfunction, derivMstep, par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), x = data[1:8100], gamma = 0.00006, epsilon = 10 ^ (-6)),
               times = 2)

emalgoritm(logLike, Mstepfunction, derivMstep, par = c(0.3, 10, 7, 3, 0.5, 2, 4), x = data, epsilon = 10 ^ (-6))

emalgoritm(logLike, Mstepfunction, derivMstep, par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), x = data, epsilon = 10 ^ (-6))



# Gaussian mixture model
# Here par[1] is p, par[2] is mu1, par[3] is mu2, par[4] is sigma1squared, par[5] is sigma2squared

EstepGau <- function(par, x){
  g1 <- par[1] * 1/(sqrt(2 * pi * par[4])) * exp(-((x - par[2])^2)/(2 * par[4])) + 10^(-6)
  g2 <- (1 - par[1]) * 1/(sqrt(2 * pi * par[5])) * exp(-((x - par[3])^2)/(2 * par[5])) + 10^(-6)
  g1/(g1 + g2)
}


Mstepfunctiongau <- function(par, x, delta){
  -sum(delta * log(par[1] * 1/(sqrt(2 * pi * par[4])) * exp(-((x - par[2])^2)/(2 * par[4]))) + 
        (1 - delta) * log((1 - par[1]) * 1/(sqrt(2 * pi * par[5])) * exp(-((x - par[3])^2)/(2 * par[5]))) + 10^(-6))
}


kf <- EstepGau(par = c(0.3, 10, 7, 3, 0.5), x = data)

#fits better if nu parameter is quite high. 

derivMstepgau <- function(par, x, delta){
  pgrad <-  -sum(delta)/par[1] + sum(1 - delta) / (1 - par[1])
  mu1grad <- -sum(delta * ((x - par[2])/par[4]))
  mu2grad <- -sum((1 - delta) * ((x - par[3])/par[5]))
  sigma1grad <- sum(delta * ((par[4] - (x - par[2])^2)/ (2 * par[4]^2)))
  sigma2grad <- sum((1 - delta) * ((par[5] - (x - par[3])^2)/ (2 * par[5]^2)))
  c(pgrad, mu1grad, mu2grad, sigma1grad, sigma2grad)
}

derivMstepgau(par = c(0.3, 10, 7, 3, 0.5), x = data, delta = kf)

MstepGau <- function(H, difH, par, x, delta, gamma0 = 0.0001, d = 0.8, c = 0.1){
  value <- H(par, x, delta)
  curgrad <- difH(par, x, delta)
  h_prime <- sum(curgrad^2)
  #if(!is.null(cb)) cb()
  ## Convergence criterion based on gradient norm
  gamma <- gamma0
  ## First proposed descent step
  par1 <- par - gamma * curgrad
  ## Backtracking while descent is insufficient
  while(H(par1, x, delta) > value - c * gamma * h_prime) {
    gamma <- d * gamma
    par1 <- par - gamma * curgrad
  }
  par <- par1
  par
}

MstepGau(H = Mstepfunctiongau, difH = derivMstepgau, delta = kf, par = c(0.3, 1, 2, 3, 4), x = c(1, 2, 3, 4))

emgautracer <- tracer(c("criteria"), N = 1)

emalgoritmgau <- function(loglike, H, difH, par, x, epsilon = 10^(-4), gamma0, cb = NULL){
  repeat{
    value <- loglike(par, x)
    delta <- EstepGau(par = par, x = x)
    par <- MstepGau(H, difH, par = par, x = x, delta = delta, gamma0 = gamma0)
    #if(!is.null(cb)) cb()
    criteria <- value - loglike(par, x) - epsilon * (loglike(par,x) + epsilon)
    if (criteria <= 0) {
      break
    }
  }
  par
}

minusloglikegau <- function(par, x){
  -sum(log(par[1] * (1/(sqrt(2 * pi * par[4])) * exp(-((x - par[2])^2)/(2 * par[4]))) +
             (1 - par[1]) * ( 1/(sqrt(2 * pi * par[5])) * exp(-((x - par[3])^2)/(2 * par[5]))) + 10^(-6)))
}


# gets values that are very small. 

minusloglikegau(par = c(0.3, 10, 7, 3, 0.5), x = data)

datanu <- datasampler(n = 10000, par = c(0.3, 10, 7, 3, 0.5, 100, 200))

microbenchmark(emalgoritmgau(minusloglikegau, H = Mstepfunctiongau, difH = derivMstepgau, par = c(0.2, 12, 6, 4, 0.7), x = datanu[1:300], 
                             epsilon = 10 ^ (-6), gamma0 = 0.00005),
               emalgoritmgau(minusloglikegau, H = Mstepfunctiongau, difH = derivMstepgau, par = c(0.2, 12, 6, 4, 0.7), x = datanu[1:900], 
                             epsilon = 10 ^ (-6), gamma0 = 0.00005),
               emalgoritmgau(minusloglikegau, H = Mstepfunctiongau, difH = derivMstepgau, par = c(0.2, 12, 6, 4, 0.7), x = datanu[1:2700], 
                             epsilon = 10 ^ (-6), gamma0 = 0.00005),
               emalgoritmgau(minusloglikegau, H = Mstepfunctiongau, difH = derivMstepgau, par = c(0.2, 12, 6, 4, 0.7), x = datanu[1:8100], 
                             epsilon = 10 ^ (-6), gamma0 = 0.00005), times = 2
  )

emalgoritmgau(minusloglikegau, H = Mstepfunctiongau, difH = derivMstepgau, par = c(0.2, 12, 6, 4, 0.7), x = datanu, 
              epsilon = 10 ^ (-6), gamma0 = 0.00005)

emalgoritmgau(minusloglikegau, H = Mstepfunctiongau, difH = derivMstepgau, par = c(0.4, 12, 5, 4, 1), x = data, 
              epsilon = 10 ^ (-6), gamma0 = 0.000001)

data2 <- datasampler(n = 1000000, par = c(0.3, 10, 7, 3, 0.5, 100, 200))


#Quite good if one knows there is a high shape parameter





emalgoritmgau(minusloglikegau, H = Mstepfunctiongau, difH = derivMstepgau, par = c(0.4, 12, 5, 4, 1), x = data2, 
              epsilon = 10 ^ (-6), gamma0 = 0.0000007, cb = emgautracer$trace)
emgaudata <- summary(emgautracer)

microbenchmark(emalgoritmgau(minusloglikegau, H = Mstepfunctiongau, difH = derivMstepgau, par = c(0.4, 12, 5, 4, 1), x = data2[1:10000], 
           epsilon = 10 ^ (-6), gamma0 = 0.00007), times = 2)
#Only better if the nu's are quite high because then a more normal distributional shape. 

ggplot(data = emgaudata, aes(.time, criteria)) + scale_y_log10() +  geom_point()


# Try to see if location parameter can be calculated from the Gaussian, scale parameter from sigma squared.
  


###################################################


datalargesgd <- summary(sgd_tracer2)

datalargeemgau <- read.table("emgautracer.txt")

str(datalargeemgau)

ggplot(data = datalargesgd, aes(.time, criteria)) + geom_point(col = "red") + 
  geom_point(data = datalargeemgau, aes(.time, criteria), col = "blue") + scale_y_log10()


sgd2(par = c(0.2, 12, 6, 4, 0.7, 2.4, 4.9), H = logLike, difH = analyticGrad, x = data2, divlength = 100, 
     gamma = 0.00005, cb = sgd_tracer2$trace)

sgd2 <- function(par, H, difH, x, divlength, gamma, epsilon = 10^(-6), cb = NULL){
  samplesize <- length(x)
  j_end <- floor(samplesize/divlength) 
  repeat{
    value <- H(par, x)
    index <- sample(1: samplesize, replace = FALSE)
    y <- x[index]
    for (j in 1: j_end){
      par[2:7] <- par[2:7] - gamma * 3 * difH(par, y[(1 + (j - 1) * divlength ): (divlength * j )])[2:7]
      par[1] <- par[1] - (gamma/10) * difH(par, y[(1 + (j - 1) * divlength ): (divlength * j )])[1]
      if(!is.null(cb)) cb()
      criteria <- value - H(par, x) - epsilon * (H(par,x) + epsilon)
      if(criteria <= 0) {
        break
      }
    }
  }
  par
}





########








##

# random initializiation could be 

c(runif(1), rexp(6))

# Random initialization. THe problem is that one might get out of bound and in the other hand one has to let 
# the learning rate go down a lot. Else one could perhaps let the learning be low in the beginning and then grow.

nasgd_tracer <- tracer(c("criteria"), N = 5)

nagd_tracer <- tracer(c("criteria"), N = 5)

nagd(par = c(0.2, 13, 5, 4, 0.8, 3, 2), H = logLike, difH = analyticGrad, alpha = 0.9 , x = data,
     gamma = 0.000001, cb = nagd_tracer$trace)

nagd(par = c(runif(1), rexp(6) + 1), H = logLike, difH = analyticGrad, alpha = 0.9 , x = data,
     gamma = 0.0000001, cb = nagd_tracer$trace)
nasgd(par = c(0.2, 13, 5, 4, 0.8, 3, 2), H = logLike, difH = analyticGrad, alpha = 0.9 , x = data,
      gamma = 0.000001, cb = nasgd_tracer$trace)

library(ggplot2)

nasgd_tracer_data <- summary(nasgd_tracer)

nagd_tracer_data <- summary(nagd_tracer)

head(nasgd_tracer_dat)

ggplot(data = nasgd_tracer_dat, aes(.time, criteria)) + geom_point(col = "red") + 
  geom_point(data = nagd_tracer_data, aes(.time, criteria)) + scale_y_log10() 


#quite good for nesterov accelarated gradient. 

#But takes 2 minutes and 28 seconds. 
gdnag(par = c(0.2, 13, 5, 4, 0.8, 3, 2), H = logLike, difH = analyticGrad, alpha = 0.9 , x = data,
      gamma = 0.000001)


# sgd(par = c(0.3, 10, 7, 3, 0.5, 2, 4), x = data, divlength = 50, maxiter = 100, gamma = 0.00001)

sgd2(par = c(0.301, 10, 7, 3, 0.5, 2, 4), x = data, divlength = 50, maxiter = 100, gamma = 0.000001)

#speed depends on which local minima it finds. 

sgd2(par = c(0.02, 13, 5, 4, 0.8, 3, 2), logLike, analyticGrad,
     x = data, divlength = 50, gamma = 0.000001)

library(numDeriv)


gd_tracer <- tracer(c("criteria"), N = 5)



