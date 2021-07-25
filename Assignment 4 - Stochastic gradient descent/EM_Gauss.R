
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
