EKernel <- function(x){
  ifelse(-1 <= x & x<= 1,(3/4)*(1-x^2), 0) 
}

bin <- function (x, lo, hi, m) {
  w <- numeric(m)
  delta <- (hi - lo)/(m - 1) # bin length
  for (i in seq_along(x)) {
    ii <- floor((x[i] - lo)/delta + 0.5) + 1 # determines which bin xi is in
    w[ii] <- w[ii] + 1 # counts the number of observations in each bin
  }
  w <- w/sum(w) #the weights are computed as the number of observatios in each bin devided by n
}

kernDens_bin <- function (x, h, m = 512, kernel) {
  rg <- range(x) + c(-3 * h, 3 * h)
  xx <- seq(rg[1], rg[2], length.out = m) # evaluation grid around x : m points 
  weights <- bin(x, rg[1], rg[2], m) # compute weights 
  kerneval <- kernel((xx - xx[1])/h)/h
  kerndif <- toeplitz(kerneval)
  y <- colSums(weights * kerndif)
  list(x = xx, y = y) # return evaluation grid and function values 
}

