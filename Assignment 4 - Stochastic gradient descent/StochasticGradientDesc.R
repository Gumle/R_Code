## Least squares objective function, gradient(s) and hessian

H <- function(beta) 
  drop(crossprod(y - X %*% beta)) / (2 * nrow(X))

grad_H <- function(beta) { 
  res <- y - X %*% beta -crossprod(X, res) / nrow(X)
}

grad_h <- function(beta, i) { 
  X <- X[i, , drop = FALSE]
  res <- y[i] - X %*% beta -crossprod(X, res) / nrow(X)
}

grad_h_slow <- function(beta, i) { 
  X <- X[i, , drop = FALSE]
  res <- y[i] - X %*% beta -t(X) %*% res / nrow(X)
}

hessian_h <- function(beta, i) {
  X <- X[i, , drop = FALSE]
  crossprod(X) / nrow(X)
}


## SDG implementation with only stopping criterion being max-iterations

SGD <- function(par, grad_h,
                n,                 ## Sample size
                mbn = 50,          ## Minibatch size
                gamma = 1e-4,      ## Learning rate
                maxiter = 100,     ## Max epoch iterations
                cb = NULL) {
  m <- floor(n / mbn)
  for(k in 1:maxiter) {
    samp <- sample(n)
    for(j in 0:(m-1)) {
      i <-  samp[(j * mbn + 1):(j * mbn + mbn)]
      par <- par - gamma * grad_h(par, i)
    }
    if(!is.null(cb)) cb()
  }
  par
}

## GD from notes for comparison

GD <- function(par, H, grad_H, d = 0.8,  c = 0.1,  gamma0 = 0.01,   epsilon = 1e-4,  cb = NULL) {
  repeat {
    value <- H(par)
    grad <- grad_H(par)
    h_prime <- sum(grad^2)
    if(!is.null(cb)) cb()
    ## Convergence criterion based on gradient norm
    if(h_prime <= epsilon) break
    gamma <- gamma0
    ## First proposed descent step
    par1 <- par - gamma * grad
    ## Backtracking while descent is insufficient
    while(H(par1) > value - c * gamma * h_prime) {
      gamma <- d * gamma
      par1 <- par - gamma * grad
    }
    par <- par1
  }
  par
}


## Nesterov implementation with only stopping criterion being max-iterations

Nesterov <- function(par, 
                     grad_h,
                     n,              ## Sample size
                     mbn = 50,       ## Minibatch size
                     gamma = 1e-4,      ## Learning rate
                     alpha = 0.9,    ## Momentum fraction
                     maxiter = 100,  ## Epoch iterations
                     cb = NULL) {
  m <- floor(n / mbn)
  v <- grad_h(par, sample(n, m))
  for(k in 1:maxiter) {
    samp <- sample(n)
    for(j in 0:(m-1)) {
      i <-  samp[(j * mbn + 1):(j * mbn + mbn)]
      par0 <- par - alpha * v
      v <- alpha * v + gamma * grad_h(par0, i)
      par <- par - v
    }
    if(!is.null(cb)) cb()
  }
  par
}
