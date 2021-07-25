CG <- function(par, H, grad_H, maxiter = 1000, d = 0.9, c = 0.1, gamma0 = 1, epsilon = 1e-7, stop = 'func', cb = NULL, ...) {
  n = m = 1; p = length(par); rho0 = numeric(p)
  repeat {
    H_val = H(par, ...) # Status of objective
    grad = grad_H(par, ...) # function and gradient
    grad_norm = sum(grad^2)
    gamma = gamma0 # Step size
    rho = - grad + grad_norm * rho0 # Descent direction
    h_prime = drop(t(grad) %*% rho) # H in descent direc.
    if(m > p || h_prime >= 0) { # Reset
      rho = - grad
      h_prime = - grad_norm
      m = 1
    }
    par1 = par + gamma * rho
    while(H(par1, ...) > H_val + c * gamma * h_prime){ # Backtracking
      gamma = d * gamma
      par1 = par + gamma * rho
    }
    H_val1 = H(par1, ...); [...]; rho0 = rho / grad_norm;
    par = par1; n = n + 1; m = m + 1 # Updating
  }
  list(par = par, H_val = H_val1, grad_norm = grad_norm, stop = stopping)
}