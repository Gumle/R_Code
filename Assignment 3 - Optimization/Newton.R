Newton = function(par, H, grad_H, hessian_H, maxiter = 50, d = 0.1, c = 0.1, gamma0 = 1, epsilon = 1e-6, stop = 'func', cb = NULL, ...) {
  n = 1
  repeat {
    # Status of objective function and derivatives
    H_val = H(par, ...)
    grad = grad_H(par, ...)
    grad_norm = sum(grad^2)
    hessian = hessian_H(par, ...)
    # Descent
    rho = -drop(solve(hessian, grad)) # Descent direction
    gamma = gamma0 # Step size
    par1 = par + gamma * rho # Proposed descent step
    h_prime = t(grad) %*% rho # H in descent direc.
    if(!is.null(cb)) cb() # Tracing
    while(H(par1, ...) > H_val + c * gamma * h_prime) { # Backtracking
      gamma = d * gamma
      par1 = par + gamma * rho
    }
    H_val1 = H(par1, ...); [...]; par = par1; n = n + 1 # Updating
  }
  list(par = par, H_val = H_val1, grad_norm = grad_norm, stop = stopping)
}

# Maximum number of iterations
if (n > maxiter) {
  stopping = 'maxiter'
  break
} else if (stop == 'func') {
  # Small relative descent for objective function
  if (H_val - H_val1 < epsilon * (H_val1 + epsilon)) {
    stopping = 'converged'
    break
  }
} else if (stop == 'grad') {
  # Small gradient
  if (grad_norm <= epsilon || n > maxiter) {
    stopping = 'converged'
    break
  }
} else if (stop == 'par') {
  # Small relative change in parameters
  if (sum((par1 - par)^2) <= epsilon * ((sum(par1)^2) + epsilon)) {
    stopping = 'converged'
    break
  }
}