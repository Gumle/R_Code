
envelope_factory <- function(x,  target_dens, logderiv = NULL,lower_support = -Inf,upper_support = Inf)
{
  
  #check if the log-derivative was supplied, otherwise we get it ourselves:
  if (is.null(logderiv))
  {
    logderiv <- function(xx) grad(function(x) log(target_dens(x)), xx)
  }
  #calculate a-vector and stop if there are integrability issues
  a <- logderiv(x)
  #check if we need to stop
  continue <- (a[1] > 0 & a[length(a)] < 0) | (a[length(a)] > 0 & lower_support > -Inf) | a[1] < 0 & upper_support < Inf
  if (!continue)
  {
    stop("Envelope is not integrable. Re-submit new x")
  }
  #now calculate b, z, Fz, Q, and const
  b <- log(target_dens(x)) - a*x
  z <- c(lower_support, -diff(b)/diff(a), upper_support)
  Fz <- numeric(length(x))
  for (i in seq_along(Fz))
  {
    Fz[i] <- exp(b[i])*( exp( a[i] * z[i+1]) - exp(a[i]*z[i])) / a[i]
  }
  
  Q <- c(0, cumsum(Fz))
  const <- Q[length(Q)]
  
  #now define simulator and envelope -- must be vectorized
  
  proposal_density <- function(x) index <- findInterval(x, z)  #Given x, find z-interval that x belongs to
  
   
  #findInterval is optimized and O(log(length(z)*length(x))
  
  #now just evaluate the function and return
  V <- a[index]*x + b[index]
  exp(V)/const

  simulator <- function(n)
  {
  #we need n uniform samples
  q <- runif(n)
  
  #find the index to which const*q belongs to. 
  index <- findInterval(const*q, Q) + 1
  
  #solve the equation for x
  log( (const*q - Q[index])*a[index ]*exp(-b[index])  + exp(a[index]*z[index]))/a[index ]
  }  

  #return the list-object that can be used directly in previous functions.
  list( tdens = dens$tdens, pdens = proposal_density, sim = simulator, alpha = 1/const)
}

