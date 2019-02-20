#' List that constains the geometric constants I, J and K that appear in the
#' Bean-Savage formula for the n-particle energy-shift.

geopar <- list(I = (-8.9136329), J = 16.532316, K = 8.4019240)

#' prediction_function `fn(par, x, ...)`. The non-linear energy-shift function to be
#' fitted to the data. Its first argument must be the fit parameters named \code{par}.
#' The second must be \code{x}, the explaining variable.

prediction_function <- function(par, x, M, n, ...) {
  f <- choose(n, 2) * ((4*pi*par[1])/(M * x^3)) * ( 1 - par[1] * (geopar$I/(pi*x)) +
        (par[1]/(pi* x))^2 * ((geopar$I)^2 + (2*n - 5)*geopar$J) -
        (par[1]/(pi*x))^3 * ((geopar$I)^3 + (2*n - 7)*geopar$I*geopar$J +
        (5*n^2 - 41*n + 63)*geopar$K)) +
        choose(n, 2) * ((8*pi^2 * par[1]^3) / (M * x^6)) * par[2]
  return(f)
}

#' Gradient of the prediction_function with respect to the scattering length.

prediction_gradient <- function(par, x, M, n, ...) {
  df1 <- choose(n, 2) * ((4*pi)/(M * x^3)) * ( 1 - 2*par[1] * (geopar$I/(pi*x)) +
          3*(par[1]/(pi* x))^2 * ((geopar$I)^2 + (2*n - 5)*geopar$J) -
          4*(par[1]/(pi*x))^3 * ((geopar$I)^3 + (2*n - 7)*geopar$I*geopar$J +
          (5*n^2 - 41*n + 63)*geopar$K)) +
          choose(n, 2) * ((24*pi^2 * par[1]^2)/(M*x^6)) * par[2]
  df2 <- choose(n, 2) * ((8*pi^2 * par[1]^3) / (M * x^6))
  res <- cbind(df1, df2)
  return(res)
}

#' Fit the n-particle energy-shift function by Bean and Savage to provided data
#'
#' Performs a bootstrap.nlsfit for the relativistic n-particle Bean-Savage formula.
#'
#' @param y the energy-shift data as a one-dimensional numerical vector to be described by
#' the fit function. This vector can be computed e.g. by means of calculate.energy.shift.
#' @param bsamples bootstrap samples of \code{y}.
#' Must be provided as array of dimensions \code{c(boot.R, n)} with \code{n}
#' equals to \code{length(y)}.
#' These samples can be computed e.g. by means of calculate.energy.shift.
#' @param l1 lower limit of the fit range.
#' @param l2 upper limit of the fit range.
#' @param n number of particles and index of the calculated energy shift.
#' @param M mass of the particles.
#' @param ain initial fit parameter for the scattering length.


fit.scattering.length <- function(bsamples,
                                  y,
                                  lvec,
                                  l1, l2,
                                  n,
                                  M,
                                  ain,
                                  useCov = FALSE,
                                  boot.fit = TRUE,
                                  fit.method = "optim",
                                  autoproceed = FALSE,
                                  every,
                                  priors = list(param = c(), p = c(), psamples = c()),
                                  cov_fn = cov,
                                  error = sd,
                                  ...) {

  fitfun <- prediction_function
  dfitfun <- prediction_gradient
  
  ## Initial guess for the fit parameter
  initial_guess = function(npar){
    par <- numeric(npar)
    par <- ain
    return (par)
    
  }
  par.guess <- initial_guess(npar = 2)
  
  Llen <- length(x)
  
  ## Number of energy shifts
    mSize <- dim(bsamples)[2] / Llen
  
  ## Index vector for timeslices to be fitted
  i <- c((l1):(l2))
  
  args <- list(fn = fitfun,
               gr = dfitfun,
               par.guess = par.guess,
               y = y[i],
               x = lvec[i],
               bsamples = bsamples[, i],
               use.minpack.lm = fit.method == 'lm',
               error = error,
               cov_fn = cov_fn,
               M = M,
               n = n,
               priors = priors,
               ...)
  
  if (useCov) {
    if(!missing(priors)){
      args$CovMatrix <- cov_fn(cbind(bsamples[, i], priors$psamples))
    }
    else{
      args$CovMatrix <- cov_fn(bsamples[, i])
    }
  }
  
  res <- do.call(bootstrap.nlsfit, args)
}