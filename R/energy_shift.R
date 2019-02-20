#' Calculate an energy shift
#'
#' Calculates the n particle energy shift \code{\Delta E_n}.
#'
#' @param psamples energy-parameter samples vector, i.e. the first column of the
#' matrix t with length \code{\length{boot.R}} that is returned from the n particle
#' fit (n > 1) with the bootstrap.nlsfit function.
#' @param msamples mass samples vector, i.e. the first column of the matrix t with
#' length \code{\length{boot.R}} that is returned from the one particle fit with
#' the bootstrap.nlsfit function.
#' @param n integer. It denotes the number of particles for which the energy shift
#' is calculated. 
#' @param error. The default error is the standard deviation, which can be applied,
#' i.e. in the case that msamples and psamples are bootstrap samples.
#' @param print.samples option for the output. If print.samples is TRUE, the energy-
#' shift samples are returned. Otherwise the energy-shift values and errors are
#' returned.
calculate.energy.shift <- function(psamples,
                                  msamples,
                                  n,
                                  error = sd,
                                  print.samples = FALSE) {
  stopifnot(!missing(psamples))
  stopifnot(!missing(msamples))
  stopifnot(!missing(n))
  
  if(n > 1){
    m <- as.matrix(psamples - n * msamples)
    DeltaE <- apply(X = m, MARGIN = 2, FUN = mean)
    dDeltaE <- apply(X = m, MARGIN = 2, FUN = error)
    
    if(print.samples == FALSE){
      res <- list(DeltaE = DeltaE, dDeltaE = dDeltaE, n = n)
      print.energyshift(object = res)
    }
    else{
      return(m)
    }
  }
  else{
    stop(sprintf('The number of particles must be greater than 1.'))
  }
}

print.energyshift <- function(object, ..., digits = 2) {
  values <- object$DeltaE
  errors <- object$dDeltaE
  tmp <- apply(X=array(c(values, errors), dim=c(length(values), 2)), MARGIN=1, FUN=tex.catwitherror, with.dollar=FALSE, digits=digits, human.readable=FALSE)
  cat("calculate energy shift\n\n")
  cat("energy shift Delta E_", object$n, "=", tmp, "\n")
}