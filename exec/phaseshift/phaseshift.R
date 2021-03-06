## omega_lm
## this needs qtilde as input!
omegalm <- function(l=0, m=0, q, gamma=1, dvec=c(0,0,0)) {
  return( LuescherZeta(qsq = q^2, l=l, m=m, gamma=gamma, dvec=dvec)/(pi^(3/2)*sqrt(2*l+1)*q^(l+1)*gamma) )
}

## center of mass energy from q^2 and Mpi
Ecmofqsq <- function(q, Mpi) {
  return(2*acosh(cosh(Mpi) +  2*sin(q/2)^2 ))
}
Ecmofqsq.cont <- function(q, Mpi) {
  return(2*sqrt(Mpi^2 + q^2))
}


## energy level in MF from Ecm
EofEcm <- function(Ecm, dvec, L) {
  return(acosh( cosh(Ecm) + 2*sum(sin(pi*dvec/L)^2) ))
}
EofEcm.cont <- function(Ecm, dvec, L) {
  return(sqrt(Ecm^2 + sum((2*pi*dvec/L)^2)))
}


## computes the boost factor as a function of q^2
gammaofq <- function(q, dvec, Mpi, L) {
  if(lat.disp) {
    Ecm <- Ecmofqsq(q, Mpi)
    return(EofEcm(Ecm, dvec, L)/Ecm)
  }
  Ecm <- Ecmofqsq.cont(q, Mpi)
  ## E/Ecm
  return(EofEcm.cont(Ecm, dvec, L)/Ecm)
}

## computes the energy level given q^2
Eofqsq <- function(q, dvec, Mpi, L) {
  if(lat.disp) {
    Ecm <- Ecmofqsq(q, Mpi)
    return(EofEcm(Ecm, dvec, L))
  }
  Ecm <- Ecmofqsq.cont(q, Mpi)
  return(EofEcm.cont(Ecm, dvec, L))
}

## center of mass frame formulae

## this function is indendet for a rough scan in q
## and it is independent of the irrep
prepdetEqCMscan <- function(q, Mpi, L) {
  W <- list()
  W$dvec <- c(0,0,0)
  W$gamma <- gammaofq(q, dvec=W$dvec, Mpi=Mpi, L=L)
  W$qt <- L*q/2/pi
  W$L <- L
  W$Mpi <- Mpi
  W$q <- q

  W$w00 <- Re(omegalm(l=0, m=0, q=W$qt, gamma=W$gamma))
  W$w40 <- Re(omegalm(l=4, m=0, q=W$qt, gamma=W$gamma))
  return(W)
}

## A1 irrep
detEqCMA1scan <- function(par, W) {
  cd0 <- par[1] + 0.5*par[2]*W$q^2
  return(cd0 - W$q*W$w00)
}

detEqCMA1 <- function(q, par, Mpi, L) {
  gamma <- rep(1., times=(length(q)))
  qt <- L*q/2/pi
  cd0 <- par[1] + 0.5*par[2]*q^2
  return(cd0 - q*Re(omegalm(l=0, m=0, q=qt, gamma=gamma)))
}

## E irrep
detEqCMEscan <- function(par, W) {
  cd2 <- par[3]
  return(cd2 - W$q^5*(W$w00 + 18/7*W$w40))
}

detEqCME <- function(q, par, Mpi, L) {
  gamma <- rep(1., times=(length(q)))
  qt <- L*q/2/pi
  w00 <- Re(omegalm(l=0, m=0, q=qt, gamma=gamma))
  w40 <- Re(omegalm(l=4, m=0, q=qt, gamma=gamma))
  cd2 <- par[3]
  return(cd2 - q^5*(w00 + 18/7*w40))
}

## T2 irrep
detEqCMT2scan <- function(par, W) {
  cd2 <- par[3]
  return(cd2 - W$q^5*(W$w00 - 12/7*W$w40))
}

detEqCMT2 <- function(q, par, Mpi, L) {
  gamma <- rep(1., times=(length(q)))
  qt <- L*q/2/pi
  w00 <- Re(omegalm(l=0, m=0, q=qt, gamma=gamma))
  w40 <- Re(omegalm(l=4, m=0, q=qt, gamma=gamma))
  cd2 <- par[3]
  return(cd2 - q^5*(w00 - 12/7*w40))
}

## moving frame with d=(0,0,1)
prepdetEqMF1scan <- function(q, Mpi, L) {
  W <- list()
  W$dvec <- c(0,0,1)
  W$gamma <- gammaofq(q, dvec=W$dvec, Mpi=Mpi, L=L)
  W$qt <- L*q/2/pi
  W$L <- L
  W$Mpi <- Mpi
  W$q <- q

  W$w00 <- Re(omegalm(l=0, m=0, q=W$qt, gamma=W$gamma, dvec=W$dvec))
  W$w20 <- Re(omegalm(l=2, m=0, q=W$qt, gamma=W$gamma, dvec=W$dvec))
  W$w40 <- Re(omegalm(l=4, m=0, q=W$qt, gamma=W$gamma, dvec=W$dvec))

  return(W)
}

## A1 irrep
detEqMF1A1scan <- function(par, W) {
  cd0 <- par[1] + 0.5*par[2]*W$q^2
  if(length(par) < 3) {
    return(cd0 - W$q*W$w00)
  }
  cd2 <- par[3]
  return( (cd0 - W$q*W$w00)*
         (cd2 - W$q^5*(W$w00 + 10/7*W$w20 + 18/7*W$w40))
         -W$q^6*5*W$w20^2
         )
}

detEqMF1A1 <- function(q, par, Mpi, L) {
  dvec <- c(0,0,1)

  gamma <- gammaofq(q, dvec=dvec, Mpi=Mpi, L=L)

  qt <- L*q/2/pi
  w00 <- Re(omegalm(l=0, m=0, q=qt, gamma=gamma, dvec=dvec))
  cd0 <- par[1] + 0.5*par[2]*q^2
  if(length(par)<3) {
    return(cd0 - q*w00)
  }
  w20 <- Re(omegalm(l=2, m=0, q=qt, gamma=gamma, dvec=dvec))
  w40 <- Re(omegalm(l=4, m=0, q=qt, gamma=gamma, dvec=dvec))

  cd2 <- par[3]
  return( (cd0 - q*w00)*
         (cd2 - q^5*(w00 + 10/7*w20 + 18/7*w40))
         -q^6*5*w20^2
         )
}

## moving frame with d=(1,1,0)
prepdetEqMF2scan <- function(q, Mpi, L) {
  W <- list()
  W$dvec <- c(1,1,0)
  W$gamma <- gammaofq(q, dvec=W$dvec, Mpi=Mpi, L=L)
  W$qt <- L*q/2/pi
  W$L <- L
  W$Mpi <- Mpi
  W$q <- q

  W$w00 <- Re(omegalm(l=0, m=0, q=W$qt, gamma=W$gamma, dvec=W$dvec))
  W$w20 <- Re(omegalm(l=2, m=0, q=W$qt, gamma=W$gamma, dvec=W$dvec))
  W$w40 <- Re(omegalm(l=4, m=0, q=W$qt, gamma=W$gamma, dvec=W$dvec))
  W$w22 <- omegalm(l=2, m=2, q=W$qt, gamma=W$gamma, dvec=W$dvec)
  W$w44 <- Re(omegalm(l=4, m=4, q=W$qt, gamma=W$gamma, dvec=W$dvec))
  W$w42 <- omegalm(l=4, m=2, q=W$qt, gamma=W$gamma, dvec=W$dvec)

  return(W)
}

## A1 irrep
detEqMF2A1scan <- function(par, W) {
  cd0 <- (par[1]/W$q + 0.5*par[2]*W$q)
  if(length(par) < 3 ) {
    return(cd0 - W$q*W$w00)
  }

  cd2 <- par[3]/W$q^5

  return(Re(W$q^11*(10/49*(10*(-2*cd0 + 2*W$w00 + 7*W$w20)*W$w22^2 +
                         3*sqrt(15)*(4*cd0 - 4*W$w00 - 7*W$w20)*W$w22*W$w42 + 27*(W$w00 - cd0)*W$w42^2)
                  +(-cd2 + W$w00 + 2/7*(5*W$w20 + 9*W$w40))*(10*W$w22^2 + 1/7*(cd0 - W$w00)*(7*cd2 - 7*W$w00 + 10*W$w20 - 3*W$w40 + 3*sqrt(70)*W$w44)) +
                  5/7*W$w20*(20*W$w22^2 - 6*sqrt(15)*W$w22*W$w42 + W$w20*(7*cd2 - 7*W$w00 + 10*W$w20 - 3*W$w40 + 3*sqrt(70)*W$w44))
                  )
            ))
}

detEqMF2A1 <- function(q, par, Mpi, L) {
  dvec=c(1,1,0)
  gamma <- gammaofq(q, dvec=dvec, Mpi=Mpi, L=L)

  qt <- L*q/2/pi
  w00 <- Re(omegalm(l=0, m=0, q=qt, gamma=gamma, dvec=dvec))
  cd0 <- (par[1]/q + 0.5*par[2]*q)
  if(length(par) < 3) {
    return(cd0 - q*w00)
  }
  w20 <- Re(omegalm(l=2, m=0, q=qt, gamma=gamma, dvec=dvec))
  w22 <- omegalm(l=2, m=2, q=qt, gamma=gamma, dvec=dvec)
  w40 <- Re(omegalm(l=4, m=0, q=qt, gamma=gamma, dvec=dvec))
  w44 <- Re(omegalm(l=4, m=4, q=qt, gamma=gamma, dvec=dvec))
  w42 <- omegalm(l=4, m=2, q=qt, gamma=gamma, dvec=dvec)

  cd2 <- par[3]/q^5
  
  return(Re(q^11*(10/49*(10*(-2*cd0 + 2*w00 + 7*w20)*w22^2 +
                      3*sqrt(15)*(4*cd0 - 4*w00 - 7*w20)*w22*w42 + 27*(w00 - cd0)*w42^2)
               +(-cd2 + w00 + 2/7*(5*w20 + 9*w40))*(10*w22^2 + 1/7*(cd0 - w00)*(7*cd2 - 7*w00 + 10*w20 - 3*w40 + 3*sqrt(70)*w44)) +
               5/7*w20*(20*w22^2 - 6*sqrt(15)*w22*w42 + w20*(7*cd2 - 7*w00 + 10*w20 - 3*w40 + 3*sqrt(70)*w44))
               )
         ))
}

## moving frame with d=(1,1,1)
prepdetEqMF3scan <- function(q, Mpi, L) {
  W <- list()
  W$dvec <- c(1,1,1)
  W$gamma <- gammaofq(q, dvec=W$dvec, Mpi=Mpi, L=L)
  W$qt <- L*q/2/pi
  W$L <- L
  W$Mpi <- Mpi
  W$q <- q

  W$w00 <- Re(omegalm(l=0, m=0, q=W$qt, gamma=W$gamma, dvec=W$dvec))
  W$w40 <- Re(omegalm(l=4, m=0, q=W$qt, gamma=W$gamma, dvec=W$dvec))
  W$w22 <- omegalm(l=2, m=2, q=W$qt, gamma=W$gamma, dvec=W$dvec)
  W$w42 <- omegalm(l=4, m=2, q=W$qt, gamma=W$gamma, dvec=W$dvec)

  return(W)
}

## A1 irrep
detEqMF3A1scan <- function(par, W) {
  cd0 <- (par[1] + 0.5*par[2]*W$q^2)
  if(length(par) < 3) {
    return(cd0 - W$q*W$w00)
  }
  cd2 <- par[3]

  return(Re(( cd0 - W$q*W$w00 )*
            ( cd2 - W$q^5*(W$w00 - 12/7*W$w40 - 12*sqrt(10)/7*1i*W$w42 - 10*sqrt(6)/7*1i*W$w22)) +
            W$q^6*30*W$w22^2
            )
         )
}

detEqMF3A1 <- function(q, par, Mpi, L) {
  dvec=c(1,1,1)
  gamma <- gammaofq(q, dvec=dvec, Mpi=Mpi, L=L)

  qt <- L*q/2/pi
  w00 <- Re(omegalm(l=0, m=0, q=qt, gamma=gamma, dvec=dvec))
  cd0 <- (par[1] + 0.5*par[2]*q^2)
  
  if(length(par) < 3) {
    return(cd0 - q*w00)
  }
  w40 <- Re(omegalm(l=4, m=0, q=qt, gamma=gamma, dvec=dvec))
  w22 <- omegalm(l=2, m=2, q=qt, gamma=gamma, dvec=dvec)
  w42 <- omegalm(l=4, m=2, q=qt, gamma=gamma, dvec=dvec)

  cd2 <- par[3]
  return(Re(( cd0 - q*w00 )*
            ( cd2 - q^5*(w00 - 12/7*w40 - 12*sqrt(10)/7*1i*w42 - 10*sqrt(6)/7*1i*w22)) +
            q^6*30*w22^2
            )
         )
}

findSignChanges <- function(fn, makeplot=FALSE, par, W, no=3, threshold=100) {
  res <- fn(par=par, W=W)
  if(makeplot) {
    if(interactive()) X11()
    plot(W$q^2, res, type="l", ylim=c(-threshold, threshold))
    abline(a=0, b=0, col="red")
  }
  
  ii <- integer(0)
  j <- 1
  for(i in c(2:(length(W$q)-1))) {
    ## find all sign-changes 
    ##cat(res[i], res[i+1], "\n")
    if(!is.na(res[i]) && !is.na(res[i+1])) {
      if(res[i]/abs(res[i]) != res[i+1]/abs(res[i+1])) {
        ## remove the poles
        if((res[i] < 0 && res[i-1] < res[i]) ||
           (res[i] > 0 && res[i-1] > res[i])) {
          ii[j] <- i
          if(j == no) break
          j <- j+1
        }
      }
    }
  }
  return(ii)
}

findZeros <- function(fn, q, ii, makeplot=FALSE, tol=1.e-14, ...) {

  ## now determine the root more precisely
  zeros <- numeric(0.)
  for(j in c(1:length(ii))) {
    z <- uniroot(fn, interval=c(q[ii[j]], q[ii[j]+1]), tol = tol, ...)
    if(z$f.root < 1) {
      zeros[j] <- z$root
      ##cat("root", j, ":", z$root, z$f.root, z$estim.prec, "\n")
    }
  }
  if(makeplot) {
    abline(v=zeros^2, col="blue")
  }

  return(zeros)
}

getEvalues <- function(par, Wlist, Lvalues, nolist, framelist, scanfnlist, fnlist, irreplist) {
  evalues <- c()
  j <- 1
  k <- 1
  for(L in Lvalues) {
    for(f in framelist) {
      for(i in irreplist) {
        if(nolist[[k]] > 0) {
          ii <- findSignChanges(fn=scanfnlist[[i]][[f]], par=par, makeplot=FALSE, W=Wlist[[j]], no=nolist[[k]])
          
          if(length(ii) != nolist[[k]]) return(NA)
          zeros <- findZeros(fn=fnlist[[i]][[f]], q=Wlist[[j]]$q, ii=ii, L=L, Mpi=Wlist[[j]]$Mpi, par=par, makeplot=FALSE)
          if(length(zeros) != nolist[[k]]) return(NA)
          evalues <- c(evalues, Eofqsq(zeros, dvec=Wlist[[j]]$dvec, Mpi=Wlist[[j]]$Mpi, L=L))
        }
        k <- k+1
      }
      j <- j+1
    }
  }
  return(evalues)
}

getEDiff <- function(x, par, index, level, scanfn, dfn, W, E) {
  par[index] <- x
  ii <- findSignChanges(fn=scanfn, par=par, makeplot=FALSE, W=W, no=level)
  if(length(ii) != level) return(NA)
  zero <- findZeros(fn=dfn, q=W$q, ii=ii[level], L=W$L, Mpi=W$Mpi, par=par, makeplot=FALSE)
  if(length(zero) != 1) return(NA)
  return(Eofqsq(zero, dvec=W$dvec, Mpi=W$Mpi, L=W$L)-E)
}

getq <- function(par, level, scanfn, dfn, W, E) {
  ii <- findSignChanges(fn=scanfn, par=par, makeplot=FALSE, W=W, no=level)
  if(length(ii) != level) return(NA)
  zero <- findZeros(fn=dfn, q=W$q, ii=ii[level], L=W$L, Mpi=W$Mpi, par=par, makeplot=FALSE)
  if(length(zero) != 1) return(NA)
  return(zero)  
}

solveforCotDelta <- function(par, index, level, scanfn, fn, W, E, interval) {
  n <- 0
  if(missing(interval)) {
    interval <- c(-pi/2.,-0.1)
    if(index == 3) interval <- c(-0.025,-0.001)
  }
  l <- getEDiff(x=interval[1], par, index=index, level=level, scanfn=scanfn, dfn=fn, W=W, E=E)
  r <- getEDiff(x=interval[2], par, index=index, level=level, scanfn=scanfn, dfn=fn, W=W, E=E)
  while(is.na(l) || is.na(r) || l/abs(l) == r/abs(r)) {
    if(n > 50*index) {
      cat("NA\n")
      return(NA)
    }
    if(index == 1)  interval[1] <- interval[1]*9/10.
    if(index == 3)  interval[1] <- interval[1]*1.1
    l <- getEDiff(x=interval[1], par, index=index, level=level, scanfn=scanfn, dfn=fn, W=W, E=E)
    r <- getEDiff(x=interval[2], par, index=index, level=level, scanfn=scanfn, dfn=fn, W=W, E=E)
    n <- n+1
  }
  qcotdelta <- uniroot(f = getEDiff, interval=interval, tol = 1.e-14, par=par, index=index, level=level, scanfn=scanfn, dfn=fn, W=W, E=E)
  cat("level", level, qcotdelta$root, qcotdelta$f.root, qcotdelta$estim.prec, "\n")
  return(qcotdelta$root)
}
