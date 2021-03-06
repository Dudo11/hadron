% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/new_matrixfit.R
\name{make_parind}
\alias{make_parind}
\title{Create a parameter index matrix for \code{matrixfit}}
\usage{
make_parind(parlist, length_time, summands = 1)
}
\arguments{
\item{parlist}{integer array. Parameter list generated with \code{make_parlist}.}

\item{length_time}{integer. Number of time slices per correlator.}

\item{summands}{integer. Number of summands in the fit model that shall be
fitted. The signal counts as one summand, each explicit pollution term with
independent amplitudes counts as its own summand.}
}
\description{
Create a parameter index matrix for \code{matrixfit}
}
