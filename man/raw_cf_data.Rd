% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/raw_cf.R
\name{raw_cf_data}
\alias{raw_cf_data}
\title{Original data mixin constructor for \code{raw_cf}}
\usage{
raw_cf_data(cf, data)
}
\arguments{
\item{cf}{\code{raw_cf} object to extend.}

\item{data}{Numeric or complex array, original data for all observables and measurements.
This should have dimensions c(Nmeas,cf$Time\emph{cf$nrObs}cf$nrStypes,cf$dim).
Having the internal dimensions innermost is not as efficient, but it allows
different transformations to be applied to different observables in the same
container more easily.}
}
\description{
Original data mixin constructor for \code{raw_cf}
}
\seealso{
Other raw_cf constructors: \code{\link{raw_cf_meta}},
  \code{\link{raw_cf}}
}
\concept{raw_cf constructors}
