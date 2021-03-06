% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/raw_cf.R
\name{idx_matrix.raw_cf}
\alias{idx_matrix.raw_cf}
\title{Construct the tensor index set for the entire raw correlator}
\usage{
idx_matrix.raw_cf(cf, component)
}
\arguments{
\item{cf}{'raw_cf' container with data and meta-data}

\item{component}{Integer vector. Optional argument to obtain a subset of the
index matrix to access a particular element of the interior
dimensions. Must of the the same length as cf$dim.}
}
\description{
Construct the tensor index set for the entire raw correlator
}
