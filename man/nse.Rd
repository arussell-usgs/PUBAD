% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nse.R
\name{nse}
\alias{nse}
\title{Compute the Nash-Sutcliffe Efficiency}
\usage{
nse(Obs, Est, MARGIN = 2, na.rm = T)
}
\arguments{
\item{Obs}{A matrix of observations.}

\item{Est}{A matrix of estimates.  Must be of the same dimensions as
\code{Obs}.}

\item{MARGIN}{(optional) The margin across which to compute the
efficiencies.  The default is to compute an efficiency for each column
(\code{MARGIN=2}).}

\item{na.rm}{(optional)  A logical indicating if \code{NA}s should be
removed.  The default behavior is to ignore \code{NA}s.}
}
\value{
A vector of Nash-Sutcliffe efficiencies.
}
\description{
The function \code{nse} computes the Nash-Sutcliffe efficiency for the rows
or columns of a matrix.
}
\details{
Lorem ipsum...
}
