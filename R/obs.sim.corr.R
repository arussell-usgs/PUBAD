#' Compute the correlation between observations and estimates
#'
#' @description
#' The function \code{obs.sim.corr} computes the correlation between
#' observations and estimates for the rows or columns of a matrix.
#'
#' @param Obs A matrix of observations.
#' @param Est A matrix of estimates.  Must be of the same dimensions as
#'  \code{Obs}.
#' @param MARGIN (optional) The margin across which to compute the
#' efficiencies.  The default is to compute an efficiency for each column
#' (\code{MARGIN=2}).
#' @param methods A character string of character vector indication which
#' correlation technique to use.  Identical to the input passed to \code{cor}.
#'
#' @details
#' Lorem ipsum...
#'
#' @return A vector or matrix of correlations.
#'@export
obs.sim.corr <- function(Obs,Est,MARGIN=2,
  methods=c("pearson","spearman","kendall"))
{
  # Function orginally developed by William Farmer, 27 May 2015
  if (MARGIN==1) {
    Obs <- t(Obs)
    Est <- t(Est)
  }
  n <- length(methods)
  corrs <- matrix(NA,nrow=ncol(Obs),n)
  for (i in 1:ncol(Obs)) {
    for (j in 1:n) {
      corrs[i,j] <- cor(Obs[,i],Est[,i],method=methods[j],
        use="na.or.complete")
    }
  }
  return(corrs)
}
