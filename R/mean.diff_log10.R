#' Compute the average differences of base-10 logs
#'
#' @description
#' The function \code{mean_diff_log10} computes the average differences of
#' the base-10 logs of the rows or columns of a matrix.
#'
#' @param Obs A matrix of observations.
#' @param Est A matrix of estimates.  Must be of the same dimensions as
#'  \code{Obs}.
#' @param MARGIN (optional) The margin across which to compute the
#' differences.  The default is to compute mean differences for each column
#' (\code{MARGIN=2}).
#'
#' @details
#' Lorem ipsum...
#'
#' @return A vector of average differences of base-10 logs.
#'@export
mean.diff_log10 <- function(Obs,Est,MARGIN=2) {
  # Function developed by T.M. Over, 31 January 2017
  # Adapted from percent.error.R by W.H. Farmer, dated 2/4/2016

  Est[which((Est<=0) | !is.finite(Est))] <- NA
  Obs[which((Obs<=0) | !is.finite(Obs))] <- NA

  if (MARGIN==1) {
    Obs <- t(Obs)
    Est <- t(Est)
  }

  mean.diff <- colMeans((log10(Est)-log10(Obs)),na.rm=T)
  return(mean.diff)
}
