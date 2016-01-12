#' Compute the average percent errors
#'
#' @description
#' The function \code{percent.error} computes the average percent errors for
#' the rows or columns of a matrix.
#'
#' @param Obs A matrix of observations.
#' @param Est A matrix of estimates.  Must be of the same dimensions as
#'  \code{Obs}.
#' @param MARGIN (optional) The margin across which to compute the
#' efficiencies.  The default is to compute an efficiency for each column
#' (\code{MARGIN=2}).
#'
#' @details
#' Lorem ipsum...
#'
#' @return A vector of average percent errors.
#'@export
percent.error <- function(Obs,Est,MARGIN=2) {
  # Function orginally developed by William Farmer, 27 May 2015

  # 26 August 2015, WHF: Added code to avoid observations that cause INF.
  Est[which(!is.finite((Obs-Est)^2))] <- NA
  Obs[which(!is.finite((Obs-Est)^2))] <- NA

  if (MARGIN==1) {
    Obs <- t(Obs)
    Est <- t(Est)
  }

  perr <- colMeans((Est-Obs)/Obs,na.rm=T)
  return(perr)
}
