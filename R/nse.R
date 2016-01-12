
#' Compute the Nash-Sutcliffe Efficiency
#'
#' @description
#' The function \code{nse} computes the Nash-Sutcliffe efficiency for the rows
#' or columns of a matrix.
#'
#' @param Obs A matrix of observations.
#' @param Est A matrix of estimates.  Must be of the same dimensions as
#'  \code{Obs}.
#' @param MARGIN (optional) The margin across which to compute the
#' efficiencies.  The default is to compute an efficiency for each column
#' (\code{MARGIN=2}).
#' @param na.rm (optional)  A logical indicating if \code{NA}s should be
#' removed.  The default behavior is to ignore \code{NA}s.
#'
#' @details
#' Lorem ipsum...
#'
#' @return A vector of Nash-Sutcliffe efficiencies.
#'@export
nse <- function(Obs,Est,MARGIN=2,na.rm=T) {
  # Function orginally developed by William Farmer, 27 May 2015

  # 26 August 2015, WHF: Added code to avoid observations that cause INF.
  Est[which(!is.finite((Obs-Est)^2))] <- NA
  Obs[which(!is.finite((Obs-Est)^2))] <- NA

  if (na.rm) {
    NA.ndx <- is.na(Obs+Est)
    Est[NA.ndx] <- NA
    Obs[NA.ndx] <- NA
  }
  if (is.vector(Obs)&is.vector(Est)) {
    NSE <- 1 - sum((Obs-Est)^2,na.rm=na.rm)/
      sum((Obs-mean(Obs,na.rm=na.rm))^2,na.rm=na.rm)
    if (sum(is.na(Obs))==length(Obs)) {NSE<-NA}
  } else if (is.matrix(Obs)&is.matrix(Est)) {
    if (MARGIN==2) {
      Obs <- t(Obs)
      Est <- t(Est)
    }
    NSE <- 1 - rowSums((Obs-Est)^2,na.rm=na.rm)/
      rowSums((Obs-rowMeans(Obs,na.rm=na.rm))^2,na.rm=na.rm)
    ndx <- which(rowSums(is.na(Obs))==ncol(Obs))
    NSE[ndx] <- NA
  }
  return(NSE)
}
