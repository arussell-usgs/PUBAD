
#' Compute the Root-Mean-Squared Error Statistics
#'
#' @description
#' The function \code{rmse_like} computes a series of statistics related to
#' the root-mean-squared error.
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
#' @return A matrix containing:
#' \item{rmse}{Root-mean-squared error.}
#' \item{rmsne}{Root-mean-square-normalized error.}
#' \item{cvrmse}{The coefficient of variation of the root-mean-squared error.}
#' \item{nrmse}{Normalized root-mean-squared error.}
#'@export
rmse_like <- function(Obs,Est,MARGIN=2) {
  # Function orginally developed by William Farmer, 27 May 2015

  # 26 August 2015, WHF: Added code to avoid observations that cause INF.
  Est[which(!is.finite((Obs-Est)^2))] <- NA
  Obs[which(!is.finite((Obs-Est)^2))] <- NA


  if (MARGIN==1) {
    Obs <- t(Obs)
    Est <- t(Est)
  }

  rmse <- colMeans((Obs-Est)^2,na.rm=T)
  rmsne <- colMeans(((Obs-Est)/Obs)^2,na.rm=T)
  cvrmse <- rmse/colMeans(Obs,na.rm=T)
  nrmse <- rmse

  for (i in 1:ncol(Obs)) {
    range <- max(Obs[,i],na.rm=T)-min(Obs[,i],na.rm=T)
    nrmse[i] <- rmse[i]/range
  }
  result <- cbind(rmse,rmsne,cvrmse,nrmse)
  return(result)
}
