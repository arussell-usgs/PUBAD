#' Collect the estiamted and observed streamflow statistics.
#'
#' @description
#' The function \code{getFlowStats} collects a wide array of streamflow
#' statistics for observed and modeled streamflow time series.
#'
#' @param modelOutput A list of model results, in the format returned by
#' \code{\link{estDAR}} or \code{\link{estQPPQ}}.
#' @param zero.val (optional) The censoring value used to indicate a
#' censored streamflow value.  The default is \code{0.001}.
#' @param cens.level (optional) The value below which a streamflow value is
#' considered censored.  (Rounded down rather than up.)  The default is
#' \code{0.005}.
#' @param min.obs (optional) The minimum non-censored value of streamflow.
#' The default is \code{0.01}.
#'
#' @details
#' See \link{\code{calcFlowStats}} for a list of statistics calculated.
#'
#' @return A list of two elements:
#' \item{Obs}{A data frame of streamflow statistics derived from observed
#' time series.}
#' \item{Est}{A data frame of streamflow statistics derived from estimated
#' time series.}
#'@export
getFlowStats <- function(modelOutput,
  zero.val=0.001,cens.level=.005,min.obs=0.01) {
  # Function developed by William Farmer, 09 June 2015

  # For each site, calculate flow statistics.
  for (i in 1:length(modelOutput)) {
    # Specific site flows
    Dates <- modelOutput[[i]]$date
    Obs <- modelOutput[[i]]$obs
    Est <- modelOutput[[i]]$est
    # Control zeros, if needed
    ndx0 <- Obs<cens.level
    ndx1 <- Obs>=cens.level&Obs<min.obs
    Obs[ndx0] <- zero.val
    Obs[ndx1] <- min.obs
    ndx0 <- Est<cens.level
    ndx1 <- Est>=cens.level&Est<min.obs
    Est[ndx0] <- zero.val
    Est[ndx1] <- min.obs
    # Only look at overlapping records
    NA.ndx <- is.na(Obs+Est+as.double(substr(Dates,1,4)))
    Est <- Est[!NA.ndx]
    Obs <- Obs[!NA.ndx]
    Dates <- Dates[!NA.ndx]
    # Calculate statistics
    if (length(Obs)==0) {next}
    a <- calcFlowStats(Obs,Dates)
    b <- calcFlowStats(Est,Dates)
    if (i == 1) {
      flowStats.Obs <- flowStats.Est <-
        matrix(NA,nrow=length(modelOutput),ncol=length(a))
    }
    flowStats.Obs[i,] <- a
    flowStats.Est[i,] <- b
  }
  flowStats.Obs <- data.frame(flowStats.Obs)
  flowStats.Est <- data.frame(flowStats.Est)
  names(flowStats.Obs) <- names(flowStats.Est) <- names(a)

  result <- list(Obs=flowStats.Obs,Est=flowStats.Est)
  return(result)
}
