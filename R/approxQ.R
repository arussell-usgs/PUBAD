#' Interpolate and Extrapolate Streamflow Probabilities along Flow Duration
#' Curve
#'
#' @description The function \code{approxQ} interpolates and extrapolates a time
#'   series of probabilities along a given set of quantiles from a flow duration
#'   curve.
#'
#' @param pts A numerical vector containing a time series of streamflow
#'   probabilities.  The probabilities can be either nonexceedance or exceedance
#'   probabilities, as long as they correspond to the probabilities in
#'   \code{fdc$p} described below.
#' @param fdc A data frame of two numerical variables.  The first, \code{p},
#'   contains streamflow probabilities corresponding to the streamflow values in
#'   the second, \code{q}.  The probabilities can be either nonexceedance or
#'   exceedance probabilities, as long as they correspond to the probabilities
#'   in \code{pts} described above.  Non-positive values of \code{q} will be
#'   removed because of the conversion to common logarithms.
#'
#' @details The interpolation is accomplished by converting streamflows to the
#' common logarithms of streamflow and converting the probabilites to standard
#' normal quantlies.  The \code{approx} is used to then linearly
#' interpolate the standard normal quantiles against the common logarithms of
#' streamflow.
#'
#' Probabilites in \code{pts} that are beyond the maximum value in \code{fdc$p}
#' are extrapolated by extending the line between the two greatest quantiles of
#' the flow duration curve.
#'
#' Similarly, probabilites in \code{pts} that are below the minimum value in
#' \code{fdc$p} are extrapolated by extending the line between the two smallest
#' quantiles of the flow duration curve.
#'
#' @return A list of three elements: \item{pts}{The input time series of
#'   streamflow probabilites.} \item{fdc}{The input flow duration curve, with
#'   addition of the common logarithm of streamflow (\code{ql}) and the standard
#'   nomral quantiles of the probabilities (\code{z}).} \item{targetFlow}{The
#'   resulting conversion of \code{pts} into streamflows by interpolating and
#'   extrapolating along the given flow duration curve.}
#'
#' @examples
#' set.seed(1)
#' pts <- runif(100)
#' fdc <- data.frame(
#'   p = seq(from = 0.01, to = 0.99, by = 0.01)
#' )
#' fdc$q <- qnorm(fdc$p, mean = 5, sd = 2)
#'
#' result <- approxQ(pts = pts, fdc = fdc)
#' @export
approxQ <- function(pts, fdc) {
  # 14 September 2016
  # William Farmer, wfarmer@usgs.gov
  #
  # Convert FDC into log10 and remove NAs and Infs
  fdc$ql <- log10(fdc$q)
  fdc$z <- qnorm(fdc$p)
  fdcIn <- fdc
  if (sum(is.infinite(fdc$ql)) > 0) {
    fdc[which(is.infinite(fdc$ql)), c('q', 'ql')] <- NA
  }
  if (sum(is.na(fdc$ql)) > 0) {
    fdc <- fdc[-which(is.na(fdc$q)), ]
  }
  # interpolate all values
  targetFlow <- 10^(approx(fdc$z, fdc$ql, qnorm(pts))$y)
  # extrapolate high-end p-values
  ndx <- which(pts > max(fdc$p) & is.na(targetFlow))
  if (length(ndx) > 0) {
    sorted <-
      sort.int(fdc$p, decreasing = TRUE, index.return = TRUE)$ix[c(1, 2)]
    intrpQ <- fdc$ql[sorted]
    intrpP <- fdc$z[sorted]
    m <- (intrpQ[1] - intrpQ[2]) / (intrpP[1] - intrpP[2])
    b <- intrpQ[1] - m * intrpP[1]
    targetFlow[ndx] <- 10 ^ (m * qnorm(pts[ndx]) + b)
  }
  # extrapolate low-end p-values
  ndx <- which(pts < min(fdc$p) & is.na(targetFlow))
  if (length(ndx) > 0) {
    sorted <-
      sort.int(fdc$p, decreasing = FALSE, index.return = TRUE)$ix[c(1, 2)]
    intrpQ <- fdc$ql[sorted]
    intrpP <- fdc$z[sorted]
    m <- (intrpQ[1] - intrpQ[2]) / (intrpP[1] - intrpP[2])
    b <- intrpQ[1] - m * intrpP[1]
    targetFlow[ndx] <- 10 ^ (m * qnorm(pts[ndx]) + b)
  }
  # output results and inputs
  output <- list(fdc = fdcIn, pts = pts, targetFlow = targetFlow)
  return(output)
}
