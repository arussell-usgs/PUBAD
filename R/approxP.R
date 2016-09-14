#' Interpolate and Extrapolate Streamflows along Flow Duration Curve to Return
#' Probabilities
#'
#' @description The function \code{approxP} interpolates and extrapolates a time
#'   series of streamflows along a given set of quantiles from a flow duration
#'   curve to produce streamflow probabilities.
#'
#' @param qts A numerical vector containing a time series of streamflows.
#'   Probabilities will not be estimated for non-positive values.
#' @param fdc A data frame of two numerical variables.  The first, \code{p},
#'   contains streamflow probabilities corresponding to the streamflow values in
#'   the second, \code{q}.  The probabilities can be either nonexceedance or
#'   exceedance probabilities; the output will be of the same nature.
#'   Non-positive values of \code{q} will be removed because of the conversion
#'   to common logarithms.
#'
#' @details The interpolation is accomplished by converting streamflows to the
#'   common logarithms of streamflow and converting the probabilites to standard
#'   normal quantlies.  The \code{approx} function is used to then linearly
#'   interpolate the standard normal quantiles against the common logarithms of
#'   streamflow.
#'
#'   Streamflows in \code{qts} that are beyond the maximum value in \code{fdc$q}
#'   are extrapolated by extending the line between the two greatest quantiles
#'   of the flow duration curve.
#'
#'   Similarly, streamflows in \code{qts} that are below the minimum value in
#'   \code{fdc$q} are extrapolated by extending the line between the two
#'   smallest quantiles of the flow duration curve.
#'
#' @return A list of three elements: \item{qts}{The input time series of
#'   streamflows.} \item{fdc}{The input flow duration curve, with addition of
#'   the common logarithm of streamflow (\code{ql}) and the standard nomral
#'   quantiles of the probabilities (\code{z}).} \item{targetP}{The resulting
#'   conversion of \code{qts} into streamflow probabilities by interpolating and
#'   extrapolating along the given flow duration curve.}
#'
#' @examples
#' set.seed(1)
#' qts <- qnorm(runif(100), mean = 5, sd = 2)
#' fdc <- data.frame(
#'   p = seq(from = 0.01, to = 0.99, by = 0.01)
#' )
#' fdc$q <- qnorm(fdc$p, mean = 5, sd = 2)
#'
#' result <- approxP(qts = qts, fdc = fdc)
#' @export
approxP <- function(qts, fdc) {
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
  targetP <- pnorm(approx(fdc$ql, fdc$z, log10(qts))$y)
  # extrapolate high-end q-values
  ndx <- which(qts > max(fdc$q) & is.na(targetP))
  if (length(ndx)>0) {
    sorted <-
      sort.int(fdc$q, decreasing = TRUE, index.return = TRUE)$ix[c(1, 2)]
    intrpQ <- fdc$ql[sorted]
    q <- 2
    while (length(unique(intrpQ)) == 1) {
      q <- q + 1
      sorted <- sort.int(fdc$q, decreasing = TRUE,
        index.return = TRUE)$ix[c(1,q)]
      intrpQ <- fdc$ql[sorted]
    }
    intrpP <- fdc$z[sorted]
    m <- (intrpP[1] - intrpP[2]) / (intrpQ[1] - intrpQ[2])
    b <- intrpP[1] - m * intrpQ[1]
    targetP[ndx] <- pnorm(m * log10(qts[ndx]) + b)
  }
  # extrapolate low-end p-values
  ndx <- which(qts < min(fdc$p) & qts > 0 & is.na(targetP))
  if (length(ndx)>0) {
    sorted <-
      sort.int(fdc$q, decreasing = FALSE, index.return = TRUE)$ix[c(1, 2)]
    intrpQ <- fdc$ql[sorted]
    q <- 2
    while (length(unique(intrpQ)) == 1) {
      q <- q + 1
      sorted <- sort.int(fdc$q, decreasing = FALSE,
        index.return = TRUE)$ix[c(1,q)]
      intrpQ <- fdc$ql[sorted]
    }
    intrpP <- fdc$z[sorted]
    m <- (intrpP[1] - intrpP[2]) / (intrpQ[1] - intrpQ[2])
    b <- intrpP[1] - m * intrpQ[1]
    targetP[ndx] <- pnorm(m * log10(qts[ndx]) + b)
  }
  # output results and inputs
  output <- list(fdc = fdcIn, qts = qts, targetP = targetP)
  return(output)
}
