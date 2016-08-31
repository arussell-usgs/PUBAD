#indexNN----
#' Rank the best index gages for each target site.
#'
#' @description
#' The function \code{indexNN} returns the ranking of potential index sites,
#' by geographic proximity, for each target location.
#'
#' @param index.gages A numeric vector of the streamgage IDs to be used as
#' potential index sites.
#' @param index.baschar A data frame with the \code{LAT_GAGE_UTM},
#' \code{LNG_GAGE_UTM} and \code{DRAIN_SQKM} for each item of
#' \code{index.gages}.
#' @param target.gages A numeric vector of the streamgage IDs to be used as
#' target sites.
#' @param target.baschar Identical to \code{index.baschar} with respect to
#' \code{target.gages}.
#'
#' @details
#' Lorem ipsum...
#'
#' @return The function returns a list as output.  This list contains:
#' \item{Dist}{A matrix of the distances between index and
#' target sites.  The columns are in the order of \code{target.gages} and
#' the rows are in the order of \code{index.gages}.}
#' \item{DAR}{A matrix of the drainage area ratios between index and
#' target sites.  The columns are in the order of \code{target.gages} and
#' the rows are in the order of \code{index.gages}.}
#' \item{Rank}{A matrix of the index IDs for each target gage.
#' The columns are in the order of \code{target.gages}.}
#' \item{index}{The input \code{index.gages}.}
#' \item{target}{The input \code{target.gages}.}
#'@export
indexNN <- function(index.gages,index.baschar,target.gages,target.baschar,
                    UTM=T, outlet=T, DAR.fact=0) {

  # Modified by Tom Over and Amy Russell, June 2016

  # Generalized from function indexNN developed by William Farmer, 03 June 2015

  # Incorporates features of indexNN_centroid and indexNN_outlet
  # created by Amy Russell, April 2016, and adds the option of including
  # the DAR in the distance measure.
  # New parameters:
  # UTM -   if true (default), compute distance between gages (outlets) using UTM
  #         locations; otherwise use Vincenty ellipsoid distance
  # outlet - if true (default) and UTM is false, use Vincenty ellipsoid distance
  #          between gages (outlets); if false and UTM is false,
  #          use Vincenty ellipsoid distance between basins centroids
  #          if UTM is true, outlet has no effect.
  # DAR.fact - factor >=0 to weight abs(log10(DAR)) in computation of generalized distance,
  #            (GD) which is computed as GD = Dist.Mat * (1 + DAR.fact*abs(log10(DAR))).
  #            DAR.fact is dimensionless and defaults to zero in which case the DAR ratio
  #            has no effect.

  Dist.Mat <- matrix(NA,ncol=length(target.gages),nrow=length(index.gages))
  DAR.Mat <- matrix(NA,ncol=length(target.gages),nrow=length(index.gages))
  Rank.Mat <- matrix(NA,ncol=length(target.gages),nrow=length(index.gages))

  for (i in 1:length(target.gages)) {
    # Compute distances between sites
    if (UTM) {
      Dist.Mat[,i] <- sqrt(
        (index.baschar$LAT_GAGE_UTM-target.baschar$LAT_GAGE_UTM[i])^2 +
          (index.baschar$LNG_GAGE_UTM-target.baschar$LNG_GAGE_UTM[i])^2)
    } else {

      # Need index gages in matrix for Vincenty formula
      #  Important that longitude is first column and latitude is second
      if (outlet){

        Index.LatLon <- data.frame(index.baschar$LNG_GAGE, index.baschar$LAT_GAGE)
        Index.LatLon.Mat <- as.matrix(Index.LatLon)
        Dist.Mat[,i] <- distVincentyEllipsoid(c(target.baschar$LNG_GAGE[i],
                                                target.baschar$LAT_GAGE[i]),
                                              Index.LatLon.Mat)
      } else {

        Index.LatLon <- data.frame(index.baschar$LONG_CENT, index.baschar$LAT_CENT)
        Index.LatLon.Mat <- as.matrix(Index.LatLon)
        Dist.Mat[,i] <- distVincentyEllipsoid(c(target.baschar$LONG_CENT[i],
                                                target.baschar$LAT_CENT[i]),
                                              Index.LatLon.Mat)
      }
    }

    # Calculate DAR between sites
    DAR.Mat[,i] <- target.baschar$DRAIN_SQKM[i]/index.baschar$DRAIN_SQKM

    # Add DAR effect to Dist.Mat
    Dist.Mat[,i] <- Dist.Mat[,i] * (1 + DAR.fact*abs(log10(DAR.Mat[,i])))

    # Rank index sites for each target
    sorted <- sort.int(Dist.Mat[,i],decreasing=FALSE,index.return=TRUE)
    Rank.Mat[,i] <- index.gages[sorted$ix]
  }
  result <- list(Dist=Dist.Mat,DAR=DAR.Mat,Rank=Rank.Mat,
                 index=index.gages,target=target.gages)
  return(result)
}
