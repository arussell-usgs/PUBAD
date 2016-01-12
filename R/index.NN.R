#index.NN----
#' Rank the best index gages for each target site.
#'
#' @description
#' The function \code{index.NN} returns the ranking of potential index sites,
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
index.NN <- function(index.gages,index.baschar,target.gages,target.baschar) {
  # Function developed by William Farmer, 03 June 2015

  Dist.Mat <- matrix(NA,ncol=length(target.gages),nrow=length(index.gages))
  DAR.Mat <- matrix(NA,ncol=length(target.gages),nrow=length(index.gages))
  Rank.Mat <- matrix(NA,ncol=length(target.gages),nrow=length(index.gages))

  for (i in 1:length(target.gages)) {
    # Compute distances between sites
    Dist.Mat[,i] <- sqrt(
      (index.baschar$LAT_GAGE_UTM-target.baschar$LAT_GAGE_UTM[i])^2 +
        (index.baschar$LNG_GAGE_UTM-target.baschar$LNG_GAGE_UTM[i])^2
    )

    # Calculate DAR between sites
    DAR.Mat[,i] <- target.baschar$DRAIN_SQKM[i]/index.baschar$DRAIN_SQKM

    # Rank index sites for each target
    sorted <- sort.int(Dist.Mat[,i],decreasing=FALSE,index.return=TRUE)
    Rank.Mat[,i] <- index.gages[sorted$ix]
  }
  result <- list(Dist=Dist.Mat,DAR=DAR.Mat,Rank=Rank.Mat,
    index=index.gages,target=target.gages)
  return(result)
}
