# estDAR ----
#' Apply drianage-area ratio for a particular target and index set
#'
#' @description
#' The function \code{estDAR} uses the ranked index network defined by
#' \code{index.network} to apply the drainage area ratio to predict streamflow
#' time series at each target location.
#'
#' @param index.network A numeric vector of the streamgage IDs to be used as
#' potential index sites.
#' @param index.baschar A data frame with the \code{LAT_GAGE_UTM},
#' \code{LNG_GAGE_UTM} and \code{DRAIN_SQKM} for each item of
#' \code{index.network$index}.
#' @param target.baschar Identical to \code{index.baschar} with respect to
#' \code{index.network$target}.
#' @param index.obs A zoo object of the observed streamflows for each
#' index gage.
#' @param target.obs A zoo object of the observed streamflows for each
#' target gage.
#' @param DAR.low The lower limit of the DAR that will be permitted for valid
#' application of DAR.  The default is \code{-Inf}.
#' @param DAR.high The upper limit of the DAR that will be permitted for valid
#' application of DAR.  The default is \code{Inf}.
#'
#' @details
#' Lorem ipsum...
#'
#' @return As output, the function returns a list of the same length as
#' \code{index.network$target}.  Each element represents a unique target
#' location.  Each element contains a data.frame consisting of
#' \item{date}{The date of the row.}
#' \item{est}{The estimated streamflow for each day.}
#' \item{ratio}{The ratio of the index and target drainage areas for each day.
#' This only changes when an alternative index is required.}
#' \item{indexflow}{The index streamflow for each day.}
#' \item{index}{The index ID used on each day.}
#' \item{metric}{The similarlity metric between the index and target.  This
#' depends on the metric using in the \code{index.network}; currently
#' functional for nearest-neighbor and map-correlation.}
#' \item{obs}{The observed streamflow at the target location on each day.
#' Used for validation.}
#'@export
#'@import zoo
estDAR <- function(index.network,index.baschar,
  target.baschar,index.obs,target.obs,
  DAR.low=-Inf,DAR.high=Inf) {
  # Based on code orginally developed by William Farmer, 31 July 2014
  # Revised by William Farmer, 03 June 2015

  # @importFrom zoo as.Date

  # Setup matrices
  result <- list()

  for (i in 1:length(index.network$target)) {
    # Get Ungaged Info
    Site.ID <- index.network$target[i]
    Site.DA <- target.baschar$DRAIN_SQKM[i]

    # Potential indices
    indices <- index.network$Rank[,i]

    # Cycle Through Dates
    Dates2Est <- index(target.obs[[i]])
    N <- length(Dates2Est)
    Data <- matrix(NA,nrow=N,ncol=7)
    # {date,est,ratio,indexflow,index,metric}
    Data[,7] <- target.obs[[i]]
    j <- 0
    while (sum(is.na(Dates2Est))<N&(j+1)<length(indices)) {
      j <- j + 1
      Index <- indices[j]
      ndx <- match(Index,index.network$index)
      Index.DA <- index.baschar$DRAIN_SQKM[ndx]
      if (Site.DA/Index.DA<DAR.low|Site.DA/Index.DA>DAR.high) {next}
      if (is.element("Dist",names(index.network))) {
        imet <- index.network$Dist[ndx,i]
      } else {
        imet <- index.network$Cor.Est[ndx,i]
      }
      # Select index data
      Index.data <- index.obs[[ndx]][!is.na(index.obs[[ndx]])]
      # Find Matching Dates
      NDX <- which(is.element(index(Index.data),Dates2Est))
      if (length(NDX)==0) {next}
      IndexDates <- index(Index.data)[NDX]
      NDX2 <- which(is.element(Dates2Est,IndexDates))
      TargetDates <- Dates2Est[NDX2]
      IndexFlow <- Index.data[NDX]

      Data[NDX2,1] <- as.Date(IndexDates) # dates (to double-check)
      Data[NDX2,2] <- IndexFlow*Site.DA/Index.DA # Q_DAR
      Data[NDX2,3] <- Site.DA/Index.DA # DAR
      Data[NDX2,4] <- IndexFlow # indexflow
      Data[NDX2,5] <- Index # index ID
      Data[NDX2,6] <- imet # selection metric

      Dates2Est[NDX2] <- NA
      length(NDX)
    }
    Data <- data.frame(Data)
    names(Data) <- c('date','est','ratio','indexflow','index','metric','obs')
    Data$date <- as.character(Data$date)
    class(Data$date) <- 'Date'
    result[[i]] <- Data
    names(result)[[i]] <- Site.ID
  }
  return(result)
}
