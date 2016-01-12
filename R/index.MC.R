
# index.MC ----
#' Rank the best index gages for each target site.
#'
#' @description
#' The function \code{index.MC} returns the ranking of potential index sites,
#' by map correlation, for each target location.
#'
#' @param index.gages A numeric vector of the streamgage IDs to be used as
#' potential index sites.
#' @param index.obs A zoo object of the observed streamflows for each
#' index gage.
#' @param index.baschar A data frame with the \code{LAT_GAGE_UTM},
#' \code{LNG_GAGE_UTM} and \code{DRAIN_SQKM} for each item of
#' \code{index.gages}.
#' @param target.gages A numeric vector of the streamgage IDs to be used as
#' target sites.
#' @param target.obs A zoo object of the observed streamflows for each
#' target gage.
#' @param target.baschar Identical to \code{index.baschar} with respect to
#' \code{target.gages}.
#' @param method A character string of either \code{"pearson"},
#' \code{"spearman"} and \code{"kendall"}, specifying which correlation
#' method to use.
#'
#' @details
#' Lorem ipsum...
#'
#' @return The function returns a list as output.  This list contains:
#' \item{Cor.Obs}{A matrix of the observed correlation between index and target
#' sites. The columns are in the order of \code{target.gages} and
#' the rows are in the order of \code{index.gages}.}
#' \item{Cor.Netw}{A matrix of the observed inter-correlation between index
#' sites. The columns are in the order of \code{target.gages} and
#' the rows are in the order of \code{index.gages}.}
#' \item{Cor.Est}{A matrix of the estimated correlation between index and target
#' sites. The columns are in the order of \code{target.gages} and
#' the rows are in the order of \code{index.gages}.}
#' \item{varPar}{A matrix of the variogram parameters (covariance parameters
#' and nugget) for the variogram built on each index gage.}
#' \item{Rank}{A matrix of the index IDs for each target gage.
#' The columns are in the order of \code{target.gages}.}
#' \item{index}{The input \code{index.gages}.}
#' \item{target}{The input \code{target.gages}.}
#'@export
index.MC <- function(
  index.gages,index.obs,index.baschar,
  target.gages,target.obs,target.baschar,
  method=c("pearson", "kendall", "spearman")) {
  # Based on code orginally developed by William Farmer, 02 February 2014
  # Revised by William Farmer, 03 June 2015

  # Reformat observations in to matrices
  target.obs.temp <- matrix(unlist(target.obs),ncol=length(target.obs))
  index.obs.temp <- matrix(unlist(index.obs),ncol=length(index.obs))

  # Calculate true correlations
  cor.tar.ind <- cor(target.obs.temp,index.obs.temp,
    use="pairwise.complete.obs",method=method)
  cor.ind.ind <- cor(index.obs.temp,index.obs.temp,
    use="pairwise.complete.obs",method=method)
  diag(cor.ind.ind) <- NA # knock-down self-correlation

  # Setup variograms
  target.Coords <- cbind(target.baschar$LAT_GAGE_UTM,
    target.baschar$LNG_GAGE_UTM)
  index.Coords <- cbind(index.baschar$LAT_GAGE_UTM,
    index.baschar$LNG_GAGE_UTM)
  distances <- dist(index.Coords,diag=T,upper=T)
  maxrange <- summary(distances)[6]
  numbins <- 10
  EstCorr <- matrix(NA,ncol=length(index.gages),nrow=length(target.gages))
  target.Coords <- cbind(target.baschar$LAT_GAGE_UTM,
    target.baschar$LNG_GAGE_UTM)
  varPar <- matrix(NA,ncol=3,nrow=length(index.gages))

  # For each index, estimate the variogram and predict target correlations.
  for (i in 1:length(index.gages)) {
    iData <- cor.ind.ind[i,]
    NDX <- !is.na(iData)
    iData <- iData[NDX]
    iCoords <- index.Coords[NDX,]
    gage.variogram <- variog(coords=iCoords,data=iData,
      breaks = seq(0,maxrange,len=numbins),
      messages=F)
    gage.variofit <- variofit(gage.variogram, cov.model = "spherical",
      fix.nugget=TRUE,messages=F)
    varPar[i,] <- c(gage.variofit$cov.pars, gage.variofit$nugget)
    estcorrval<-krige.conv(coords=iCoords,data=iData,locations=target.Coords,
      krige=krige.control(obj.m=gage.variofit),
      output=output.control(messages=F))
    EstCorr[,i] <- estcorrval$predict
  }

  Rank.Mat <- matrix(NA,ncol=length(target.gages),nrow=length(index.gages))
  for (i in 1:length(target.gages)) {
    # Rank index sites for each target
    sorted <- sort.int(EstCorr[i,],decreasing=TRUE,index.return=TRUE)
    Rank.Mat[,i] <- index.gages[sorted$ix]
  }

  result <- list(Cor.Obs=t(cor.tar.ind),Cor.Netw=cor.ind.ind,Cor.Est=t(EstCorr),
    varPar = varPar, Rank = Rank.Mat, index=index.gages,target=target.gages)
  return(result)

}
