

# krigeLUD ----
#' Apply kriging of drainage area ratio
#'
#' @description
#' The function \code{krigeLUD} uses daily variograms of the logarithm of the
#' unit discharge to predict streamflows at the target locations.
#'
#' @param index.gages A numeric vector of the streamgage IDs to be used as
#' potential index sites.
#' @param index.baschar A data frame with the \code{LAT_GAGE_UTM},
#' \code{LNG_GAGE_UTM} and \code{DRAIN_SQKM} for each item of
#' \code{index.gages}.
#' @param index.obs A zoo object of the observed streamflows for each
#' index gage.
#' @param target.gages A numeric vector of the streamgage IDs to be used as
#' target sites.
#' @param target.baschar Identical to \code{index.baschar} with respect to
#' \code{target.gages}.
#' @param target.obs A zoo object of the observed streamflows for each
#' target gage.
#' @param zero.val (optional) A value to replace zeros.  The deafult is
#' \code{NA}.
#' @param FixNug (optional) A logical value indicating if the nugget value
#' should be fixed.  The deafult is \code{FALSE}.
#' @param FixKap (optional) A logical value indicating if the kappa value
#' should be fixed.  The deafult is \code{TRUE}.
#' @param numbins (optional) The number of bins to be used in each variogram.
#' The deafult is \code{10}.
#' @param distperc (optional) The percentile of the distances that should
#' constrain the range of the empirical variograms.  The defaul is \code{1}.
#' @param CovMod (optional) A character string indicating which theoretical
#' variogram model should be used.  The deafult is the \code{'spherical'} model.
#'
#'
#' @details
#' Lorem ipsum...
#'
#' @return A list of two elements:
#' \item{OKLUD}{\itemize{A list with each element corresponding to \code{target.gages},
#' with each element of that list containg a data frame of
#' \item{date}{The date of the row.}
#' \item{est}{The estimated streamflow values for each day.}
#' \item{est.LUD}{The raw kriged logarithm of the unit discharge for each day.}
#' \item{est.LUD.var}{The variance of the raw kriged logarithm of the unit
#' discharge for each day.}
#' \item{extrap}{A logical indicating if the estimate was extrapolated using
#' nearest-neighbor drainage area ratios.}
#' \item{ratio}{If extrapolation was used, the ratio between target and index
#' locations.}
#' \item{indexflow}{If extrapolation was used, the index streamflow.}
#' \item{index}{If extrapolation was used, the index ID.}
#' \item{metric}{If extrapolation was used, the similarity metric between the
#' index and target.}
#' \item{obs}{The observed streamflow.}}}
#' \item{varPar}{A data frame of daily variogram parameters.  These include
#' covariance parameters, nugget values, kappa values, the variogram value,
#' the covariance model, the maximum distance, the distance percentile and
#' the number of bins.}
#'
#'@export
#'@import zoo
#'@import geoR
krigeLUD <- function(index.gages,index.baschar,index.obs,
  target.gages,target.baschar,target.obs,zero.val=NA,
  FixNug=F,FixKap=T,numbins=10,distperc=1,CovMod="spherical") {
  # Based on code orginally developed by William Farmer, 31 July 2014
  # Revised by William Farmer, 05 June 2015

  # @importFrom zoo index as.Date
  # @importFrom geoR variog variofit krige.conv

  # Assume all TS are the same length... (with NA for missing)
  TargetDates <- index(target.obs[[1]])

  # Reformat observations in to matrices
  target.obs.temp <- matrix(unlist(target.obs),ncol=length(target.obs))
  index.obs.temp <- matrix(unlist(index.obs),ncol=length(index.obs))
  target.obs.temp[target.obs.temp==0] <- zero.val
  index.obs.temp[index.obs.temp==0] <- zero.val

  # Empty matrices for saving
  varPar <- matrix(NA,nrow=length(TargetDates),ncol=9)
  result <- list()
  for (i in 1:length(target.obs)) {
    result[[i]] <- matrix(NA,ncol=10,nrow=length(TargetDates))
    # {date,est,est.raw,est.raw.var,extrap,ratio,indexflow,index,metric,obs}
    result[[i]] <- data.frame(result[[i]])
    names(result[[i]]) <- c('date','est','est.LUD','est.LUD.var','extrap',
      'ratio','indexflow','index','metric','obs')
    class(result[[i]]$extrap) <- 'logical'
    result[[i]]$date <- as.Date(result[[i]]$date)
  }

  IndexLatLon <- cbind(index.baschar$LAT_GAGE_UTM,index.baschar$LNG_GAGE_UTM)
  TargetLatLon <- cbind(target.baschar$LAT_GAGE_UTM,target.baschar$LNG_GAGE_UTM)
  distances <- dist(IndexLatLon,diag=T,upper=T)
  if (distperc==1){
    maxrange <- summary(distances)[6]
  } else {
    maxrange <- quantile(distances,probs=distperc)
  }
  # Daily Model
  for (j in 1:length(TargetDates)) {
    #cat(paste(j,"-",Sys.time(),"\n"))
    SubSetLUD <- log(index.obs.temp[j,]/index.baschar$DRAIN_SQKM)
    SubSetLUD[is.infinite(SubSetLUD)] <- NA
    NDX <- !is.na(SubSetLUD)
    if (sum(NDX)==1|sum(NDX)==2|sum(NDX)==3) {
      for (i in 1:length(result)) {
        result[[i]][j,1] <- TargetDates[j]
        result[[i]][j,2] <- exp(SubSetLUD[NDX][1])*target.baschar$DRAIN_SQKM[i]
        result[[i]][j,3] <- NA
        result[[i]][j,4] <- NA
        result[[i]][j,5] <- TRUE
        result[[i]][j,10] <- target.obs.temp[j,i]
      }
      next
    } else if (sum(NDX)==0) {
      for (i in 1:length(result)) {
        result[[i]][j,1] <- TargetDates[j]
        result[[i]][j,2] <- NA
        result[[i]][j,3] <- NA
        result[[i]][j,4] <- NA
        result[[i]][j,5] <- FALSE
        result[[i]][j,10] <- target.obs.temp[j,i]
      }
      next
    }
    gage.variogram <- variog(coords=IndexLatLon[NDX,],data=SubSetLUD[NDX],
      breaks = seq(0,maxrange,len=numbins),messages=F)
    gage.variofit <- variofit(gage.variogram, cov.model = CovMod,
      fix.nugget=FixNug,fix.kappa=FixKap,messages=F)
    if (gage.variofit$nugget<0) {gage.variofit$nugget<-0}
    # correction for negative nuggets because
    #   "micro.scale must be in the interval [0, nugget]"
    estLUD<-krige.conv(coords=IndexLatLon[NDX,],data=SubSetLUD[NDX],
      locations=TargetLatLon,krige=krige.control(obj.m=gage.variofit),
      output=output.control(messages=F))
    varPar[j,] <- c(gage.variofit$cov.pars,gage.variofit$nugget,
      gage.variofit$kappa,gage.variofit$value,gage.variofit$cov.model,
      gage.variofit$max.dist,distperc,numbins)

    # {date,est,est.LUD,est.LUD.var,extrap,ratio,indexflow,index,metric,obs}
    for (i in 1:length(result)) {
      result[[i]][j,1] <- TargetDates[j]
      result[[i]][j,2] <- exp(estLUD$predict[i])*target.baschar$DRAIN_SQKM[i]
      result[[i]][j,3] <- estLUD$predict[i]
      result[[i]][j,4] <- estLUD$krige.var[i]
      if (is.na(estLUD$predict[i])) {
        result[[i]][j,5] <- TRUE
      } else {
        result[[i]][j,5] <- FALSE
      }
      result[[i]][j,10] <- target.obs.temp[j,i]
    }
  }

  # Fill extrapolations with NN-DAR
  NN.ranking <- index.NN(index.gages,index.baschar,target.gages,target.baschar)
  NNDAR <- estDAR(index.network=NN.ranking,index.baschar,
    target.baschar,index.obs,target.obs)
  for (i in 1:length(result)) {
    result[[i]][,10] <- NNDAR[[i]]$obs
    if (sum(result[[i]][,5])==0) {next}
    ndx <- which(result[[i]][,5])
    result[[i]][ndx,c(2,6,7,8,9)] <- NNDAR[[i]][ndx,c(2,3,4,5,6)]
  }
  names(result) <- names(NNDAR)

  OutPut <- list(OKLUD=result,varPar=varPar)
  return(OutPut)
}
