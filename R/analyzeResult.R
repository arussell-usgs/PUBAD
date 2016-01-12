#' Analyze performance of estimated time series.
#'
#' @description
#' The function \code{analyzeResult} conducts a first level analysis of
#' estimated time series.
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
#' @param ErrLevs (optional) A vector of percent error values that should be
#' used to construct the cumulative distribution of absolute percent errors.
#' The default is \code{seq(0.01,1.5,0.01)}.
#' @param bins (optional) A vector of the bin limits that should be used to
#' characterize the the percent errors along the flow duration curve.  The
#' default is \code{seq(0,1,0.01)}.
#' @param FlowStat (optional) A logical indicating if the streamflow
#' statistics should be evaluated.  This will increase run time.  The default
#' is \code{FALSE}.
#'
#' @details
#' The single-value performance metric calcuated include: The Nash-Sutcliffe
#' efficieny of daily streamflow, the Nash-Sutcliffe efficiency of the
#' logarithms of the daily streamflow, the root-mean-square error statistics
#' from \link{\code{rmse.like}} for natural and logarithm streamflows, the
#' average percent errors of the natural and logarithm streamflows and the
#' Pearson and Spearman correlations between observed and simulated streamflow.
#'
#' @return A list of several elements:
#' \item{PerfMat}{A data frame of single-value performance metrics for each
#' target location.  See details.}
#' \item{CumDistErr}{A list of two elements:
#' \item{Levels}{The input \code{ErrLevs}.}
#' \item{CumFreq}{A matrix summarizing the quantiles of the absolute percent
#' errors at each site.}}
#' \item{ErrAlongFDC}{A list of several elements:
#' \item{bins}{The input \code{bins}.}
#' \item{centers}{The plotting centers of each bin.}
#' \item{meanErr}{The mean percent error for each bin for each target location.}
#' \item{medianErr}{The median percent error for each bin for each target
#' location.}}
#'@export
analyzeResult <- function(modelOutput,
  zero.val=0.001,cens.level=.005,min.obs=0.01,
  ErrLevs=seq(0.01,1.5,0.01),
  bins=seq(0,1,0.01),FlowStat=F) {
  # Function developed by William Farmer, 08 June 2015
  # 13 October 2015: Corrected to screen for WYs. WHF.

  # Water Years
  Dates <- modelOutput[[1]]$date
  Year <- as.numeric(substr(Dates,1,4))
  Month <- as.numeric(substr(Dates,6,7))
  WY <- ifelse(Month>9,Year+1,Year)
  FullYears <- cbind(as.numeric(names(table(WY))),table(WY))

  # Clean data, screen for complete water years
  Obs <- Est <- matrix(NA,ncol=length(modelOutput),nrow=nrow(modelOutput[[1]]))
  for (i in 1:length(modelOutput)) {
    iObs <- modelOutput[[i]]$obs
    iEst <- modelOutput[[i]]$est
    nandx <- which(!is.na(rowSums(cbind(iObs,iEst))))
    iWY <- WY[nandx]
    YD <- cbind(as.numeric(names(table(iWY))),table(iWY))
    YD <- cbind(YD,FullYears[match(YD[,1],FullYears[,1]),2])
    iWYs <- as.numeric(YD[YD[,2]==YD[,3],1])
    wyndx <- which(is.element(WY,iWYs))
    Obs[wyndx,i] <- iObs[wyndx]
    Est[wyndx,i] <- iEst[wyndx]
  }
  ndx0 <- Obs<cens.level
  ndx1 <- Obs>=cens.level&Obs<min.obs
  Obs[ndx0] <- zero.val
  Obs[ndx1] <- min.obs
  ndx0 <- Est<cens.level
  ndx1 <- Est>=cens.level&Est<min.obs
  Est[ndx0] <- zero.val
  Est[ndx1] <- min.obs

  # Calculate basic performance metrics
  NSE <- nse(Obs,Est)
  NSEL <- nse(log(Obs),log(Est))
  RMSE <- rmse.like(Obs,Est)
  RMSEL <- rmse.like(log(Obs),log(Est))
  PErr <- percent.error(Obs,Est)
  PErrL <- percent.error(log(Obs),log(Est))
  Corr <- obs.sim.corr(Obs,Est,methods=c("pearson","spearman"))
  Perf <- data.frame(NSE,NSEL,RMSE,RMSEL,PErr,PErrL,Corr)
  names(Perf) <- c("nse","nsel","rmse","rmsne","nrmse","cvrmse","rmsel",
    "rmsnel","nrmsel","cvrmsel","perr","perrl","cor.p","cor.s")

  # Cumulative distribution of absolute percent errors
  dailyPErr <- (Est-Obs)/Obs
  DisErr <- matrix(NA,nrow=ncol(Obs),ncol=length(ErrLevs))
  for (i in 1:length(ErrLevs)) {
    DisErr[,i] <- colMeans(abs(dailyPErr)<=ErrLevs[i],na.rm=T)
  }

  # Percent error along the flow duration curve
  Probs <- matrix(NA,ncol=ncol(Obs),nrow=nrow(Obs))
  for (i in 1:ncol(Probs)) {
    ndx <- which(!is.na(Obs[,i]))
    Probs[ndx,i] <- rank(Obs[ndx,i])/(length(ndx)+1)
  }
  Binned.PErr.mean <- matrix(NA,nrow=ncol(Obs),ncol=length(bins)-1)
  Binned.PErr.median <- matrix(NA,nrow=ncol(Obs),ncol=length(bins)-1)
  centers <- rep(NA,length(bins)-1)
  for (i in 1:ncol(Obs)) {
    for (j in 1:(length(bins)-1)) {
      ndx <- which(Probs[,i]<=bins[j+1]&Probs[,i]>bins[j])
      Binned.PErr.mean[i,j] <- mean(dailyPErr[ndx,i],na.rm=T)
      Binned.PErr.median[i,j] <- median(dailyPErr[ndx,i],na.rm=T)
    }
  }

  if (FlowStat) {
    # Get flow statistics
    FlowStats <- getFlowStats(modelOutput, zero.val=zero.val,
      cens.level=cens.level, min.obs=min.obs)
    FlowStats.PErr <- (FlowStats$Est-FlowStats$Obs)/FlowStats$Obs

    # Compile output
    result <- list(PerfMat=Perf,
      CumDistErr=list(Levels=ErrLevs,CumFreq=DisErr),
      ErrAlongFDC=list(bins=bins,centers=centers,
        meanErr=Binned.PErr.mean,medianErr=Binned.PErr.median),
      FlowStats=list(Value=FlowStats,PerErr=FlowStats.PErr))
    return(result)
  }

  # Compile output
  result <- list(PerfMat=Perf,
    CumDistErr=list(Levels=ErrLevs,CumFreq=DisErr),
    ErrAlongFDC=list(bins=bins,centers=centers,
      meanErr=Binned.PErr.mean,medianErr=Binned.PErr.median))
  return(result)
}
