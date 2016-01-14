#' Calculate streamflow statistics.
#'
#' @description
#' The function \code{calcFlowStats} calculates a wide array of streamflow
#' statistics for a given time series.
#'
#' @param Flow A vector of streamflow values.
#' @param Dates A corresponding vector of dates.
#'
#' @details
#' Calculated statistics include: The daily median, the daily mean (FDSS 1),
#' the daily coefficient of variation, the annual coefficient of variation,
#' several streamflow percentiles (\code{percs=c(0.01,0.05,0.1,0.25,0.75,
#' 0.9,0.95,0.99)}), the 90th percentile of annual maximums, the L-CV (FDSS 2),
#' the L-skew (FDSS 3), the L-kurtosis (FDSS 4), the lag-one autocorrelation (
#' FDSS 5), the amplitude and phase (FDSS 6 and 7) of the sinusoidal seasonal
#'  signal, and the 50th and 10th percentiles of the annual minimum seven-day
#'  average streamflow.
#'
#' @return A named vector of streamflow statistics.
#'@export
#'@import zoo
calcFlowStats <- function(Flow,Dates) {
  # Based on analysis from SIR 2014-5231 (circa 16 June 2014)
  # FDSS calculations based on code from Stacey Archfield
  # Function redeveloped by William Farmer, 09 June 2015

  # @importFrom zoo as.Date rollmean aggregate.zoo

  # Find and index complete water years
  Wyrs <- min(as.double(substr(Dates,1,4))):max(as.double(substr(Dates,1,4)))
  comp.Wyrs <- logical(length(Wyrs))
  for (yr in 1:length(Wyrs)) {
    beg.date <- as.Date(paste(Wyrs[yr]-1,10,01,sep="-"))
    end.date <- as.Date(paste(Wyrs[yr],09,30,sep="-"))
    date.inds <- (Dates>=beg.date & Dates<=end.date)
    nDays.Wyr.site <- sum(!is.na(Flow[date.inds]))
    nDays.Wyr <- end.date - beg.date + 1
    if (nDays.Wyr.site==nDays.Wyr) {
      comp.Wyrs[yr] <- T
    } else {
      comp.Wyrs[yr] <- F
    }
  }
  Wyr.site <- ifelse(as.double(substr(Dates,6,7))>9,
    as.double(substr(Dates,1,4))+1,as.double(substr(Dates,1,4)))
  ndxCompWyr <- which(is.element(Wyr.site,Wyrs[comp.Wyrs]))

  # Find  and index complete climate years
  Cyrs <- min(as.double(substr(Dates,1,4))):max(as.double(substr(Dates,1,4)))
  comp.Cyrs <- logical(length(Cyrs))
  for (yr in 1:length(Cyrs)) {
    beg.date <- as.Date(paste(Cyrs[yr],04,01,sep="-"))
    end.date <- as.Date(paste(Cyrs[yr]+1,03,31,sep="-"))
    date.inds <- (Dates>=beg.date & Dates<=end.date)
    nDays.Cyr.site <- sum(!is.na(Flow[date.inds]))
    nDays.Cyr <- end.date - beg.date + 1
    #Water year is complete if nvals.pWyr matches glb.Wyr.ndays
    if (nDays.Cyr.site==nDays.Cyr) {
      comp.Cyrs[yr] <- T
    } else {
      comp.Cyrs[yr] <- F
    }
  }
  Cyr.site <- ifelse(as.double(substr(Dates,6,7))<4,
    as.double(substr(Dates,1,4))-1,as.double(substr(Dates,1,4)))
  ndxCompCyr <- which(is.element(Cyr.site,Cyrs[comp.Cyrs]))

  flowStats <- NA
  nms <- NA
  # Calculate statistics on complete water years
  wyQ <- Flow[ndxCompWyr]
  wyD <- Wyr.site[ndxCompWyr]
  wyDf <- Dates[ndxCompWyr]
  ## Daily Median
  i <- 1
  flowStats[i] <- median(wyQ)
  nms[i] <- "dailymean"
  ## Daily Mean (FDSS 1)
  i <- i + 1
  flowStats[i] <- mean(wyQ)
  nms[i] <- "dailymedian"
  ## Daily CV
  i <- i + 1
  flowStats[i] <- sd(wyQ)/mean(wyQ)
  nms[i] <- "dailyCV"
  ## Annual CV
  i <- i + 1
  flowStats[i] <- sd(aggregate(wyQ, list(wyD), mean)$x)/
    mean(aggregate(wyQ, list(wyD), mean)$x)
  nms[i] <- "annualCV"
  ## Percentiles
  i <- i + 1
  percs <- c(0.01,0.05,0.1,0.25,0.75,0.9,0.95,0.99)
  flowStats[i:(i+length(percs)-1)] <- quantile(wyQ, probs = percs)
  nms[i:(i+length(percs)-1)] <- paste0("percs",round(percs,digits=2))
  ## 90th percentile of annual maximum
  i <- i + length(percs)
  flowStats[i] <- quantile(aggregate(wyQ, list(wyD), max)$x,probs=0.9)
  nms[i] <- "annmax90"
  ## L-CV (FDSS 2)
  complmom<-lmom.ub(wyQ)
  i <- i + 1
  flowStats[i] <- complmom$LCV
  nms[i] <- "LCV"
  ## L-Skew (FDSS 3)
  i <- i + 1
  flowStats[i] <- complmom$TAU3
  nms[i] <- "Lskew"
  ## L-kurtosis (FDSS 4)
  i <- i + 1
  flowStats[i] <- complmom$TAU4
  nms[i] <- "Lkurt"
  ## Auto-regressive Lag-1 Correlation (FDSS 5)
  ### Deseasonalize
  decimal.month <- as.double(substr(wyDf,1,4)) +
    ((as.double(substr(wyDf,6,7)))/12)
  monthlymean <- aggregate(wyQ,list(decimal.month),mean)
  wyQ.ds <- wyQ
  wyQ.ds[] <- NA
  for (j in 1:nrow(monthlymean)) {
    NDX<-which(decimal.month==monthlymean$Group.1[j])
    wyQ.ds[NDX]<- (wyQ[NDX]-monthlymean$x[j])
  }
  #### Estimate AR1
  armdl<-ar(scale(wyQ.ds, center = TRUE, scale = TRUE),
    aic = FALSE, order.max = 1, method="yule-walker")
  i <- i + 1
  flowStats[i] <- armdl$ar
  nms[i] <- "AR1"
  ## Seasonal amplitude (FDSS 6)
  wyQ.s <- scale(wyQ, center = TRUE, scale = TRUE)
  decimal_year<- as.double(substr(wyDf,1,4)) +
    ((as.POSIXlt(wyDf,"yyyy-mm-dd")$yday+1)/365.25)
  seasonfit<-lm(wyQ.s~cos(2*pi*decimal_year)+sin(2*pi*decimal_year))
  seasonA<-as.vector(seasonfit$coefficients[2])
  seasonB<-as.vector(seasonfit$coefficients[3])
  i <- i + 1
  flowStats[i] <- sqrt((seasonA^2)+(seasonB^2))
  nms[i] <- "seasAmp"
  ## Seasonal phase (FDSS 7)
  i <- i + 1
  flowStats[i] <- atan((-seasonB)/seasonA)
  nms[i] <- "seasPhs"

  # Calculate statistics on complete climate years
  cyQ <- Flow[ndxCompCyr]
  cyD <- Cyr.site[ndxCompCyr]
  ## 50th percentile of seven-day average annual minimum
  i <- i + 1
  flowStats[i] <- quantile(
    aggregate(
      rollmean(cyQ, 7, align = "right",na.pad = TRUE),
      list(cyD), min, na.rm = TRUE)$x
    ,probs=0.5)
  nms[i] <- "l750"
  ## 10th percentile of seven-day average annual minimum
  i <- i + 1
  flowStats[i] <- quantile(
    aggregate(
      rollmean(cyQ, 7, align = "right",na.pad = TRUE),
      list(cyD), min, na.rm = TRUE)$x
    ,probs=0.1)
  nms[i] <- "l710"

  names(flowStats) <- nms
  return(flowStats)
}
