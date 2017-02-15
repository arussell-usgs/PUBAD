# estQPPQ ----
#' Apply QPPQ for a particular target and index set
#'
#' @description
#' The function \code{estQPPQ} uses the ranked index network defined by
#' \code{index.network} to apply the nonlinear spatial interpolation using
#' flow duration curves to predict streamflow time series at each target
#' location.
#'
#' @param index.network A numeric vector of the streamgage IDs to be used as
#' potential index sites.
#' @param index.obs A zoo object of the observed streamflows for each
#' index gage.
#' @param index.empFDC A matrix of the empirical quantiles for each site in
#' the \code{index.network$index}.
#' @param zero.val (optional) A value to replace zeros.  The deafult is
#' \code{NA}.
#' @param target.obs A zoo object of the observed streamflows for each
#' target gage.
#' @param target.empFDC A matrix of the empirical quantiles for each site in
#' the \code{index.network$target}.
#' @param target.estFDC A matrix of the estimated quantiles for each site in
#' the \code{index.network$target}.
#' @param pvals A vector of the p-values corresponding to the quantiles in
#' \code{index.empFDC}, \code{target.empFDC} and \code{target.estFDC}.
#'
#' @details
#' Lorem ipsum...
#'
#' @return As output, the function returns a list of the same length as
#' \code{index.network$target}.  Each element represents a unique target
#' location.  Each element contains a data.frame consisting of
#' \item{date}{The date of the row.}
#' \item{est}{The estimated streamflow for each day.}
#' \item{indexflow}{The index streamflow for each day.}
#' \item{index}{The index ID used on each day.}
#' \item{index.extrap}{A logical value indicating if the p-value of the index
#' station had to be extrapolated beyond the limits of the index flow duration
#' curve.}
#' \item{target.extrap}{A logical value indicating if the streamflow at the
#' target location had to be extrapolated beyond the limits of the target
#' flow duration curve.}
#' \item{metric}{The similarlity metric between the index and target.  This
#' depends on the metric using in the \code{index.network}; currently
#' functional for nearest-neighbor and map-correlation.}
#' \item{obs}{The observed streamflow at the target location on each day.
#' Used for validation.}
#' \item{obs.pval}{The observed p-value at the target location on each day.
#' Used for validation.}
#'@export
#'@import zoo
estQPPQ <- function(index.network,index.obs,index.empFDC,zero.val=NA,
  target.obs,target.empFDC,target.estFDC,pvals) {
  # Based on code orginally developed by William Farmer, 19 May 2014
  # Revised by William Farmer, 05 June 2015

  # @importFrom zoo as.Date index

  # Setup matrices
  result <- list()
  row.names(index.empFDC) <- NULL
  row.names(target.empFDC) <- NULL
  row.names(target.estFDC) <- NULL
  index.empFDC[index.empFDC==0] <- zero.val
  target.empFDC[target.empFDC==0] <- zero.val
  target.estFDC[target.estFDC==0] <- zero.val
  if (is.vector(target.empFDC)) {
    target.empFDC <- as.matrix(target.empFDC)
  }
  if (is.vector(target.estFDC)) {
    target.estFDC <- as.matrix(target.estFDC)
  }

  for (i in 1:length(index.network$target)) {
    #cat(paste(i,"-",Sys.time()))
    flag <- FALSE
    # Get Ungaged Info
    Site.ID <- index.network$target[i]
    Site.empFDC <- target.empFDC[,i]
    if (sum(is.na(Site.empFDC))==length(Site.empFDC)) {
      Site.empFDC[is.na(Site.empFDC)] <-
        seq(9999,10000,length.out=length(Site.empFDC))
    }
    ndx <- !is.na(Site.empFDC)
    Site.empFDC <- Site.empFDC[ndx]
    Site.empFDC.p <- pvals[ndx]
    Site.estFDC <- target.estFDC[,i]
    if (sum(is.na(Site.estFDC))==length(Site.estFDC)) {
      flag<-TRUE
      Site.empFDC[is.na(Site.empFDC)] <-
        seq(9999,10000,length.out=length(Site.empFDC))
    }
    Site.estFDC <- Site.estFDC[ndx]
    Site.estFDC.p <- pvals[ndx]

    # Potential indices
    indices <- index.network$Rank[,i]

    # Setup empty matrices
    Dates2Est <- index(target.obs[[i]])
    N <- length(Dates2Est)
    Data <- matrix(NA,nrow=N,ncol=10)
    if (flag) {
      Data <- data.frame(Data)
      names(Data) <- c('date','est','indexflow','index.pval','index',
        'index.extrap','target.extrap','metric','obs','obs.pval')
      result[[i]] <- Data
      names(result)[[i]] <- Site.ID
      next
    }
    # {date,est,indexflow,index.pval,index,extrap,metric,obs,obs.pval}
    Data[,c(6,7)] <- FALSE
    Data[,9] <- target.obs[[i]]
    Data[Data[,9]==0,9] <- zero.val
    Data[,10] <- pnorm(approx(log10(Site.empFDC),qnorm(Site.empFDC.p),
      log10(Data[,9]))$y)
    ndx <- which(is.na(Data[,10])&!is.na(Data[,9])&Data[,9]>max(Site.empFDC))
    if (length(ndx)>0) { # upper outliers
      sorted <- sort.int(Site.empFDC,decreasing=T,index.return=T)$ix[c(1,2)]
      intrpQ <- log10(Site.empFDC[sorted])
      intrpP <- qnorm(Site.empFDC.p[sorted])
      m <- (intrpP[1]-intrpP[2])/(intrpQ[1]-intrpQ[2])
      b <- intrpP[1] - m*intrpQ[1]
      Data[ndx,10] <- pnorm(m*log10(Data[ndx,9]) + b)
    }
    ndx <- which(is.na(Data[,10])&!is.na(Data[,9])&Data[,9]<min(Site.empFDC))
    if (length(ndx)>0) { # lower outliers
      sorted <- sort.int(Site.empFDC,decreasing=F,index.return=T)$ix[c(1,2)]
      intrpQ <- log10(Site.empFDC[sorted])
      intrpP <- qnorm(Site.empFDC.p[sorted])
      m <- (intrpP[1]-intrpP[2])/(intrpQ[1]-intrpQ[2])
      b <- intrpP[1] - m*intrpQ[1]
      Data[ndx,10] <- pnorm(m*log10(Data[ndx,9]) + b)
    }
    j <- 0
    temp <- rep(NA,7)
    while (sum(is.na(Dates2Est))<N&j<length(indices)) {
      j <- j + 1
      #cat(paste(j,":",N-sum(is.na(Dates2Est)),"\n"))
      Index <- indices[j]
      ndx <- match(Index,index.network$index)
      Index.empFDC <- index.empFDC[,ndx]
      Index.empFDC[Index.empFDC==0] <- zero.val
      if (sum(is.na(Index.empFDC))==length(Index.empFDC)) {next}
      ndx2 <- !is.na(Index.empFDC)
      Index.empFDC <- Index.empFDC[ndx2]
      Index.empFDC.p <- pvals[ndx2]
      if (is.element("Dist",names(index.network))) {
        imet <- index.network$Dist[ndx,i]
      } else {
        imet <- index.network$Cor.Est[ndx,i]
      }
      # Select index data
      Index.data <- index.obs[[ndx]]
      Index.data[Index.data==0] <- zero.val
      Index.data <- Index.data[!is.na(Index.data)]
      # Find Matching Dates
      NDX <- which(is.element(index(Index.data),Dates2Est))
      if (length(NDX)==0) {next}
      IndexDates <- index(Index.data)[NDX]
      NDX2 <- which(is.element(Dates2Est,IndexDates))
      TargetDates <- Dates2Est[NDX2]
      IndexFlow <- Index.data[NDX]

      # Convert IndexFlows into p-values
      IndexPValue <- pnorm(approx(log10(Index.empFDC),qnorm(Index.empFDC.p),
        log10(IndexFlow))$y)
      ndx.index.extrap <- which(is.na(IndexPValue))
      ndx <- which(is.na(IndexPValue)&IndexFlow>max(Index.empFDC))
      if (length(ndx)>0) { # upper outliers
        sorted <- sort.int(Index.empFDC,decreasing=T,index.return=T)$ix[c(1,2)]
        intrpQ <- log10(Index.empFDC[sorted])
        q <- 2
        while (length(unique(intrpQ))==1) {
          q <- q + 1
          sorted <- sort.int(Index.empFDC,decreasing=T,
            index.return=T)$ix[c(1,q)]
          intrpQ <- log10(Index.empFDC[sorted])
        }
        intrpP <- qnorm(Index.empFDC.p[sorted])
        m <- (intrpP[1]-intrpP[2])/(intrpQ[1]-intrpQ[2])
        b <- intrpP[1] - m*intrpQ[1]
        IndexPValue[ndx] <- pnorm(m*log10(IndexFlow[ndx]) + b)
      }
      ndx <- which(is.na(IndexPValue)&IndexFlow<min(Index.empFDC))
      if (length(ndx)>0) { # lower outliers
        sorted <- sort.int(Index.empFDC,decreasing=F,index.return=T)$ix[c(1,2)]
        intrpQ <- log10(Index.empFDC[sorted])
        q <- 2
        while (length(unique(intrpQ))==1) {
          q <- q + 1
          sorted <- sort.int(Index.empFDC,decreasing=F,
            index.return=T)$ix[c(1,q)]
          intrpQ <- log10(Index.empFDC[sorted])
        }
        intrpP <- qnorm(Index.empFDC.p[sorted])
        m <- (intrpP[1]-intrpP[2])/(intrpQ[1]-intrpQ[2])
        b <- intrpP[1] - m*intrpQ[1]
        IndexPValue[ndx] <- pnorm(m*log10(IndexFlow[ndx]) + b)
      }

      # Convert p-values at ungaged site into streamflows
      TargetFlow <- 10^(approx(qnorm(Site.estFDC.p),log10(Site.estFDC),
        qnorm(IndexPValue))$y)
      ndx.target.extrap <- which(is.na(TargetFlow))
      ndx <- which(is.na(TargetFlow)&IndexPValue>max(Site.estFDC.p))
      if (length(ndx)>0) { # upper outliers
        sorted <- sort.int(Site.estFDC.p,decreasing=T,index.return=T)$ix[c(1,2)]
        intrpQ <- log10(Site.estFDC[sorted])
        intrpP <- qnorm(Site.estFDC.p[sorted])
        m <- (intrpQ[1]-intrpQ[2])/(intrpP[1]-intrpP[2])
        b <- intrpQ[1] - m*intrpP[1]
        TargetFlow[ndx] <- 10^(m*qnorm(IndexPValue[ndx]) + b)
      }
      ndx <- which(is.na(TargetFlow)&IndexPValue<min(Site.estFDC.p))
      if (length(ndx)>0) { # upper outliers
        sorted <- sort.int(Site.estFDC.p,decreasing=F,index.return=T)$ix[c(1,2)]
        intrpQ <- log10(Site.estFDC[sorted])
        intrpP <- qnorm(Site.estFDC.p[sorted])
        m <- (intrpQ[1]-intrpQ[2])/(intrpP[1]-intrpP[2])
        b <- intrpQ[1] - m*intrpP[1]
        TargetFlow[ndx] <- 10^(m*qnorm(IndexPValue[ndx]) + b)
      }

      #if (sum(!is.finite(TargetFlow))>0) {break}

      # Save result
      Data[NDX2,1] <- as.Date(IndexDates) # dates (to double-check)
      Data[NDX2,2] <- TargetFlow # Q_QPPQ
      Data[NDX2,3] <- IndexFlow # indexflow
      Data[NDX2,4] <- IndexPValue # index p-value
      Data[NDX2,5] <- Index # index ID
      Data[NDX2[ndx.index.extrap],6] <- TRUE
      Data[NDX2[ndx.target.extrap],7] <- TRUE
      Data[NDX2,8] <- imet # selection metric

      ndx <- which(!is.finite(TargetFlow))
      NDX2.fin <- NDX2
      if (length(ndx)>0) {
        NDX2.fin <- NDX2[-ndx]
      }
      Dates2Est[NDX2.fin] <- NA
    }
    Data[,9] <- target.obs[[i]]
    Data <- data.frame(Data)
    names(Data) <- c('date','est','indexflow','index.pval','index',
      'index.extrap','target.extrap','metric','obs','obs.pval')
    Data$date <- as.character(Data$date)
    result[[i]] <- Data
    names(result)[[i]] <- Site.ID
  }
  return(result)
}
