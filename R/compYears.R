#' Function to determine complete water years
#'
#' @description
#' The function \code{compYears} looks through the data to determine which
#' years are full of real values.
#'
#' @param dates A vector of \code{class(dates)='Date'}.
#' @param data Either a vector of the same length as \code{dates} or a
#' two-dimensional array with one dimension equal to the length of \code{dates}.
#' @param year A character string from \code{c('water','calendar','climate')}
#' indicating which type of year should be evaluated.  Water years run from
#' October 01 through September 30.  Calendar years run from January 01
#' through December 31.  Climate years run from April 01 through March 31.
#'
#' @details
#' Lorem ipsum...
#'
#' @return The function retruns a list of objects:
#' \item{logiY}{A logical data.frame of the same dimensions as \code{data}.
#' \code{TRUE} indicates that the element is a member of a complete wate year.}
#' \item{numY}{The number of complete years from each site.}
#' \item{sY}{A numeric vector naming the nominal year of
#' each element in \code{dates}.}
#' \item{year}{The same as the input, for validation.}
#'@export
#'@import zoo
compYears <- function(dates,data,year=c('water','calendar','climate')) {
  # Developed by William Farmer, USGS, 15 January 2016

  # Check that dates and data share a dimension
  if (sum(is.na(dates))>0) stop('DATES cannot contain NA')
  if (is.vector(data)) {
    data <- as.matrix(data)
  }
  if (ncol(data)!=length(dates)&nrow(data)!=length(dates)) {
    stop('DATES and DATA  must have one identical dimension')
  }
  data1 <- data
  if (nrow(data)!=length(dates)) {
    data <- t(data)
  }
  i1 <- ncol(data)

  # Pull out years and months
  Y <- as.double(substr(dates,1,4))
  M <- as.double(substr(dates,6,7))
  # determine full range of dates
  if (year=='water') {
    minDate <- as.Date(paste(min(Y)-1,'10','01',sep='-'))
    maxDate <- as.Date(paste(max(Y)+1,'09','30',sep='-'))
  } else if (year=='calendar') {
    minDate <- as.Date(paste(min(Y),'01','01',sep='-'))
    maxDate <- as.Date(paste(max(Y),'12','31',sep='-'))
  } else if (year=='climate') {
    minDate <- as.Date(paste(min(Y)-1,'04','01',sep='-'))
    maxDate <- as.Date(paste(max(Y)+1,'03','31',sep='-'))
  }
  fullDates <- seq(minDate,maxDate,by='days')
  #fullDates <- seq(min(dates),max(dates),by='days')
  fY <- as.double(substr(fullDates,1,4))
  fM <- as.double(substr(fullDates,6,7))

  # Assign the special year values
  if (year=='water') {
    fsY <- ifelse(fM>9,fY+1,fY)
    sY <- ifelse(M>9,Y+1,Y)
  } else if (year=='calendar') {
    fsY <- fY
    sY <- Y
  } else if (year=='climate') {
    fsY <- ifelse(fM>3,fY,fY-1)
    sY <- ifelse(M>3,Y,Y-1)
  }

  # table of days in full years (to handle leap years)
  fYtable <- cbind(as.numeric(names(table(fsY))),table(fsY))
  row.names(fYtable) <- NULL
  ndx <- fYtable[,2]<365 # partial year from input
  fYtable <- fYtable[!ndx,]

  # Loop through all sites
  returnNDX <- matrix(NA,nrow=length(dates),ncol=i1)
  numY <- rep(NA,i1)
  for (s in 1:i1) {
    isY <- sY[!is.na(data[,s])]
    iYtable <- cbind(as.numeric(names(table(isY))),table(isY))
    row.names(iYtable) <- NULL
    ifYtable <- fYtable[which(is.element(fYtable[,1],iYtable[,1])),]
    isYs <- as.numeric(iYtable[iYtable[,2]==ifYtable[,2],1])
    numY[s] <- length(isYs)
    returnNDX[,s] <- is.element(sY,isYs)
  }

  if (nrow(data1)!=length(dates)) {
    returnNDX <- t(returnNDX)
  }

  # format result for output
  output <- list(logiY=returnNDX,numY=numY,sY=sY,year=year)
  return(output)
}
