#' Apply two version of Drainage Area Ratio
#'
#' @description
#' The function \code{applyDAR} conducts a cross-validated evaluation of
#' nearest-neighbor and map-correlation drainage area ratio.
#'
#' @param listofgages A dataframe of sites to be used for analysis.  Sites are
#' character strings with at least 8 numeric characters.
#' @param startDate Specifies a particular date from which to begin retrieving
#' the streamflow data.  Format should be "YYYY-MM-DD".  If blank, all available
#' data will be retrieved.
#' @param endDate Specifies a particular date up to which to retrieve
#' the streamflow data.  Format should be "YYYY-MM-DD".  If blank, all available
#' data will be retrieved.
#' @param limWY (optional) The minimum number of water years required for a
#' site to be included. The default is \code{1}.
#' @param limDAR (optional) A numeric vector with two elements, indicating the
#' minimum and maximum drainage area ratio to be allowed for the applicaiton
#' of drainage area ratio. The default is \code{c(-Inf,Inf)}.
#' @param cv (optional) A list of two elements defining the parameters of
#' cross validation.
#' \itemize{
#' \item{fold}{An integer indicating how many folds should be used for
#' evaluation.  The default, \code{NA} or \code{1}, requests a leave-one-out
#' analysis.}
#' \item{dense}{A logical value indicating if the smaller or larger half of the
#' cross-validation fold should be used as the index.  The defulat, \code{TRUE},
#'  calls for the dense scenario using the larger half of the fold as
#'  the index network.}
#' }
#' @param flowStat (optional) A logical indicating if the streamflow
#' statistics should be evaluated.  This will increase run time.  The default
#' is \code{FALSE}.
#'
#' @details
#' Lorem ipsum...
#'
#' @return A list of two elements:
#' \item{cvList}{Corss-validated estiamted and observed streamflows.}
#' \item{anList}{basic analysis of estiamted streamflows.}
#'@export
applyDAR <- function(listofgages,startDate,endDate,
  limWY=1,limDAR=c(-Inf,Inf),cv=list(fold=NA,dense=TRUE),flowStat=FALSE) {


  # Load the streamflow data
  Streamflow <- readCompileFlow(listofgages,units="cfs",
    startDate=startDate,endDate=endDate,checkmiss=FALSE,keepNAdays=TRUE)

  # Load basin characteristics
  expVars <- getBasinChar(listofgages)

  # Check for complete water years
  empFDC.calc <- calcEmpFDCs(flow_data=Streamflow)
  FDC.sites <- as.double(colnames(empFDC.calc$empFDC))
  listofgages.new <- listofgages
  ndx <- which(empFDC.calc$PORstats.compWyr[,2]>=limWY)
  ndx2 <- which(is.element(as.double(as.character(listofgages$gages)),
    FDC.sites[ndx]))
  listofgages.new <- data.frame(listofgages[ndx2,])
  names(listofgages.new) <- names(listofgages)
  listofgages.old <- listofgages
  listofgages <- listofgages.new # revised listofgages

  # set up CV scenario (LOO default)
  if (is.na(cv$fold)|cv$fold==1|cv$fold>length(listofgages$gages)) {
    cv$fold=length(listofgages$gages)
  }
  if (cv$fold == length(listofgages$gages)) {
    scenNums <- 1:length(listofgages$gages)
    cv$dense <- TRUE
  } else {
    sampSize <- round(length(listofgages$gages)/cv$fold)
    scenNums <- rep(NA,length(listofgages$gages))
    for (i in 1:(cv$fold-1)) {
      ndx <- sample(which(is.na(scenNums)),sampSize)
      scenNums[ndx] <- i
    }
    scenNums[is.na(scenNums)] <- cv$fold
  }

  # Looping through CV scenarios
  cat('\n')
  cvList <- list()
  for (i in 1:cv$fold) {
    stime<-Sys.time()
    # Set up index and target according to 'dense' specification.
    if (cv$dense) {
      index.gages <- as.double(as.character(listofgages$gages[scenNums!=i]))
      target.gages <- as.double(as.character(listofgages$gages[scenNums==i]))
    } else {
      index.gages <- as.double(as.character(listofgages$gages[scenNums==i]))
      target.gages <- as.double(as.character(listofgages$gages[scenNums!=i]))
    }
    # Get target information
    ndx <- match(target.gages,expVars$DescBCs$STAID)
    target.baschar <- expVars$DescBCs[ndx,]
    ndx <- match(target.gages,as.double(names(Streamflow)))
    target.obs <- Streamflow[ndx]

    # Get index information
    ndx <- match(index.gages,expVars$DescBCs$STAID)
    index.baschar <- expVars$DescBCs[ndx,]
    ndx <- match(index.gages,as.double(names(Streamflow)))
    index.obs <- Streamflow[ndx]

    # Determine Nearest Neighbor ####
    NN.ranking <- indexNN(index.gages,index.baschar,
      target.gages,target.baschar)

    # Determine Map Correlation ####
    MC.ranking <- indexMC(index.gages,index.obs,index.baschar,
      target.gages,target.obs,target.baschar,
      method="pearson")

    # NN-DAR
    NNDAR <- estDAR(index.network=NN.ranking,index.baschar,
      target.baschar,index.obs,target.obs,DAR.low=min(limDAR),DAR.high=max(limDAR))

    # MC-DAR
    MCDAR <- estDAR(index.network=MC.ranking,index.baschar,
      target.baschar,index.obs,target.obs,DAR.low=min(limDAR),DAR.high=max(limDAR))

    # Save results
    cvList[[i]] <- list(NNDAR=NNDAR,MCDAR=MCDAR)
    ltime <- Sys.time()
    etime <- ltime - stime
    cat(paste0("CV.fold ",i,": Started ",stime,
      "; Completed ",ltime,"(",round(etime,digits=2),")\n"))
  }

  # Analyze results
  sets <- c("NNDAR","MCDAR")
  AnList <- list()
  for (i in 1:length(sets)) {
    stime<-Sys.time()
    data <- eval(parse(text=
        paste0("c(",paste0("cvList[[",c(1:length(cvList)),
          "]]$",sets[i],collapse=","),")")))
    AnList[[i]] <- analyzeResult(data,FlowStat=flowStat)
    ltime <- Sys.time()
    etime <- ltime - stime
    cat(paste0("Completed analysis of ",sets[i],": Started ",stime,
      "; Completed ",ltime,"(",round(etime,digits=2),")\n"))
  }
  names(AnList) <- sets

  # return results
  result <- list(cvList=cvList,anList=AnList)
  return(result)
}
