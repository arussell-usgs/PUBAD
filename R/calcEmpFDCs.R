# calcEmpFDCs ----
#' Calculate empirical FDCs for all sites
#'
#' @description
#' The function \code{nse} computes the Nash-Sutcliffe efficiency for the rows
#' or columns of a matrix.
#'
#' @param flow_data A zoo object of streamflow data.
#' @param quant.type (optional) An integer, from 1 to 9, defining the
#' formula used to compute quantiles. The default is \code{9}.
#' @param cens (optional) A character string, either \code{'filled'} or
#' \code{'unfilled'}, specifying how censored values will be handled.
#' The default is \code{'unfilled'}.
#' @param cens_level (optional)  The number specifying the censoring level.
#' The default is \code{0.005}.
#' @param begWyr (optional) The lower limit of water years to consider.
#' The default is \code{1900}.
#' @param endWyr (optional) The upper limit of water years to consider.
#' The default is \code{2014}.
#' @param saveName (optional) A character string of a file name, without
#' extension, in which to save results.  An empty string, the default, requests
#' that nothing be saved.
#' @param FDC.probs (optional) THe probabilities of the quantiles to be
#' calculated.  The default is \code{c(0.0002,0.0005,0.001,0.002,0.005,0.01,
#' 0.02,0.05,0.10,0.20,0.25,0.30,0.40,0.50,0.60,0.70,0.75,0.80,0.90,0.95,
#' 0.98,0.99,0.995,0.998,0.999,0.9995,0.9998)}.
#'
#' @details
#' This function computes flow duration curves and summary statistics for a
#' given time series (aka zoo object) of stream flows.
#' It outputs the FDCs directly, but can also write to external files
#' via \code{saveName}.
#'
#' If \code{cens} is set to \code{'filled'}, \code{library(smwrQW)} is required.
#' (Install via instructions found here: https://github.com/USGS-R/smwrQW
#' as the package is not on CRAN.)
#'
#' @return A large list of many internal parameters.  The most salient are:
#' \item{empFDC}{A matrix of the empircal flow duration curve quantiles.}
#' \item{PORstats.compWyr}{The number of complete water
#' years available at each site.}
#'@export
calcEmpFDCs <- function(flow_data, quant.type=9, cens='unfilled',
  cens_level=.005, begWyr=1900, endWyr=2014, saveName='',
  FDC.probs=c(0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.10,0.20,0.25,
    0.30,0.40,0.50,0.60,0.70,0.75,0.80,0.90,0.95,0.98,0.99,0.995,0.998,0.999,
    0.9995,0.9998))
{
  # Function orginially designed by Thomas M. Over and Mike Olsen, 05 June 2015.
  # Modified by William Farmer, 08 June 2015.

  # List inputs for saving later
  inputs <- list(flow_data=flow_data, quant.type=quant.type, cens=cens,
    cens_level=cens_level, begWyr=begWyr, endWyr=endWyr, saveName=saveName,
    FDC.probs=FDC.probs)

  #Compute number of days per water year
  glb.Wyrs = begWyr:endWyr
  glb.Wyr.ndays = numeric(length(glb.Wyrs))
  names(glb.Wyr.ndays) = glb.Wyrs
  for (i in 1:length(glb.Wyrs)){
    glb.Wyr.ndays[i] = as.Date(paste(glb.Wyrs[i],10,01,sep="-")) -
      as.Date(paste(glb.Wyrs[i]-1,10,01,sep="-"))
  }

  ngages = length(flow_data)
  gage.names = names(flow_data)

  #Create various storage arrays
  #Arrays containing statistics for all water years
  nvals.pWyr = matrix(nrow=length(glb.Wyrs), ncol=ngages,
    dimnames=list(glb.Wyrs,gage.names))
  comp.Wyr = nNAs.pWyr = nzeroes.pWyr = nnegs.pWyr = nvals.pWyr
  #Arrays containing statistics for complete water years
  nzeroes.pcompWyr = nnegs.pcompWyr = nvals.pcompWyr = nvals.pWyr
  nvals.compWyr = max_val.compWyr = min_val.compWyr = frac_negs.compWyr =
    frac_zeroes.compWyr = numeric(ngages)
  names(nvals.compWyr) = names(max_val.compWyr) = names(min_val.compWyr) =
    names(frac_negs.compWyr) = names(frac_zeroes.compWyr) = gage.names
  nzeroes.pcompWyr[] = nnegs.pcompWyr[] = nvals.pcompWyr[] = nvals.pWyr[]
  nvals.compWyr[] = max_val.compWyr[] = min_val.compWyr[] =
    frac_negs.compWyr[] = frac_zeroes.compWyr[] = NA

  #Set up FDC arrays
  nFDC.probs = length(FDC.probs)
  FDCs = matrix(nrow=nFDC.probs, ncol=ngages,
    dimnames=list(FDC.probs,gage.names))

  #Loop through stations and do computations
  for (j in 1:ngages) {
    #Extract date and streamflow information
    dates = as.Date(index(flow_data[[j]]))
    ndays = length(dates)
    flows = flow_data[[j]]
    #Set up FDC.data array that will be filled each complete water year
    FDC.data = numeric(sum(glb.Wyr.ndays))
    FDC.data[] = NA
    FDC.cnt = 0
    for (yr in 1:length(glb.Wyrs)) {
      beg.date = as.Date(paste(glb.Wyrs[yr]-1,10,01,sep="-"))
      end.date = as.Date(paste(glb.Wyrs[yr],09,30,sep="-"))
      date.inds = (dates>=beg.date & dates<=end.date)
      wyr.flows = flows[date.inds]
      #Compute number of days with data in series this water year
      #(Perhaps could use zoo function dwi() but it does not work
      #when applied as dwi(series_f[[j]]) for some reason)
      NA.cnt = sum(is.na(flows[date.inds]))
      zero.cnt = sum(wyr.flows==0)
      neg.cnt = sum(wyr.flows<0)
      flow.cnt = sum(date.inds) - NA.cnt
      nvals.pWyr[yr,j] = flow.cnt
      nNAs.pWyr[yr,j] = NA.cnt
      nzeroes.pWyr[yr,j] = zero.cnt
      nnegs.pWyr[yr,j] = neg.cnt
      #Water year is complete if nvals.pWyr matches glb.Wyr.ndays
      if (glb.Wyr.ndays[yr]==nvals.pWyr[yr,j]) {
        comp.Wyr[yr,j] = T
      } else {
        comp.Wyr[yr,j] = F
      }
      #If water year is complete
      if (comp.Wyr[yr,j]) {
        #Compile data for FDC computation if WY is complete
        FDC.data[(FDC.cnt+1):(FDC.cnt+flow.cnt)] = wyr.flows
        FDC.cnt = FDC.cnt+flow.cnt
        #Other statistics
        nvals.pcompWyr[yr,j] = flow.cnt
        nzeroes.pcompWyr[yr,j] = zero.cnt
        nnegs.pcompWyr[yr,j] = neg.cnt
      }
    }
    #If there were any complete WYs
    if (FDC.cnt>0) {
      #FDC computation
      FDC.data = FDC.data[1:FDC.cnt] #all the daily flows from the complete WYs
      #Other statistics over complete water year data
      nvals.compWyr[j] = FDC.cnt
      max_val.compWyr[j] = max(FDC.data)
      min_val.compWyr[j] = min(FDC.data)
      frac_negs.compWyr[j] = sum(FDC.data<0)/FDC.cnt
      frac_zeroes.compWyr[j] = sum(FDC.data==0)/FDC.cnt
      # If all one value, ignore site...
      if (length(unique(FDC.data))==1) {
        nvals.compWyr[j] = 0
        next
      }
      #Fill in censored quantile values
      if(cens=='filled') {
        quants = quantile(as.lcens(FDC.data, cens_level),FDC.probs,
          type=quant.type, method='log MLE', alpha=3/8)
      } else if(cens=='unfilled') {
        quants = quantile(FDC.data,FDC.probs,names=F,type=quant.type)
      }
      #quantile() returns values for all probs even if
      #  extrapolation is required.
      #go back and insert NA values where extrapolation was applied
      if (quant.type==4) {
        NA.subs = which(FDC.probs < 1/FDC.cnt)
        if (length(NA.subs)>0) quants[NA.subs] = NA
      } else if (quant.type==5) {
        NA.subs = c(which(FDC.probs < 0.5/FDC.cnt),
          which(FDC.probs > (FDC.cnt-0.5)/FDC.cnt))
        if (length(NA.subs)>0) quants[NA.subs] = NA
      } else if (quant.type==6) {
        NA.subs = c(which(FDC.probs < 1/(FDC.cnt+1)),
          which(FDC.probs > FDC.cnt/(FDC.cnt+1)))
        if (length(NA.subs)>0) quants[NA.subs] = NA
      } else if (quant.type==7) {
        #No need to do anything: min and max are 0 and 1 respsectively.
      } else if (quant.type==8) {
        NA.subs = c(which(FDC.probs < (2/3)/(FDC.cnt+1/3)),
          which(FDC.probs > (FDC.cnt-1/3)/(FDC.cnt+1/3)))
        if (length(NA.subs)>0) quants[NA.subs] = NA
      } else if (quant.type==9) {
        NA.subs = c(which(FDC.probs < (5/8)/(FDC.cnt+1/4)),
          which(FDC.probs > (FDC.cnt-3/8)/(FDC.cnt+1/4)))
        if (length(NA.subs)>0) quants[NA.subs] = NA
      }
      FDCs[,j] = quants
    } else {
      nvals.compWyr[j] = 0
    }
  }

  #1-d arrays containing statistics for complete water years combined
  nWys.compWyr = round(nvals.compWyr/365.25)
  PORstats.compWyr <- cbind(nvals.compWyr,nWys.compWyr,max_val.compWyr,
    min_val.compWyr,frac_negs.compWyr,frac_zeroes.compWyr)

  if (saveName!='') {# User requests saved output.
    saveFile <- paste0(saveName,'.RData')
    varList <- c("nvals.pWyr","nvals.pWyr","nNAs.pWyr","nzeroes.pWyr",
      "nnegs.pWyr","comp.Wyr","nvals.pcompWyr","nzeroes.pcompWyr",
      "nnegs.pcompWyr","PORstats.compWyr","inputs","FDCs")
    save(file=saveFile,list=varList)
  }

  result <- list(empFDC=FDCs,nvals.pWyr=nvals.pWyr,nvals.pWyr=nvals.pWyr,
    nNAs.pWyr=nNAs.pWyr,nzeroes.pWyr=nzeroes.pWyr,nnegs.pWyr=nnegs.pWyr,
    comp.Wyr=comp.Wyr,nvals.pcompWyr=nvals.pcompWyr,
    nzeroes.pcompWyr=nzeroes.pcompWyr,nnegs.pcompWyr=nnegs.pcompWyr,
    PORstats.compWyr=PORstats.compWyr,inputs=inputs)
  return(result)
}
