#' Compile models across quantiles and regimes
#'
#' @description
#' The function \code{compile.vars} takes output from
#' \code{\link{compute.leaps.for}} and compiles the candidate models across
#' specified quantiles and regimes.
#'
#' @param leaps_list Output from \code{\link{compute.leaps.for}}
#' @param nvmax The maximum number of variables that will be
#' considered in a regression.  Prior to 9/2016, the default was \code{6}.
#' @param n The maximum number of models subsetted.
#' Prior to 9/2016, the default was \code{20}.
#' @param regimes (optional) A list of two or more elements (\code{names(regimes) =
#' c('lowflow','medflow','highflow')}), where each element is a character
#' vector indicating the frequencies of each quantile. The default is \code{
#' list(lowflow = c("0.0002","0.0005","0.001","0.002","0.005","0.01","0.02",
#' "0.05","0.1"),medflow = c("0.2", "0.25", "0.3", "0.4", "0.5","0.6", "0.7",
#'  "0.75", "0.8", "0.9"),highflow = c("0.95", "0.98", "0.99", "0.995",
#'  "0.998", "0.999", "0.9995", "0.9998"))}.
#'
#' @details
#' Lorem ipsum...
#'
#' @return A list ranking all models across the number of variables
#' for each regime. Higher ranks are better.
#'@export
compile.vars <- function(leaps_list, nvmax, n,
  regimes = list(
    lowflow = c("0.0002","0.0005","0.001","0.002",
      "0.005","0.01","0.02","0.05","0.1"),
    medflow = c("0.2", "0.25", "0.3", "0.4", "0.5",
      "0.6", "0.7", "0.75", "0.8", "0.9"),
    highflow = c("0.95", "0.98", "0.99", "0.995",
      "0.998", "0.999", "0.9995", "0.9998")),
  unfilled_FDCs=NA)
{
  # Function originally designed by Tom Over and Mike Olson, 03 July 2015.
  # Modified by William Farmer, 06 July 2015.
  # Revised by TMO, June 2016
  # Implemented by Amy Russell, 9/2016


  # Following if statement deleted by AMR, 9/2016
  #  due to new requirement that BCs be passed to the regimes function.
  # When not using default regimes, suggest user call regimes function first then
  #  pass its result to this function (compile.vars) rather than calling regimes
  #  from within the code here.


  # # if (is.na(regimes)) { #modified by TMO, 6/2016
  # if (sum(is.na(regimes))>0) {
  #   regimes.ret <- regimes(unfilled_FDCs,BCs,zero_val=zero.val)
  #   regimes <- regimes.ret$regimes
  # }

  # compute range from regimes (added by TMO, 6/2016)
  range <- regimes
  rngcnt <- 0
  for(i in 1:length(regimes)) {
    range[[i]] = (rngcnt+1):(rngcnt+length(regimes[[i]]))
    rngcnt = rngcnt+length(regimes[[i]])
  }

  top_n_list=list()
  var_list=list()
  length(top_n_list) = length(regimes)
  length(var_list)=nvmax


  for(i in 1:length(regimes)){
    for(j in 1:min(nvmax,length(leaps_list[[i]]))){
      if (ncol(leaps_list[[i]][[j]])==1) {next}
      #creates a list of 3: the Predictor variables in each model,
      #each model's rank, and all the explanatory variables

      #Revised by TMO, 6/2016:
      ret = tolist(regimes = regimes, regm.indx=i, vars = j, ranges = range,
        leaps_list = leaps_list)
      y = ret[[1]]
      ranks = ret[[2]]
      expl = ret[[3]]

      out = array(" ", dim = c((length(expl)+ 3),1))
      out[,1] = c(as.character(expl), "Rank", "Num", "Quantile")
      #compares all models against each other and sums the ranks between
      # quantiles for each model, recording how many quantiles
      # the model appears in
      output = compare(out, y, ranks, vars=j, names = regimes[[i]])
      #order models by their summed ranks, higher ranks are better
      ord = c(0,sort.list(as.numeric(output[nrow(output)-2,
        2:length(output[1,])]), decreasing=T)) + 1
      top=output[,c(ord)]
      #subset top n models, while accounting for the flow regimes that might
      # have less models than the top number specified
      if (is.null(top)) {
        top_n=NA
      #} else if(ncol(top) >= n)  {
      #Replacing line above which seems to ignore that top has a column of variable names on left,
      #so only ncol(top)-1 columns contain the actual "data"
      } else if((ncol(top)-1) >= n)  {
        top_n=top[,1:(n+1)]
      } else {
        top_n=top
      }
      var_list[[j]]=top_n
    }
    top_n_list[[i]]=var_list
  }
  return(top_n_list)
}
