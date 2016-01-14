#' Compile models across quantiles and regimes
#'
#' @description
#' The function \code{compile.vars} takes output from
#' \code{\link{compute.leaps.for}} and compiles the candidate models across
#' specified quantiles and regimes.
#'
#' @param leaps_list Output from \code{\link{compute.leaps.for}}
#' @param nvmax (optional)  The maximum number of variables that will be
#' considered in a regression.  The default is \code{6}.
#' @param nb (optional) The maximum number of models subsetted.
#' The deafult is \code{20}.
#' @param range A list of three elements (\code{names(range) = c('low_range',
#' 'med_range','high_range')}), where each element is a vector of indices
#' specify the quantiles belonging to each regime.  The default is \code{
#' list(low_range=c(1:9), med_range=c(10:19), high_range=c(20:27))}.
#' @param regimes (optional) A list of three elements (\code{names(regimes) =
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
#' for each regime.
#'@export
compile.vars <- function(leaps_list,nvmax = 6, n=20,
  range = list(low_range=c(1:9), med_range=c(10:19), high_range=c(20:27)),
  regimes = list(
    lowflow = c("0.0002","0.0005","0.001","0.002",
      "0.005","0.01","0.02","0.05","0.1"),
    medflow = c("0.2", "0.25", "0.3", "0.4", "0.5",
      "0.6", "0.7", "0.75", "0.8", "0.9"),
    highflow = c("0.95", "0.98", "0.99", "0.995",
      "0.998", "0.999", "0.9995", "0.9998")),
  unfilled_FDCs=NA)
{
  # Function orginially designed by Thomas M. Over and Mike Olsen, 03 July 2015.
  # Modified by William Farmer, 06 July 2015.
  if (is.na(regimes)) {
    regimes.ret <- regimes(unfilled_FDCs,zero_val=zero.val)
    regimes <- regimes.ret$regimes
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
      ret = tolist(names = regimes[[i]], vars = j, range = range[[i]],
        leaps_list=leaps_list)
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
      } else if(ncol(top) >= n)  {
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
