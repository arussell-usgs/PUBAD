#' Determine Optimal Regimes Along FDC
#'
#' @description
#' The function \code{regimes} explores the length of the FDC and breaks it
#' into several regimes, if appropriate.
#'
#' @param unfilled_FDCs A matrix of the raw FDC quantiles for each site.
#' This is derived from the output of \code{\link{calcEmpFDCs}}.
#' @param zero_val (optional)  The value to which zeroes or negative quantiles
#' will be set to.  The default is \code{0.001}.
#'
#' @details
#' Lorem ipsum...
#'
#' @return A matrix containing:
#' \item{rmse}{Root-mean-squared error.}
#' \item{rmsne}{Root-mean-square-normalized error.}
#' \item{cvrmse}{The coefficient of variation of the root-mean-squared error.}
#' \item{nrmse}{Normalized root-mean-squared error.}
#'@export
#'@import pastecs
regimes = function(unfilled_FDCs,zero_val=0.001){

  # Function orginially designed by Thomas M. Over and
  #   Mike Olsen, 06 October 2015.
  # Modified by William Farmer, 09 October 2015.
  #     Re-written to prevent writing key files; returns output instead.


  # @importFrom pastecs turnpoints

  gages = colnames(unfilled_FDCs)
  probs <- as.double(row.names(empFDCs.all))
  BCs_path <- file.path("Data","Raw",
    "gagesII_PUBairdropchar_full_UTM.csv")
  BCs = read.csv(BCs_path)
  BCs = BCs[as.numeric(BCs[,1]) %in% as.numeric(gages),]

  red.FDCs = unfilled_FDCs; red.BCs = BCs

  unit.red.FDCs = red.FDCs
  unit.red.FDCs[unit.red.FDCs <= 0] = zero_val
  unit.conv = 24*3600*1.60934^2*12*25.4/5280^2 #converts ft^3/s/km^2 to mm/day
  for (z in 1:length(gages)) {
    unit.red.FDCs[,z] = unit.conv*red.FDCs[,z]/red.BCs$DRAIN_SQKM[z]
  }
  #fill in zeroes with designated positive value so y-axis can be log-transformed
  pos.subs = NULL
  pos.subs = which(unit.red.FDCs > .005*24*3600*1.60934^2*12*25.4/5280^2)
  min.pos.val = min(unit.red.FDCs[pos.subs])
  full.unit.FDCs = unit.red.FDCs
  zero.subs = which(unit.red.FDCs<=.005*24*3600*1.60934^2*12*25.4/5280^2)
  full.unit.FDCs[zero.subs] = 10^(floor(log10(min.pos.val))-1)

  poly_model = lm(apply(log10(full.unit.FDCs), 1, sd) ~ poly(qnorm(probs),
    3, raw=TRUE))

  x <- qnorm(probs)
  coefs1 = summary(poly_model)$coefficients[2,1]
  coefs2 = summary(poly_model)$coefficients[3,1]
  coefs3 = summary(poly_model)$coefficients[4,1]
  d2 = deriv3(~ coefs1*x + coefs2*x^2 + coefs3*x^3, "x")

  sd_fit = fitted(poly_model)

  t_points = turnpoints(sd_fit)

  if(any(attributes(eval(d2))$hessian <= 0 &
      any(attributes(eval(d2))$hessian >= 0))) {
    if(attributes(eval(d2))$gradient[which(abs(attributes(eval(d2))$hessian)==
        min(abs(attributes(eval(d2))$hessian)))] < 0){
      low = approx(attributes(eval(d2))$hessian, probs, xout = 0)$y
      if(any(t_points$pits)) {
        high = approx(attributes(eval(d2))$gradient[probs > low],
          probs[probs > low], xout = 0)$y
      }
    }
  } else if(any(t_points$pits)) {
    high = approx(attributes(eval(d2))$gradient, probs, xout = 0)$y
  }
  if(exists('low') & exists('high')) {
    regimes = list(lowflow = probs[probs < low],
      medflow = probs[probs >= low & probs < high],
      highflow = probs[probs >= high])
  }
  if(exists('low') & !exists('high')) {
    regimes = list(lowflow = probs[probs < low],
      highflow = probs[probs >= low])
  }
  if(!exists('low') & exists('high')) {
    regimes = list(lowflow = probs[probs < high],
      highflow = probs[probs >= high])
  }
  if(!exists('low') & !exists('high')) {
    regimes = list(lowflow = probs)
  }

  hi = ifelse(exists("high"), high, "-")
  lo = ifelse(exists("low"), low, "-")
  summ = rbind(summary(poly_model)$coefficients, rep("", 4),
    "R^2" = c(signif(summary(poly_model)$r.squared, 4), rep("", 3)),
    "Sqrt Mean Variance" = c(signif(sqrt(mean(apply(log10(full.unit.FDCs),
      1, var))), 4), rep("", 3)),
    "Low Break" = lo, "High Break" = hi)

  result <- list(regimes=regimes,summary=summ)
  return(result)
}
