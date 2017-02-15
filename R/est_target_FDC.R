#' Estimate a flow duration curve for an ungaged site.
#'
#' @description
#' A sub-function for construction of FDC for an ungaged site.
#'
#' @param best.mods A list of the best models for each flow regime.
#' @param regSelect Statistic used to select best model. The options
#' are \code{adjR2} or \code{AIC}.
#' @param target.empFDC.filled
#' This is derived from the output of \code{\link{calcEmpFDCs}}.
#' @param target.regvar
#' This is derived from the output of \code{\link{getBasinChar}}.
#'
#' @details
#' Lorem ipsum...
#'
#' @return Unclear.
#'@export
#'
#'
est_target_FDC <- function(best.mods, regSelect, target.empFDC.filled, target.regvar)

{
  # Estimate ungauged (target) FDC ####
  # Adapted by Tom Over, 7/2016, from code in runGageList.R written by Will Farmer
  # Modified 1/2017 by AMR to keep mod.mat a matrix after first column deleted.
  #  This was causing errors when custom flow regime was comprised of a single flow quantile.

  target.estFDC <- as.matrix(target.empFDC.filled)
  target.estFDC[] <- NA
  js <- 0
  regimes <- best.mods$regimes
  for (r in 1:length(regimes)) {
    j0 <- max(js+1)
    j1 <- j0 + length(regimes[[r]]) - 1
    js <- j0:j1
    if (regSelect=="adjR2") {# Get best adjusted R2.
      m <- 1
    } else if (regSelect=="AIC") {# Get best AIC.
      m <- 4
    }
    mod.mat <- best.mods$best[[r]][[m]]
    if (sum(is.na(as.double(mod.mat[nrow(mod.mat)-1,-1])))>0
        & regSelect=="adjR2") {
      # If no valid model, switch to AIC.
      m <- 4
    } else if (sum(is.na(as.double(mod.mat[nrow(mod.mat),-1])))>0
               & regSelect=="AIC") {
      # If no valid model, switch to adjusted R2
      m <- 1
    }
    mod.mat <- best.mods$best[[r]][[m]]
    ndx <- !is.element(mod.mat[,1],c("SD","VIF","R^2","AIC"))
    mod.dat <- data.frame(t(apply(as.matrix(mod.mat[ndx,-1]),2,as.double)))
    names(mod.dat) <- mod.mat[ndx,1]
    ndx <- which(is.element(names(target.regvar),names(mod.dat)))
    r.target.regvar <- target.regvar[,ndx]
    target.estFDC[js,] <- mod.dat[,1] +
      as.matrix(mod.dat[,-1])%*%t(as.matrix(r.target.regvar))
  }
  target.estFDC <- as.numeric(10^target.estFDC)
  dim(target.estFDC) <- dim(target.empFDC.filled)
  names(target.estFDC) <- names(target.empFDC.filled)

  return(target.estFDC)
}
