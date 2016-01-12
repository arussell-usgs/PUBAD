#' Determine Best Subsets Models
#'
#' @description
#' The function \code{compute.leaps.for} explores the wide
#' set of candidate best-subsets models for each quantile.
#' The function \code{compute.leaps.foreach} accomplishes the same thing using
#' a parallel \code{for} loop.
#'
#' @param unfilled_FDCs A matrix of the raw FDC quantiles for each site.
#' This is derived from the output of \link{\code{calcEmpFDCs}}.
#' @param filled_FDCs A matrix of the censor-filled FDC quantiles
#' for each site.  This is derived from the output of \link{\code{calcEmpFDCs}}.
#' @param comp.wys A vector indicating how many complete water years
#' were present for each site.  This is derived from the output of
#' \link{\code{calcEmpFDCs}}.
#' @param expl A data frame of the potential explanatory variables,
#' derived from the output of \link{\code{getBasinChar}}.
#' @param WY.lim (optional)  The minimum number of water years required to be
#'  included in the formation of regional regressions.  Only reference-quality
#'  sites with at least this many complete water years will be used.
#'  The default is \code{10}.
#' @param zero.val (optional)  The value to which zeroes or negative quantiles
#' will be set to.  The default is \code{0.001}.
#' @param cens_level (optional)  The value establishing left censoring.
#' The deafult is \code{0.005}.
#' @param max.VIF (optional)  The maximum variance inflation factor permitted
#' in a candidate model.  The default is \code{10}.
#' @param nvmax (optional)  The maximum number of variables that will be
#' considered in a regression.  The default is \code{6}.
#' @param nb (optional) The maximum number of models that will be considered
#' for each quantile.  The deafult is \code{3}.
#' @param plunge.ahead (optional) A logical vector indicating if zero values
#' should or should not be handled with the \code{zero.val}.
#' The default is \code{TRUE}.
#'
#' @details
#' This function creates and ranks regression models. It runs regsubsets on
#' the filled FDCs, and then uses those models in \link{\code{censReg}} on the
#' unfilled FDCs. All of the FDCs are transformed using the common logarithm.
#' Any model with a variable that has a VIF higher than the user defined
#' limit \code{max.VIF} will be dropped. The models are automatically ranked
#' by regsubsets' weighted least-sqaures method, but the reported R^2 and BIC
#' are calculated from \link{\code{censReg}}.
#'
#' @return A list of lists, one for each quantile, with each of those
#' containing one array of models per variable level. The models are
#' automatically ranked by regsubsets' weighted least-sqaures method,
#' but the reported R^2 and BIC are calculated from \link{\code{censReg}}.
#'@export
compute.leaps.for=function(unfilled_FDCs,filled_FDCs,comp.wys,
  expl,WY.lim=10,zero.val=0.001, cens_level=.005,
  max.VIF=10, nvmax=6, nb=3, plunge.ahead=T){

  # Function orginially designed by Thomas M. Over and Mike Olsen, 03 July 2015.
  # Modified by William Farmer, 06 July 2015.

  dep = unfilled_FDCs
  dep = t(dep)
  nquants = ncol(dep)

  filled = filled_FDCs
  filled = t(filled)

  #Read in basin characteristics
  expl=expl[as.numeric(rownames(expl)) %in%  as.numeric(rownames(dep)),]
  class=as.character(expl$CLASS)
  #keep only reference gages with >= 10 complete water years
  screen.ndx <- comp.wys>=WY.lim & class=='Ref'
  dep=dep[screen.ndx,]

  #if NA, NaN, or negative, set to .001 (log10(.001)=-3) or
  # return error that takes you out of the function
  # WHF: This could be a problem.  Shouldn't the censor be based on user input?
  if(plunge.ahead==T) {dep[dep <= 0 | is.na(dep) | is.nan(dep)]=zero.val}
  if(plunge.ahead==F & (any(dep <= 0 | is.na(dep) | is.nan(dep)))) {
    stop(paste0("Negative Values, NA's, or NaN's present/n",
      "Please use clean data or set plunge.ahead to ",
      "TRUE to censor these values."))}

  dep=log10(dep)
  cens_level=log10(cens_level)

  filled=log10(filled[screen.ndx,])
  expl=expl[screen.ndx,-1]
  #remove class variable
  expl=as.matrix(expl)
  expl=apply(expl,2,as.numeric)
  #normalize huge variables
  comp.wys=comp.wys[screen.ndx]

  #   cl=makeCluster(detectCores()-1)
  #   registerDoParallel(cl)

  leaps_list = list()
  #loop through quantiles
  for (i in 1:nquants) {
    #create regsubsets object using filled FDCs
    x=regsubsets(x=expl, y = filled[,i], really.big = T,
      weights = comp.wys/mean(comp.wys), nbest=nb, nvmax = nvmax)
    #TRUE/FALSE grid showing BC's and their models
    w_all=summary(x)$which
    w=list()
    length(w)=min(ncol(expl),nvmax)
    #separate model by variable subsets in a list
    var_list=list()
    for (k in 1:length(w)) {
      #split w_all into a list containing models for each variable level
      w[[k]] = w_all[rownames(w_all)==as.character(k),]
      out = array(dim = c((3+3*(ncol(expl)+1)),1))
      out[(1:((length(out)-3)/3))*3 + 1] = c("(Intercept)", colnames(expl))
      out[(1:((length(out)-3)/3))*3 + 2] = "SE"
      out[(1:((length(out)-3)/3))*3 + 3] = "VIF"
      out[1:3] = c("Rank(filled-WLS)", "R^2(censored-MLE)",
        "BIC(censored-MLE)")
      if (is.vector(w[[k]])) {
        model=censReg(as.lcens(dep[,i], cens_level)~
            expl[,colnames(expl)%in%names(which(w[[k]]==TRUE))],
          weights = comp.wys/mean(comp.wys))
        r_2=summary(model)$R2

        if (length(which(smwrStats::vif(model)>max.VIF)) == 0){
          out = cbind(out, rep('', length(out[,1])))
          out[which(out[,1]%in%names(which(w[[k]] == TRUE))==TRUE),
            length(out[1,])] = coef(model) #coeffs
          out[(which(out[,1]%in%names(which(w[[k]] == TRUE))==TRUE)+1),
            length(out[1,])] = sqrt(diag(vcov(model))) #SE
          out[(which(out[,1]%in%names(which(w[[k]] == TRUE))==TRUE)+2),
            length(out[1,])] = c(" ", smwrStats::vif(model)) #VIF
          out[2:3,length(out[1,])] = c(r_2, AIC(model,
            k = log(length(dep[,1])))) #R^2, BIC
        }
      } else {
        for (j in 1:length(w[[k]][,1])) {
          #Use regsubsets models for censReg with unfilled FDCs
          model=censReg(as.lcens(dep[,i], cens_level)~
              expl[,colnames(expl)%in%names(which(w[[k]][j,]==TRUE))],
            weights = comp.wys/mean(comp.wys))
          r_2=summary(model)$R2

          if (length(which(smwrStats::vif(model)>max.VIF)) == 0){
            out = cbind(out, rep('', length(out[,1])))
            out[which(out[,1]%in%names(which(w[[k]][j,] == TRUE))==TRUE),
              length(out[1,])] = coef(model) #coeffs
            out[(which(out[,1]%in%names(which(w[[k]][j,] == TRUE))==TRUE)+1),
              length(out[1,])] = sqrt(diag(vcov(model))) #SE
            out[(which(out[,1]%in%names(which(w[[k]][j,] == TRUE))==TRUE)+2),
              length(out[1,])] = c(" ", smwrStats::vif(model)) #VIF
            out[2:3,length(out[1,])] = c(r_2, AIC(model,
              k = log(length(dep[,1])))) #R^2, BIC
          }
        }
      }
      out[1,][-1] = c(1:(ncol(out)-1))
      var_list[[k]] <- out
    }
    leaps_list[[i]] <- var_list
  }
  return(leaps_list)
}
