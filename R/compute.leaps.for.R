#' Determine Best Subsets Models
#'
#' @description
#' The function \code{compute.leaps.for} explores the wide
#' set of candidate best-subsets models for each quantile.
#' The function \code{compute.leaps.foreach} accomplishes the same thing using
#' a parallel \code{for} loop.
#'
#' @param unfilled_FDCs A matrix of the raw FDC quantiles for each site.
#' This is derived from the output of \code{\link{calcEmpFDCs}}.
#' @param filled_FDCs A matrix of the censor-filled FDC quantiles
#' for each site.  This is derived from the output of \code{\link{calcEmpFDCs}}.
#' @param comp.wys A vector indicating how many complete water years
#' were present for each site.  This is derived from the output of
#' \code{\link{calcEmpFDCs}}.
#' @param expl A data frame of the potential explanatory variables,
#' derived from the output of \code{\link{getBasinChar}}.
#' @param nvmax The maximum number of variables that will be
#' considered in a regression.  When specifying this number, the user should
#' take into consideration the number of sites. Prior to 9/2016 the default was \code{6}.
#' @param nb The maximum number of models that will be considered
#' for each quantile.  Prior to 9/2016 the default was \code{3}.
#' @param WY.lim (optional)  The minimum number of water years required to be
#'  included in the formation of regional regressions.  Only reference-quality
#'  sites with at least this many complete water years will be used.
#'  The default is \code{10}.
#' @param zero.val (optional)  The value to which zeroes or negative quantiles
#' will be set to.  The default is \code{0.001}.
#' @param cens_level (optional)  The value establishing left censoring.
#' The default is \code{0.005}.
#' @param max.VIF (optional)  The maximum variance inflation factor permitted
#' in a candidate model.  The default is \code{10}.
#' @param plunge.ahead (optional) A logical vector indicating if zero values
#' should or should not be handled with the \code{zero.val}.
#' The default is \code{TRUE}.
#' @param forced.BC (optional) A single basin characteristic that candidate models
#' are required to include.
#' The default is \code{NULL}.
#'
#' @details
#' This function creates and ranks regression models. It runs regsubsets on
#' the filled FDCs, and then uses those models in \code{\link{censReg}} on the
#' unfilled FDCs. All of the FDCs are transformed using the common logarithm.
#' Any model with a variable that has a VIF higher than the user defined
#' limit \code{max.VIF} will be dropped. The models are automatically ranked
#' by regsubsets' weighted least-sqaures method, but the reported R^2 and BIC
#' are calculated from \code{\link{censReg}}.
#'
#' @return A list of lists, one for each quantile, with each of those
#' containing one array of models per variable level. The models are
#' automatically ranked by regsubsets' weighted least-sqaures method,
#' but the reported R^2 and BIC are calculated from \code{\link{censReg}}.
#'@export
#'@import leaps
#'@import smwrStats
#'@import smwrQW
compute.leaps.for=function(unfilled_FDCs,filled_FDCs,comp.wys,
  expl, nvmax, nb, WY.lim=10, zero.val=0.001, cens_level=0.005,
  max.VIF=10, plunge.ahead=T, forced.BC=NULL){

  # Function orginially designed by Tom Over and Mike Olson, 03 July 2015.
  # based on prior work by Elizabeth Cleveland and Riten Patel.
  # Modified by William Farmer, 06 July 2015.

  # Modified by TMO and Amy Russell to allow forcing in a
  # basin characteristic (e.g., drainage area) using the added parameter "forced.BC".

  # @importFrom leaps regsubsets
  # @importFrom smwrStats vif
  # @importFrom smwrQW censReg summary.censReg

  dep = unfilled_FDCs
  dep = t(dep)
  nquants = ncol(dep)

  filled = filled_FDCs
  filled = t(filled)

  #Read in basin characteristics
  expl=expl[as.numeric(rownames(expl)) %in%  as.numeric(rownames(dep)),]
  #class=as.character(expl$CLASS) #removed hard-coding 9/2016 by AMR
  #keep only reference gages with >= 10 complete water years
  screen.ndx <- comp.wys>=WY.lim #& class=='Ref' #removed hard-coding 9/2016 by AMR
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

  if (!is.null(forced.BC)) {# if forcing in a Basin Characteristic
    #determine column indices of BCs to be forced into regression models
    forced.BC.indx = (1:ncol(expl))[colnames(expl) %in% forced.BC]
  }

  leaps_list = list()
  #loop through quantiles
  for (i in 1:nquants) {
    #create regsubsets object using filled FDCs
    #x=regsubsets(x=expl, y = filled[,i], really.big = T,
    #  weights = comp.wys/mean(comp.wys), nbest=nb, nvmax = nvmax)
    #TMO edits 6/2016, 7/2016
    if (is.null(forced.BC)) {
      x=regsubsets(x=expl, y = filled[,i], really.big = T,
                   weights = comp.wys/mean(comp.wys), nbest=nb, nvmax = nvmax)
    } else {
      # if forcing in a Basin Characteristic
      if (nvmax>1) {
        x=regsubsets(x=expl, y = filled[,i], really.big = T,
                     weights = comp.wys/mean(comp.wys), nbest=nb, nvmax = nvmax,
                     force.in=forced.BC.indx)
      } else {
        #If nvmax=1 do a "fake" call to regsubsets to set up the output file
        x=regsubsets(x=expl, y = filled[,i], really.big = T,
                     weights = comp.wys/mean(comp.wys), nbest=nb, nvmax = nvmax)
      }
    }

    #TRUE/FALSE grid showing BC's and their models
    if (is.null(forced.BC)) {

      w_all=summary(x)$which

    } else {

      # if forcing in a Basin Characteristic
      w_all_temp=summary(x)$which;
      #Since when a variable is forced in, the model containing it does not appear in the regsubsets output,
      #append a first row to w_all for the regression against the forced in variable alone
      #(Implicitly only works for 1 forced in variable)
      first_row<-logical(ncol(expl)+1) #defaults to FALSE
      first_row[which(colnames(w_all_temp)=="(Intercept)")]<-T #entry for intercept
      first_row[which(colnames(w_all_temp)==forced.BC)]<-T #entry for forced in BC
      w_all<-rbind(first_row,w_all_temp); rownames(w_all)<-c(1,rownames(w_all_temp))
      #Remove all but first row of w_all when nvmax==1  but maintain matrix form
      if (nvmax==1) {
        w_all_temp2<-w_all
        w_all<-matrix(nrow=1,ncol=ncol(w_all_temp2),dimnames=list(1,colnames(w_all_temp2)))
        w_all[1,]<-w_all_temp2[1,]
      }
    }

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

        if (length(which(vif(model)>max.VIF)) == 0){
          out = cbind(out, rep('', length(out[,1])))
          out[which(out[,1]%in%names(which(w[[k]] == TRUE))==TRUE),
              length(out[1,])] = coef(model) #coeffs
          out[(which(out[,1]%in%names(which(w[[k]] == TRUE))==TRUE)+1),
              length(out[1,])] = sqrt(diag(vcov(model))) #SE
          out[(which(out[,1]%in%names(which(w[[k]] == TRUE))==TRUE)+2),
              length(out[1,])] = c(" ", vif(model)) #VIF
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

          if (length(which(vif(model)>max.VIF)) == 0){
            out = cbind(out, rep('', length(out[,1])))
            out[which(out[,1]%in%names(which(w[[k]][j,] == TRUE))==TRUE),
                length(out[1,])] = coef(model) #coeffs
            out[(which(out[,1]%in%names(which(w[[k]][j,] == TRUE))==TRUE)+1),
                length(out[1,])] = sqrt(diag(vcov(model))) #SE
            out[(which(out[,1]%in%names(which(w[[k]][j,] == TRUE))==TRUE)+2),
                length(out[1,])] = c(" ", vif(model)) #VIF
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
