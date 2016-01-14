#' Find top models for each flow regime.
#'
#' @description Takes output from \code{\link{compile.vars}} and determines
#' the best average model for each regime.
#'
#' @param top_n_list Output from \code{\link{compile.vars}}
#' @param unfilled_FDCs A matrix of the raw FDC quantiles for each site.
#' This is derived from the output of \code{\link{calcEmpFDCs}}.
#' @param expl A data frame of the potential explanatory variables,
#' derived from the output of \code{\link{getBasinChar}}.
#' @param comp.wys A vector indicating how many complete water years
#' were present for each site.  This is derived from the output of
#' \code{\link{calcEmpFDCs}}.
#' @param save.dir (optional) A directory to which results are to be saved.
#' The default behavior \code{save.dir=NULL} does not save results.
#' @param nvmax (optional)  The maximum number of variables that will be
#' considered in a regression.  The default is \code{6}.
#' @param n (optional) The maximum number of models subsetted.
#' The deafult is \code{20}.
#' @param cens_level (optional)  The value establishing left censoring.
#' The deafult is \code{0.005}.
#' @param zero.val (optional)  The value to which zeroes or negative quantiles
#' will be set to.  The default is \code{0.001}.
#' @param WY.lim (optional)  The minimum number of water years required to be
#'  included in the formation of regional regressions.  Only reference-quality
#'  sites with at least this many complete water years will be used.
#'  The default is \code{10}.
#' @param regimes (optional) A list of three elements (\code{names(regimes) =
#' c('lowflow','medflow','highflow')}), where each element is a character
#' vector indicating the frequencies of each quantile. The default is \code{
#' list(lowflow = c("0.0002","0.0005","0.001","0.002","0.005","0.01","0.02",
#' "0.05","0.1"),medflow = c("0.2", "0.25", "0.3", "0.4", "0.5","0.6", "0.7",
#'  "0.75", "0.8", "0.9"),highflow = c("0.95", "0.98", "0.99", "0.995",
#'  "0.998", "0.999", "0.9995", "0.9998"))}.
#'
#' @details
#' This function outputs the top 3 models for each flow regime by both
#' adjusted R^2 and AIC (for a total of 6 models). It takes the average
#' adjusted R^2 and AIC across quantiles within a flow regime by each number
#' of variables, and whichever number of variables has the highest R^2 or
#' lowest AIC is the "best" type of model for the flow regime.
#'
#' @return Top three models (by AIC and adjustes R2) for each regime.
#'@export
best.models = function(top_n_list,unfilled_FDCs,expl,comp.wys,save.dir=NULL,
  nvmax=6, n=20, cens_level = .005,zero.val=0.001,WY.lim=10,
  regimes = list(lowflow = c("0.0002","0.0005","0.001","0.002","0.005",
    "0.01","0.02","0.05","0.1"),
    medflow = c("0.2", "0.25", "0.3", "0.4", "0.5", "0.6", "0.7", "0.75",
      "0.8", "0.9"),
    highflow = c("0.95", "0.98", "0.99", "0.995", "0.998", "0.999",
      "0.9995", "0.9998"))
) {
  # Function orginially designed by Thomas M. Over and Mike Olsen, 03 July 2015.
  # Modified by William Farmer, 06 July 2015.

  if (is.na(regimes)) {
    regimes.ret <- regimes(unfilled_FDCs,zero_val=zero.val)
    regimes <- regimes.ret$regimes
  }

  raw.input <- list(top_n_list=top_n_list,
    unfilled_FDCs=unfilled_FDCs,expl=expl,comp.wys=comp.wys,save.dir=save.dir,
    nvmax=nvmax, n=n, cens_level = cens_level,
    regimes = regimes,WY.lim=WY.lim,zero.val=zero.val)
  flows=names(regimes)

  dep = unfilled_FDCs
  dep = t(dep)
  nquants = ncol(dep)

  #Read in basin characteristics
  expl=expl[as.numeric(rownames(expl)) %in%  as.numeric(rownames(dep)),]
  class=as.character(expl$CLASS)
  #keep only reference gages with >= 10 complete water years
  screen.ndx <- comp.wys>=WY.lim & class=='Ref'
  dep=dep[screen.ndx,]
  # WHF, 07/06/2015: I'm not crazy about this hard-coded pre-processing.

  #if NA, NaN, or negative, set to .001 (log10(.001)=-3)
  dep[dep <= 0 | is.na(dep) | is.nan(dep)]=zero.val

  #log10 transform dependent variables and censoring level
  dep=log10(dep)
  cens_level=log10(cens_level)

  #apply same subsetting to explanatory variables
  expl=expl[screen.ndx,-1]
  #remove class variable
  expl=as.matrix(expl)
  comp.wys=comp.wys[screen.ndx]

  #create empty lists and matrices
  best=list(); length(best)=3
  best_by_flow=list(); length(best_by_flow)=6
  adjr2_quants = aic_quants = matrix(nrow=3, ncol=27)
  best_nvars_adjr2 = best_nvars_aic = c()
  length(best_nvars_adjr2) = length(best_nvars_aic) = 3

  #loop through flow regimes
  for(i in 1:length(flows)){
    #set size of adjusted r^2 and aic matrices
    adjr2_mat=aic_mat=matrix(nrow=nvmax, ncol=length(regimes[[i]]))
    if(flows[[i]] == "lowflow") {
      js = c(1:9)
    } else if (flows[[i]] == "medflow") {
      js = c(10:19)
    } else {
      js = c(20:27)
    }
    for(j in 1:nvmax){
      if (j > length(raw.input$top_n_list[[i]])) {next}
      if (is.null(raw.input$top_n_list[[i]][[j]])) {next}
      if (is.na(raw.input$top_n_list[[i]][[j]])) {next}
      #read in R^2 and aic values from compare.censReg output
      temp <- comp(flow=flows[i], numvars=j,
        unfilled_FDCs=raw.input$unfilled_FDCs,
        comp.wys=raw.input$comp.wys,top_n_list=raw.input$top_n_list,
        expl=raw.input$expl,WY.lim=raw.input$WY.lim,
        zero.val=raw.input$zero.val,cens_level = raw.input$cens_level,
        n=raw.input$n, plot.dir=NULL)

      adjr2 = temp$r2adjg
      n_models=length(adjr2)/length(regimes[[i]])
      aic = temp$aicg
      #loop through each quantile to fill in top 3 r^2's
      # and aic's in r^2 and aic matrix
      if (n_models<1) {next}
      if(n > n_models) top=n_models else top=n
      for(k in 1:length(regimes[[i]])){
        adjr2_mat[j,k]= adjr2[seq(k,top*length(regimes[[flows[[i]]]]),
          length(regimes[[flows[[i]]]]))][1]
        aic_mat[j,k] = aic[seq(k,top*length(regimes[[flows[[i]]]]),
          length(regimes[[flows[[i]]]]))][1]
      }
    }
    #apply mean across quantiles
    mean_adjr2 = apply(adjr2_mat, 1, mean)
    mean_aic = apply(aic_mat, 1, mean)
    #find number of variables with best average GOF stats
    best_nvars_adjr2[i] = which(mean_adjr2==max(mean_adjr2,na.rm=T))
    best_nvars_aic[i] = which(mean_aic==min(mean_aic,na.rm=T))
    #find "best" models in list from compile.vars
    nc <- min(ncol(as.matrix(top_n_list[[i]][[best_nvars_adjr2[i]]])),4)
    adjr2_vars=top_n_list[[i]][[best_nvars_adjr2[i]]][-1,1:nc]
    nc <- min(ncol(as.matrix(top_n_list[[i]][[best_nvars_aic[i]]])),4)
    aic_vars=top_n_list[[i]][[best_nvars_aic[i]]][-1,1:nc]
    #create a list containing the coefficients, sd's, VIF's, r^2, and
    #  aic for each of the top 6 models, read from compare.censReg output
    #fill in  the second set of R^2 and aic matrices across quantiles for
    #  the larger overall plots
    temp.ar2 <- comp(flow= flows[[i]], numvars=best_nvars_adjr2[i],
      unfilled_FDCs=raw.input$unfilled_FDCs,
      comp.wys=raw.input$comp.wys,top_n_list=raw.input$top_n_list,
      expl=raw.input$expl,WY.lim=raw.input$WY.lim,
      zero.val=raw.input$zero.val,cens_level = raw.input$cens_level,
      n=raw.input$n, plot.dir=NULL)
    temp.aic <- comp(flow= flows[[i]], numvars=best_nvars_aic[i],
      unfilled_FDCs=raw.input$unfilled_FDCs,
      comp.wys=raw.input$comp.wys,top_n_list=raw.input$top_n_list,
      expl=raw.input$expl,WY.lim=raw.input$WY.lim,
      zero.val=raw.input$zero.val,cens_level = raw.input$cens_level,
      n=raw.input$n, plot.dir=NULL)
    for(l in 1:6){
      if(l < 4){
        if (l>length(temp.ar2$summaries)) {
          best_by_flow[[l]] <- NA
          adjr2_quants[l,js] <- NA
          next
        }
        #top 3 adjr2 models for each flow regime (using first quantile)
        best_by_flow[[l]] =
          temp.ar2$summaries[[l]][-c(3,(((best_nvars_adjr2[i]+1)*3)+2)),]
        #colnames(best_by_flow[[l]])=substr(colnames(best_by_flow[[l]]), 2, 10)
        colnames(best_by_flow[[l]])=temp.ar2$summaries.names[[l]]
        adjr2_quants[l,js]=
          as.numeric(best_by_flow[[l]][nrow(best_by_flow[[l]]),][-1])
        best_by_flow[[l]] = temp.ar2$summaries[[l]][-3,]
      }
      if(l > 3){
        if ((l-3)>length(temp.aic$summaries)) {
          best_by_flow[[l]] <- NA
          aic_quants[l-3,js] <- NA
          next
        }
        #top 3 aic models for each flow regime (using first quantile)
        best_by_flow[[l]] =
          temp.aic$summaries[[l-3]][-c(3,(((best_nvars_aic[i]+1)*3)+1)),]
        colnames(best_by_flow[[l]])=temp.aic$summaries.names[[l-3]]
        aic_quants[l-3,js]=
          as.numeric(best_by_flow[[l]][nrow(best_by_flow[[l]]),][-1])
        best_by_flow[[l]] = temp.aic$summaries[[l-3]][-3,]
      }
    }
    #names the models in the list by their GOF statstic and model rank
    names(best_by_flow)=c("Adj R2 1", "Adj R2 2", "Adj R2 3",  "AIC 1",
      "AIC 2", "AIC 3")
    best[[i]]=best_by_flow
  }
  names(best)=flows

  #Plot top 3 adjusted R^2's and top 3 AIC's across quantiles
  if(!is.null(save.dir)){
    pdf(file.path(save.dir,"ADJR2_across_quants.pdf"), paper= "USr",
      width=11, height=17)
    plot(adjr2_quants[1,], col='red', type='l',
      xlab="Non-Exceedance Probabilities", ylab="Adj R^2", xaxt='n',
      ylim=range(adjr2_quants)+c(-sd(adjr2_quants), sd(adjr2_quants)),
      main="Top 3 Adj R^2 Models",
      sub=paste("lowflow # of variables = ", best_nvars_adjr2[1],
        " | medflow # of variables = ",
        best_nvars_adjr2[2], " | highflow # of variables = ",
        best_nvars_adjr2[3]), cex.sub=.75)
    lines(adjr2_quants[2,], col='blue')
    lines(adjr2_quants[3,], col='green')
    axis(1, at=1:27, labels=as.character(unlist(regimes)),
      tick=T, cex.axis=.55)
    legend("topright", legend=c(1:3), col=c('red','blue','green'),
      cex=.6, lty=rep(1,3), lwd=rep(2,3), title="Model Rank")
    dev.off()

    pdf(file.path(save.dir,"AIC_across_quants.pdf"), paper= "USr",
      width=11, height=17)
    plot(aic_quants[1,], col='red', type='l',
      xlab="Non-Exceedance Probabilities", ylab="AIC", xaxt='n',
      ylim=range(aic_quants)+c(-sd(aic_quants), sd(aic_quants)),
      main="Top 3 AIC Models",
      sub=paste("lowflow # of variables = ", best_nvars_aic[1],
        " | medflow # of variables = ",
        best_nvars_aic[2], " | highflow # of variables = ",
        best_nvars_aic[3]), cex.sub=.75)
    lines(aic_quants[2,], col='blue')
    lines(aic_quants[3,], col='green')
    axis(1, at=1:27, labels=as.character(unlist(regimes)),
      tick=T, cex.axis=.55)
    legend("topright", legend=c(1:3), col=c('red','blue','green'),
      cex=.6, lty=rep(1,3), lwd=rep(2,3), title="Model Rank")
    dev.off()
    save(list=c("best","aic_quants","adjr2_quants"),
      file=file.path(save.dir("bestmodels.RData")))
  }
  result <- list(best=best,adjr2_quants=adjr2_quants,aic_quants=aic_quants,
    regimes=regimes)
  return(result)
}
