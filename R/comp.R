#' Compare the models with a specifc number of variables.
#'
#' @description
#' A sub-function for single comparison of models in a flow regime.
#'
#' @param flow A character string indicating the name of the flow regime.
#' @param numvars The maximum number of variables considered in
#' the regression.
#' @param unfilled_FDCs A matrix of the raw FDC quantiles for each site.
#' This is derived from the output of \link{\code{calcEmpFDCs}}.
#' @param comp.wys A vector indicating how many complete water years
#' were present for each site.  This is derived from the output of
#' \link{\code{calcEmpFDCs}}.
#' @param top_n_list Output from \link{\code{compile.vars}}
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
#' @param n (optional) The maximum number of models subsetted.
#' The deafult is \code{20}.
#' @param plot.dir A directory to which plots are to be saved.
#' The default behavior \code{save.dir=NULL} does not create plots.
#'
#' @details
#' Lorem ipsum...
#'
#' @return Unclear.  Used by \link{\code{best.models}}, among others.
#'@export
comp<-function(flow, numvars, unfilled_FDCs,comp.wys,top_n_list,expl,WY.lim=10,
  zero.val=0.001,cens_level = 0.005, n=20,plot.dir=NULL) {
  # Function orginially designed by Thomas M. Over and Mike Olsen, 03 July 2015.
  # Modified by William Farmer, 06 July 2015.
  #     Re-written to prevent writing key files; returns output instead.

  make.plots <- FALSE
  if (!is.null(plot.dir)) {make.plots <- TRUE}

  dep = unfilled_FDCs
  dep = t(dep)
  nquants = ncol(dep)

  #Read in basin characteristics
  expl=expl[as.numeric(rownames(expl)) %in%  as.numeric(rownames(dep)),]
  class=as.character(expl$CLASS)
  #keep only reference gages with >= 10 complete water years
  screen.ndx <- comp.wys>=WY.lim & class=='Ref'
  dep=dep[screen.ndx,]

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

  flows=c('lowflow', 'medflow', 'highflow')
  #read in best n models
  vars = top_n_list[[which(flows==flow)]][[numvars]]
  if (is.null(vars)){
    result <- NULL
    return(result)
  }
  if(flow == "lowflow") {
    js = c(1:9)
  } else if (flow == "medflow") {
    js = c(10:19)
  } else {
    js = c(20:27)
  }
  #set up summary statistics arrays
  r2=array(0,dim=c(2,numvars)); r2g=rep(0,2)
  aic=array(0,dim=dim(r2)); aicg=rep(0,2)
  sig=array(0,dim=dim(r2)); sigg=rep(0,2)
  r2adj=array(0,dim=dim(r2)); r2adjg=rep(0,2)
  evfreq = array(0, dim = c((length(vars[,1])-3), 2))
  #output2 is summary of coefficients
  output2 = array(" ", dim = c(20*(as.numeric(numvars)+2), (length(js)+1)))
  #zoplot is binary plot of whether or not a variable was used
  zoplot = array(0, dim = c((length(vars[,1])-3),22)) #WHF: Why hard-code 22?
  zoplot[,1] = evfreq[,1] =
    as.character(vars[,1])[-((length(vars[,1])-2):length(vars[,1]))]
  #frequency of EV used
  evs = c()
  result <- names.result <- list()
  for(i in 1:(ncol(vars)-1)){
    if (sum(vars[,i+1]=="*",na.rm=T)==0) {
      output <- NULL
      next
    }
    evfreq[which(evfreq%in%as.character(vars[which(vars[,i+1]=="*"),1])==TRUE),
      2] = as.numeric(evfreq[which(evfreq%in%as.character(
        vars[which(vars[,i+1]=="*"),1])==TRUE),2]) + 1
    evs[(length(evs)+1):(length(evs)+as.numeric(numvars))] =
      as.character(vars[which(vars[,i+1]=="*"),1])
    zoplot[which(evfreq%in%as.character(vars[which(vars[,i+1]=="*"),
      1])==TRUE),(i+1)] = as.numeric(1)
    #coefficients & summary stats for each model
    output = array(" ", dim = c((3*((1+as.numeric(numvars)))+2),
      (1+length(js))))
    output[((1:((1+as.numeric(numvars)))*3)-2),1] = c("(Intercept)",
      as.character(vars[which(vars[,i+1]=="*"),1]))
    output[((1:((1+as.numeric(numvars)))*3)-1),1] = "SD"
    output[((1:((1+as.numeric(numvars)))*3)),1] = "VIF"
    output[(length(output[,1])-1):length(output[,1]),1] = c("R^2", "AIC")
    #loop through quantiles in flow regime
    for(j in 1:length(js)){
      #create censored model
      if(numvars == 1){
        cens_mod = censReg(as.lcens(dep[,js[j]],
          cens_level) ~expl[,colnames(expl)%in%vars[which(vars[,i+1]=="*"),
            1]==TRUE],weights = comp.wys/mean(comp.wys))
      } else{
        cens_mod = censReg(as.lcens(dep[,js[j]], cens_level) ~.,
          as.data.frame(expl[,colnames(expl)%in%vars[which(vars[,i+1]=="*"),
            1]==TRUE]),
          weights = comp.wys/mean(comp.wys))
      }
      #coefficient values
      output[((1:(1+as.numeric(numvars))*3)-2),(j+1)] =
        as.numeric(coefficients(cens_mod))
      output2[(2+(i-1)*(as.numeric(numvars)+2)):((2+as.numeric(numvars))*i),
        (j+1)] = as.numeric(coefficients(cens_mod))
      #Standard Deviations
      output[((1:((1+as.numeric(numvars)))*3)-1),(j+1)] =
        cens_mod$STDDEV[1:(as.numeric(numvars)+1)]
      #VIF
      output[((1:((1+as.numeric(numvars)))*3)),(j+1)] =
        if(numvars!=1) {c(" ", smwrStats::vif(cens_mod))} else {0}
      #R^2
      output[(length(output[,1])-1),(j+1)] = summary(cens_mod)$R2
      #AIC
      output[length(output[,1]),(j+1)] = cens_mod$AIC
      r2g = cbind(r2g, c(i, summary(cens_mod)$R2))
      #couldn't figure out how to extract sigma, so this calculates it
      sigg = cbind(sigg, c(i, sqrt(sum(residuals(cens_mod)^2)/
          cens_mod$survreg$df.residual)))
      aicg = cbind(aicg, c(i, cens_mod$AIC))
      #calculate adjusted r^2
      r2adjg = cbind(r2adjg, c(i, 1-((1-summary(cens_mod)$R2)*(nrow(dep)-1)/
          (nrow(dep)-(as.numeric(numvars))-1))))

      #Diagnostic Plots
      stzd_res= residuals(cens_mod)/(sqrt((sum(residuals(cens_mod)^2)/
          cens_mod$survreg$df.residual)*
          (1-summary(cens_mod)$diagstats$leverage)))
      if (make.plots) {
        if (j==1) {
          #plot diagnostics for each quantile
          pdf(file.path(plot.dir,paste("model", i,
            "_quantile_diagnostics.pdf", sep = "")), onefile=T)
          par(mfrow = c(2,2))
        }
        #Residuals vs. fitted
        plot(fitted(cens_mod), residuals(cens_mod), xlab="Fitted Values",
          ylab="Residuals", main=paste("Quantile ", colnames(dep)[js[j]],
            " | Residuals vs. Fitted", sep=''))
        lines(lowess(fitted(cens_mod), residuals(cens_mod)), col='red')
        #QQ
        qqnorm(stzd_res)
        #Scale-Location
        plot(fitted(cens_mod), sqrt(abs(stzd_res)), xlab="Fitted Values",
          ylab=expression(sqrt(abs("Standardized Residuals"))),
          main=paste("Quantile ", colnames(dep)[js[j]]," | Scale-Location",
            sep=''))
        lines(lowess(fitted(cens_mod), sqrt(abs(stzd_res))), col='red')
        #Res vs. Lev
        plot(summary(cens_mod)$diagstats$leverage, stzd_res, xlab="Leverage",
          ylab="Standardized Residuals", main=paste("Quantile ",
            colnames(dep)[js[j]]," | Residuals vs. Leverage", sep=''))
        lines(lowess(summary(cens_mod)$diagstats$leverage, stzd_res),
          col='red')
      }
    }
    if (make.plots) {
      dev.off()
    }
    var_names= c("Intercept", colnames(expl)[
      colnames(expl)%in%vars[which(vars[,i+1]=="*"),1]==TRUE])
    output2[(2+(i-1)*(as.numeric(numvars)+2)):((2+as.numeric(numvars))*i),
      1] = var_names
    result[[i]] <- output
    names.result[[i]] <- c(" ", colnames(dep)[(js)])
    #plot coefficients across quantiles
    if (make.plots) {
      pdf(file.path(plot.dir,paste("model", i, "_coefs.pdf", sep = "")),
        onefile=T)
      par(mfrow = c(min((as.numeric(numvars)+1),3),1))
      par(mar = c(4.4,  4.1, 1.5, .5))
      for(p in 1:(1+as.numeric(numvars))){
        plot(as.numeric(output[((p-1)*3+1),(2:length(output[1,]))]),
          xlab = "Quantile",
          ylab = "Coefficient Value", pch = 16, cex = .5, xaxt = 'n', type='b',
          ylim=range(as.numeric(output[((p-1)*3+1),(2:length(output[1,]))]))+
            1.5*c(-min(as.numeric(output[((p-1)*3+2),(2:length(output[1,]))])),
              max((as.numeric(output[((p-1)*3+2),(2:length(output[1,]))])))))
        axis(1, at = (js - min(js) +1), as.numeric(colnames(dep))[js]/100)
        title(output[((p-1)*3+1),1])
        #SE lines
        lines(as.numeric(output[((p-1)*3+1),(2:length(output[1,]))])+
            as.numeric(output[((p-1)*3+2),(2:length(output[1,]))]), lty=2)
        lines(as.numeric(output[((p-1)*3+1),(2:length(output[1,]))])-
            as.numeric(output[((p-1)*3+2),(2:length(output[1,]))]), lty=2)
      }
      dev.off()
      #for higher number of variables, make text for pairs plot smaller
      if(as.numeric(numvars) == 6){ cl = .56
      } else if(as.numeric(numvars) == 5) { cl = .8
      } else {cl = 1}
      #pairs plot
      if(numvars!=1){
        pdf(paste("model", i, "_pairs.pdf", sep = ""))
        pairs(expl[,which(colnames(expl)%in%var_names)], label =
            paste((var_names)[-1], '\nVIF: ',
              round(smwrStats::vif(cens_mod),4)), cex.labels = cl)
        dev.off()
      }
    }
  }

  if (identical(r2g,rep(0,2))) {
    r2g <- sigg <- aicg <- r2adjg <- matrix(NA,nrow=2,ncol=1)
  } else {
    r2g = r2g[,-1]; sigg = sigg[,-1]; aicg = aicg[,-1]; r2adjg = r2adjg[,-1]
  }
  zoplot[,22] = evfreq[,2]
  colnames(output2) <- c(" ",colnames(dep)[(js)])
  big.result <- list(summaries=result, summaries.names=names.result,
    coefs=output2, r2g=r2g[2,], r2adjg=r2adjg[2,], sigg=sigg[2,],
    aicg=aicg[2,], zoplot=zoplot)

  if (make.plots) {
    pdf(paste(numvars, "vars_", flow, "_r2.pdf", sep = ""))
    plot.scores(n = n, flow = flow, numvars = numvars, names = colnames(dep),
      value = r2g, type= "R^2", js = js)
    dev.off()

    pdf(paste(numvars, "vars_", flow, "_adj_r2.pdf", sep = ""))
    plot.scores(n = n, flow = flow, numvars = numvars, names = colnames(dep),
      value = r2adjg, type = "Adj R^2", js = js)
    dev.off()

    pdf(paste(numvars, "vars_", flow, "_sigma.pdf", sep = ""))
    plot.scores(n = n, flow=flow, numvars=numvars, names=colnames(dep),
      value = sigg, type = "Sigma", js = js)
    dev.off()

    pdf(paste(numvars, "vars_", flow, "_AIC.pdf", sep = ""))
    plot.scores(n = n, flow=flow, numvars=numvars, names=colnames(dep),
      value = aicg, type = "AIC", js = js)
    dev.off()

    pdf(paste(numvars, "vars_", flow, "_freqev.pdf", sep = ""))
    par(mar = c(7,4.1,4.1,2.1))
    barplot(table(evs), las = 2, cex.names = .5, ylab="Frequency")
    title(paste(numvars, " variables ", flow, " Explanatory Variables",
      sep = ""))
    dev.off()
  }
  return(big.result)
}
