
# CODE BELOW HERE ----
# Has not been fully developed and integrated into the main Airdrop.

# compare.censReg ----
#GENERAL INFO
#This function compares models across the qauntiles within a flow regime for each number of variables by their ranks
#(models re-ranked by Pseudo R^2 from censReg)
#The first function takes the output from compiles.vars and creates summary statistics and diagnostic plots that help determine which models have the best GOF stats
#The second function simply implements the first over each flow regime and each number of variables
#The outputs are as follows: plots and tables showing R^2, Adj R^2, Sigma, and AIC for each model, frequency statsistcs for the BC's, and summary
#statistics for each of the censored regression models
#The outputs can be found within the "compare.censReg" folder in the base directory (base.dir)
#Written by Mike Olson on 7/1/2015 at USGS ILWSC
#Modified from 'Z:\SE_Scenarios\full.199.1981-2010\daily\se_basin_chars3\updates8\leaps.VIF3.WHF_NoNHD.NoRoads.NoPop.NoDev\Rscripts\comparelm.R'

#PARAMETERS (for the second function, not first)
#flows is a vector containing the names for each of the flow regimes, only one at a time is used within the first function
#nvmax is the number of variable levels the models will have
#base.dir is the location where the "compare.censReg" folder will be, which contains all the output from this function
#unfilled_FDCs_path
#expl_path
#PORstats_path
#cens_level
#top_n_list_path
#n is number of top models, it should be the same n as for compare.censReg

#EXAMPLE CALLS
#compare.censReg()


compare.censReg=function(flows=c("lowflow", "medflow", "highflow"), nvmax=6,
  save.dir, unfilled_FDCs,expl,comp.wys,top_n_list,WY.lim,
  zero.val=0.001,cens_level = 0.005, n=20) {
  # Function orginially designed by Thomas M. Over and Mike Olsen, 03 July 2015.
  # Modified by William Farmer, 06 July 2015.
  #   Appears to just wrap a sub function for producing figures.
  for(i in 1:length(flows)){
    for(j in 1:nvmax){
      temp <- comp(flow=flows[i], numvars=j, unfilled_FDCs=unfilled_FDCs,
        expl=expl, comp.wys=comp.wys, top_n_list=top_n_list, plot.dir=save.dir,
        zero.val=zero.val,cens_level=cens_level, n=n,WY.lim=WY.lim)
    }
  }
}



plot.scores<-function(value, type, js, names, numvars, flow, n){
  if (n > max(value[1,])) top=max(value[1,]) else top=n
  col_ramp=rainbow(length(js), start=.1, end=.9)
  plot(value[1,(1:top)*length(js)], value[2,(1:top)*length(js)], type = "S",
    xlim=c(.5,n), ylim=range(value[2,])+c(-sd(value[2,]), sd(value[2,])),
    xlab = "Model Number", ylab = type, col=col_ramp[length(js)])
  points(value[1,(1:top)*length(js)], value[2,(1:top)*length(js)], pch = 16,
    cex = .5, col=col_ramp[length(js)])
  text((value[1,length(js)]-.5), value[2,length(js)], names[max(js)], cex=.6,
    col=col_ramp[length(js)])
  for(i in 1:(length(js))-1){
    points(value[1,((1:top)*length(js)-i)], value[2,((1:top)*length(js)-i)],
      type = "S", col=col_ramp[length(js)-i])
    points(value[1,((1:top)*length(js)-i)], value[2,((1:top)*length(js)-i)],
      pch = 16, cex = .5, col=col_ramp[length(js)-i])
    text((value[1,(length(js)-i)]-.5), value[2,(length(js)-i)],
      names[max(js)-i], cex=.6, col=col_ramp[length(js)-i])
  }
  legend("top", names[js], col=col_ramp, cex=.6, lty=rep(1, length(js)),
    lwd=rep(2, length(js)), horiz=T)
  if(numvars > 1) {title(paste(numvars, "Variables", flow, type,  sep = " "),
    sub=paste("Max VIF = 5 | Top ", top, " Models"), cex.sub=.75)}
  if(numvars == 1) {title(paste(numvars, "Variable", flow, type,  sep = " "),
    sub=paste("Max VIF = 5 | Top ", top, " Models"), cex.sub=.75)}
}

# freqev.plot ----
#GENERAL INFO
#This function creates a color coded plot for each flow regime, detailing the frequency of each BC in all the models for each number of variables
#Written by Mike Olson on 7/2/2015 at USGS ILWSC
#Modified from 'Z:\SE_Scenarios\full.199.1981-2010\daily\se_basin_chars3\updates8\leaps.VIF3.WHF_NoNHD.NoRoads.NoPop.NoDev\Rscripts\imageplotFreqofEv.TMOmod.R'

#PARAMETERS
#base.dir the location where the "csv's will be output"compare.censReg" folder will be, which contains the rest of the output
#nvmax is the number of variable levels the models will have
#flows is a vector containing the names for each of the flow regimes, only one at a time is used within the first function
#n is number of top models, it should be the same n as for compile.vars and compare.censReg

#EXAMPLE CALLS
#freqev.plot()

freqev.plot=function(base.dir="Z:/UpperMissouri/MikeOlson/My_Flow_Data/leaps",
  nvmax=6, flows=c("lowflow", "medflow", "highflow"), n=20){
  setwd(base.dir)
  dir.create("freqev.plot", showWarnings=F)
  setwd(paste(base.dir, "/freqev.plot", sep=''))
  #source("Z:/UpperMissouri/MikeOlson/My_Flow_Data/
  #  leaps functions/freqev.plot/sub_functions.R")
  for(i in 1:3) {
    #BC names
    labs = as.character(read.csv(paste(base.dir, '/compare.censReg/Frequency/',
      (nvmax-1), 'vars_', flows[i], '_freqplot.csv', sep=''),
      header = TRUE)[,1])
    #frequency data for each BC from compare.censReg output
    freq = as.numeric(read.csv(paste(base.dir, '/compare.censReg/Frequency/',
      nvmax, 'vars_', flows[i], '_freqplot.csv', sep=''), header = TRUE)[,n+2])
    #loops through remaining variables levels and combines frequency
    #  data across number of variables
    for(j in (nvmax-1):1){
      freq = cbind(freq, as.numeric(read.csv(paste(base.dir,
        '/compare.censReg/Frequency/', j, "vars_",
        flows[i], "_freqplot.csv", sep = ""), header = TRUE)[,n+2]))
    }
    #remove observations with 0 frequency
    k = 1
    while(k<length(freq[,1])){
      if(sum(as.numeric(freq[k,]))==0){
        freq = freq[-k,]
        labs = labs[-k]
      }
      else { k = k+1 }
    }
    #plot frequency data
    pdf(paste( flows[i], " Freq.pdf"), paper='USr', height=8.5,
      width=round(length(names)/5), onefile=T)
    freqev_ImagePlot(t(freq), xLabels = labs, yLabels = c(nvmax:1),
      title = c(paste("Freq of EV", flows[i])))
    dev.off()
  }
}

freqev_ImagePlot <- function(x, ...){
  library(Hmisc)
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }

  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(5.5,1), heights=c(1,1))

  # Colormap

  ColorRamp <- rainbow(21, start = .1, end = .9)
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]

  # Data Map
  par(mar = c(10,4,2.5,2))
  #color-coded image
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
    ylab="Number of Variables in Model", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  #gridlines
  grid(nx = length(xLabels), ny = 0, col = "black", lty = "solid")
  #custom axes
  mgp.axis(side=1, at=1:length(xLabels), mgp=c(3,0,0), labels=xLabels,
    cex.axis=0.75, las = 2, tick=F)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
    cex.axis=0.8)

  # Color Scale on right hand side of plot
  par(mar = c(3,2,2.5,2))
  image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels),
    nrow=1), col=ColorRamp, xlab="",ylab="", xaxt="n", cex.main=.75)

  layout(1)
}

# flow.numvars.plot ----
#GENERAL INFO
#This function creates a grid detailing all the BC's in the top n models for each flow regime and number of variables
#Written by Mike Olson on 7/2/2015 at USGS ILWSC
#Modified from 'Z:\SE_Scenarios\full.199.1981-2010\daily\se_basin_chars3\updates8\leaps.VIF3.WHF_NoNHD.NoRoads.NoPop.NoDev\Rscripts\imagePlotbyflowandnumvars.TMOmod.R'

#PARAMETERS
#base.dir the location where the "csv's will be output"compare.censReg" folder will be, which contains the rest of the output
#nvmax is the number of variable levels the models will have
#flows is a vector containing the names for each of the flow regimes, only one at a time is used within the first function
#n is number of top models, it should be the same n as for compile.vars and compare.censReg

#EXAMPLE CALLS
#freqev.plot()

flow.numvars.plot=function(
  base.dir="Z:/UpperMissouri/MikeOlson/My_Flow_Data/leaps", nvmax=6,
  flows=c("lowflow", "medflow", "highflow"), n=20){
  setwd(base.dir)
  dir.create("flow.numvars.plot", showWarnings=F)
  setwd(paste(base.dir, "/flow.numvars.plot", sep=''))
  #source("Z:/UpperMissouri/MikeOlson/My_Flow_Data/
  #  leaps functions/flow.numvars.plot/sub_functions.R")
  #loop through flow regimes
  for(i in 1:3) {
    #loop through variables levels
    for(j in 1:nvmax) {
      #read in frequency data for BC's from compare.censReg output
      freq = read.csv(paste(base.dir, '/compare.censReg/Frequency/', j,
        "vars_", flows[i], "_freqplot.csv", sep = ""), header = TRUE)
      k = 1
      #remove observations with 0 frequency
      while(k<length(freq[,1])) {
        if(sum(as.numeric(freq[k,2:(n+1)]))==0) {
          freq = freq[-k,]
        }
        else { k = k+1 }
      }
      #remove header rows and set labels
      lab = freq[,1]; freq=freq[,-1]; freq=freq[,-(n+1)]
      #plot graphs
      pdf(paste(j, "vars", flows[i], "Flow Freq.pdf"))
      flow_numvars_ImagePlot(t(freq), xLabels = as.character(lab),
        yLabels = c(1:n), title = c(paste("Freq of EV", flows[i], j, "Vars")))
      dev.off()
    }
  }
  setwd(base.dir)
}

flow_numvars_ImagePlot <- function(x, ...){
  min <- min(x); max <- max(x)
  yLabels <- rownames(x); xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]; max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ) yLabels <- c(Lst$yLabels)
    if( !is.null(Lst$xLabels) ) xLabels <- c(Lst$xLabels)
    if( !is.null(Lst$title) ) title <- Lst$title
  }
  # check for null values
  if( is.null(xLabels) ) xLabels <- c(1:ncol(x))
  if( is.null(yLabels) ) yLabels <- c(1:nrow(x))

  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))

  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0,1,length=256),  # Red
    seq(0,1,length=256),  # Green
    seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))

  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]

  # Data Map
  par(mar = c(8,4,2.5,0.5))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=((0:32)/32), xlab="",
    ylab="Model Number (black means included in model)", axes=FALSE,
    zlim=c(min,max))
  if( !is.null(title) ) title(main=title)

  grid(nx = length(xLabels), ny = 0, col = "red")
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.5, las = 2)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
    cex.axis=0.7)
  layout(1)
}

# byvars ----
#GENERAL INFO
#This function creates graphs for GOF statistics across the top n models for each quantile
#Right now it gets the GOF by calling censReg, but should be edited to read in output from censReg at some point so it goes faster
#Written by Mike Olson on 7/2/2015 at USGS ILWSC
#Modified from 'Z:\SE_Scenarios\full.199.1981-2010\daily\se_basin_chars3\updates8\leaps.VIF3.WHF_NoNHD.NoRoads.NoPop.NoDev\Rscripts\comparelm.R'

#PARAMETERS
#n is number of top models, it should be the same n as for compile.vars and compare.censReg
#cens_level is the censoring level used in censReg and should be inputted in the same units as the FDC's
#with1.flag is a logical parameter than determines whether the plots should include the 1 variable models
#base.dir is the location where the "byvars" folder is created, which is where pdf's will be saved
#nvmax is the number of variable levels the models will have
#gof is a vector of character values detailing the GOF statistics to be plotted
#names is a vector of non-exceedance probilities to be used in labeling the plots
#unfilled_FDCs_path is the location of the unfilled in FDCs .RData file created from the compute.FDCs function
#PORstats_path is the location of water year statistics for each gage created from the compute.FDCs function
#expl_path is the location of the BCs
#top_n_list_path is the location of the list of lists RData object outputted from compile.vars

#EXAMPLE CALLS
#byvars()

byvars=function(n=20, cens_level=.005, with1.flag=T,
  base.dir="Z:/UpperMissouri/MikeOlson/My_Flow_Data/leaps",
  gof=c("r2","r2adj","sigma","aic"),
  names=c(0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,
    0.05,0.10,0.20,0.25,0.30,0.40,0.50,
    0.60,0.70,0.75,0.80,0.90,0.95,0.98,0.99,0.995,0.998,0.999,0.9995,0.9998),
  unfilled_FDCs_path=paste0('Z:/UpperMissouri/MikeOlson/My_Flow_Data/output/',
    'FDC output/FDCs.cens_nofill.RData'),
  PORstats_path=paste0("Z:/UpperMissouri/MikeOlson/My_Flow_Data/",
    "output/PORstats.compWyr.csv"), numvars=6,
  expl_path=paste0("Z:/UpperMissouri/BasinCharacteristics/keep_cod6.max_",
    "oneval_frac0.5.max_divide/finalBCs.keep_code6.max_oneval_frac0.5.csv"),
  top_n_list_path = paste0("Z:/UpperMissouri/MikeOlson/My_Flow_Data/",
    "leaps/compile.vars/top_n_list.RData")) {

  #source("Z:/UpperMissouri/MikeOlson/My_Flow_Data/
  #  leaps functions/byvars/sub_functions.R")

  library(smwrQW)
  #create empty arrays
  r2=array(0,dim=c(n,numvars)); r2g=rep(0,n)
  aic=array(0,dim=dim(r2)); aicg=rep(0,n)
  sig=array(0,dim=dim(r2)); sigg=rep(0,n)
  r2adj=array(0,dim=dim(r2)); r2adjg=rep(0,n)

  load(unfilled_FDCs_path)
  load(top_n_list_path)

  #read in file with complete water years data
  PORstats = read.csv(PORstats_path,row.names=1)
  comp.wys=PORstats$nWys.compWyr

  dep = unfilled_FDCs
  dep = t(dep)
  nquants = ncol(dep)

  #Read in basin characteristics
  expl = read.csv(expl_path,row.names=1,
    colClasses=c("character","factor",rep("numeric",(62))))
  expl=expl[as.numeric(rownames(expl)) %in%  as.numeric(rownames(dep)),]
  class=as.character(expl$CLASS)
  #keep only reference gages with >= 10 complete water years
  dep=dep[comp.wys>=10 & class=='Ref',]

  #if NA, NaN, or negative, set to .001 (log10(.001)=-3)
  dep[dep <= 0 | is.na(dep) | is.nan(dep)]=.001

  #log10 transform dependent variables and censoring level
  dep=log10(dep)
  cens_level=log10(cens_level)

  #apply same subsetting to explanatory variables
  expl=expl[comp.wys>=10 & class=='Ref',-1]
  #remove class variable
  expl=as.matrix(expl)
  comp.wys=comp.wys[comp.wys>=10 & class=='Ref']

  #loop through gof statistics
  for(g in gof) {
    setwd(base.dir)
    dir.create("byvars", showWarnings=F)
    setwd(paste(base.dir, "/byvars", sep=''))
    if (with1.flag) {
      pdf(paste(g,"_byvars_with1vars.pdf",sep=""),onefile=T)
    } else {
      pdf(paste(g,"_byvars.pdf",sep=""), onefile=T)
    }
    flows=c("lowflow","medflow","highflow")
    #loop through flow regimes
    for(flow in c("lowflow","medflow","highflow")) {
      if (flow == "lowflow") {
        js = c(1:9)
      } else if (flow == "medflow") {
        js = c(10:19)
      } else {
        js = c(20:27)
      }
      #loop through quantiles in that flow regime
      for(j in 1:length(js)) {
        if (with1.flag) min.nvars = 1 else min.nvars = 2
        #loop through number of variables
        for(k in min.nvars:numvars){
          #single regression uses slightly different formula
          #  than multiple regression
          if(k == 1){
            for(i in 1:n) {
              #find corresponding entry in compile.vars list
              vars = top_n_list[[which(flow==flows)]][[k]]
              #create model
              cens_mod = censReg(as.lcens(dep[,js[j]],
                cens_level) ~ expl[,colnames(expl)%in%vars[
                  which(vars[,i+1]=="*"),1]==TRUE],
                weights = comp.wys/mean(comp.wys))
              #fill in GOF arrays
              r2g[i]=summary(cens_mod)$R2
              sigg[i]=sqrt(sum(residuals(cens_mod)^2)/
                  cens_mod$survreg$df.residual)
              r2adjg[i]=1-((1-summary(cens_mod)$R2)*(nrow(dep)-1)/
                  (nrow(dep)-(as.numeric(numvars))-1))
              aicg[i]=cens_mod$AIC
            }
          }
          else {
            for(k in 2:numvars) {
              vars = top_n_list[[which(flow==flows)]][[k]]
              for(i in 1:n) {
                cens_mod = censReg(as.lcens(dep[,js[j]], cens_level) ~ .,
                  as.data.frame(expl[,colnames(expl)%in%vars[
                    which(vars[,i+1]=="*"),1]==TRUE]),
                  weights = comp.wys/mean(comp.wys))
                r2[i,k]=summary(cens_mod)$R2
                sig[i,k]=sqrt(sum(residuals(cens_mod)^2)/
                    cens_mod$survreg$df.residual)
                r2adj[i,k]=1-((1-summary(cens_mod)$R2)*(nrow(dep)-1)/
                    (nrow(dep)-(as.numeric(numvars))-1))
                aic[i,k]=cens_mod$AIC
              }
            }
          }
        }
        #combine 1 variable GOF stats with others
        r2[,1]=r2g;aic[,1]=aicg;sig[,1]=sigg;r2adj[,1]=r2adjg
        #plot GOF stats
        if (with1.flag) {
          if(g==gof[1]) {
            plot.across.vars1(names=names, n=n, value=r2,type="R^2",js=js,
              numvars=numvars, j=j)
          } else if(g==gof[2]) {
            plot.across.vars1(names=names, n=n, value=r2adj,type="Adj R^2",
              js=js, numvars=numvars, j=j)
          } else if(g==gof[3]) {
            plot.across.vars1(names=names, n=n, value=sig,type="Sigma",js=js,
              numvars=numvars, j=j)
          } else if(g==gof[4]) {
            plot.across.vars1(names=names, n=n, value=aic,type="AIC",js=js,
              numvars=numvars, j=j)
          }
        } else {
          if(g==gof[1]) {
            plot.across.vars(names=names, n=n, value=r2,type="R^2",js=js,
              numvars=numvars, j=j)
          } else if(g==gof[2]) {
            plot.across.vars(names=names, n=n, value=r2adj,type="Adj R^2",
              js=js, numvars=numvars, j=j)
          } else if(g==gof[3]) {
            plot.across.vars(names=names, n=n, value=sig,type="Sigma",js=js,
              numvars=numvars, j=j)
          } else if(g==gof[4]) {
            plot.across.vars(names=names, n=n, value=aic,type="AIC",js=js,
              numvars=numvars, j=j)
          }
        }
      }
    }
    dev.off()
  }
}

#the only differnece between these functions is the plot.across.vars1 includes the 1 variable models
plot.across.vars <- function(value,type,js,numvars, n, j, names) {
  #legend labels
  name = sapply(c(2:numvars), function(x) paste(x, "vars", sep=''))
  color = c("black","red","orange","green","purple","blue")
  name = name[1:(numvars-1)]; color = color[1:(numvars-1)]
  #plot GOF stats
  plot(1:n,value[,2],type="S", ylim=range(value[,-1])+c(-sd(value[,-1]),
    sd(value[,-1])), ylab=type,xlab="Model Number",
    main=names[(js[j])],col=color[1])
  points(1:n,value[,2],pch=16,cex=.5,col=color[1])
  #add other variables GOF stats to plot
  for (k in 3:numvars) {
    points(1:n,value[,k],type="S",col=color[k-1])
    points(1:n,value[,k],pch=16,cex=.5,col=color[k-1])
  }
  legend("topright", lty=rep(1,numvars-1), col=color, legend=name,
    cex=0.6, lty=rep(1, length(js)), lwd=rep(2, length(js)))
  #leg("topleft",name,numvars)
}

plot.across.vars1 <- function(value,type,js,numvars, n, j, names) {
  name1 = c(1:numvars)
  color = c("black","red","orange","green","purple","blue")
  color = color[1:numvars]
  plot(1:n,value[,1],type="S", ylim=range(value)+.5*c(-sd(value),sd(value)),
    ylab=type,xlab="Model Number",main=names[(js[j])],col=color[1])
  points(1:n,value[,1],pch=16,cex=.5,col=color[1])
  for (k in 2:numvars) {
    points(1:n,value[,k],type="S",col=color[k])
    points(1:n,value[,k],pch=16,cex=.5,col=color[k])
  }
  legend("topright", lty=rep(1,numvars), col=color, legend=name1,
    horiz=T, cex=0.6, title="Number of Variables")
}
