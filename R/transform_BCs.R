# transform_BCs ----
#' Develop transformations of basin characteristics.
#'
#' @description
#' The function \code{transform_BCs} develops transformations of basin
#' characteristics to induce linearity.
#'
#' @param BCs A data frame of basin characteristics from
#' \code{\link{winnow_BCs}}.
#' @param debug.flg (optional) A logical indicating if verbose outputs are
#' required from variable transformation.  The default is \code{FALSE}.
#' @param BC_sfx (optional) A character string appending certain condition
#' onto saved file names.
#' @param destination (optional) A character string specfying a file name,
#' without an extension, to which to write intermediate results.  The default
#' is not to save results.
#'
#' @details
#' Receives the "winnowed" basin characteristics (BCs) produced by
#' \code{\link{winndow_BCs}} and then determines need for transformation of
#' BCs by computing some summary statistics and then skewness for original and
#' centered data (if needed) and different transformations (currently powers
#' of 1/10, 1/5, 1/3, 1/2, 2/3, 3/2, 2, 3, 5 and 10 plus log10 and 10^)
#' following addition of constants as needed to remove negative and
#' non-positive values as needed, depending on the transformation.
#' Then skewness information is used to choose the preferred transformation
#' as the one which gives the skewness nearest to zero.
#' Following selection of best transformation, post-transformation centering is
#' applied as needed.
#'
#' @return A list of four elements:
#' \item{finalBCs}{A data frame of the final, transformed basin
#' characteristics.}
#' \item{transfBC.info}{A data fram of information to recreate variable
#' transformations.  This includes the constant added before transformation,
#' the transformation used and the constant added after transformation.}
#' \item{finalBC.VIFs}{The variance inflation factors between basin
#' characteristics.}
#' \item{finalBC.corr}{The cross-correlation matrix of final, transformed
#' basin characteristics.}
#'@export
#'@import moments
transform_BCs <- function(BCs, debug.flg=F, BC_sfx="", destination="") {
  # Function orginially designed by Thomas M. Over and Mike Olsen, 30 June 2015.
  # Modified by William Farmer, 30 June 2015.
  # Revised by TMO 02 July 2015, Implemented by WHF 06 July 2015

  # @importFrom moments skewness

  #Compute number of BCs in input data frame
  class_col = which(names(BCs)=="CLASS")
  col_order = c(class_col,(1:ncol(BCs))[-class_col])
  BCs = BCs[,col_order]
  nBCs = ncol(BCs)-1
  #(continue to assume first row is CLASS and station ID is row name)

  silent <- TRUE
  if (destination!="") {
    silent <- FALSE
    out_folder <- destination
  }

  #Compute summary statistics
  sum.stat.names = c("min","1st_quart","median","mean","3rd_quart","max",
    "stdev","CV","skew")
  origBC.sum_stats = data.frame(matrix(nrow=length(sum.stat.names),ncol=nBCs,
    dimnames=list(sum.stat.names,names(BCs)[-1])))
  CVs = numeric(nBCs)
  for (i in 1:nBCs) {
    origBC.sum_stats[1:6,i] = as.numeric(summary(BCs[,i+1]))
    origBC.sum_stats[7,i] = sd(BCs[,i+1])
    origBC.sum_stats[9,i] = skewness(BCs[,i+1])
    if (mean(BCs[,i+1])==0) {
      CVs[i] = NA
    } else {
      CVs[i] = sd(BCs[,i+1])/mean(BCs[,i+1])
    }
    origBC.sum_stats[8,i] = CVs[i]
  }
  if (!silent) write.csv(origBC.sum_stats,paste(out_folder,
    ".origBC.sum_stats.",BC_sfx,".csv",sep=""))

  #Compute skewness of each BC as given and after various transformation
  #create BC arrays for transformed versions
  BCs.exp10 = BCs.power5 = BCs.power3 = BCs.power10 = BCs.power2 =
    BCs.log10 = BCs
  BCs.1over10 = BCs.1over5 = BCs.1over3 = BCs.1over2 = BCs.2over3 =
    BCs.3over2 = BCs.center = BCs;
  tran_names = c("orig","center","3over2","2over3","1over2","1over3","1over5"
    ,"1over10","log10","power2","power10","power3","power5","exp10")

  BC_skews = data.frame(matrix(nrow=length(tran_names),ncol=nBCs,
    dimnames=list(tran_names,names(BCs)[-1])))
  BC_added = BC_skews; BC_added[,] = 0 #initialize BC_added arrays to zeroes.
  center.minCV = 0.10 #Center a BC whose abs(CV) (sd/mean) < center.minCV

  #Iterate over each column (BC)
  for (i in 2:ncol(BCs)) {
    #Skewness of original values
    type.cnt = 1
    BC_skews[type.cnt,i-1] = skewness(BCs[,i])

    #Centered values
    type.cnt = 2
    if (abs(CVs[i-1]) < center.minCV) {
      temp.BCs = BCs[,i] - mean(BCs[,i])
      BC_added[type.cnt,i-1] = -mean(BCs[,i])
    } else temp.BCs = BCs[,i]
    BC_skews[type.cnt,i-1] = skewness(temp.BCs); BCs.center[,i] = temp.BCs

    #NB: From here on out, all transformed values will
    #be based on centered values;
    #additional shifts are applied as needed to obtain
    #non-negative or positive values.

    #Sqrt-transform: To take sqrt transform,
    #need to make all values non-negative
    #For additional "resolution", do also 3/2, 2/3, 1/3, 1/5, 1/10 powers
    #For all these need to make all values non-negative.
    if (min(BCs.center[,i])<0) {
      temp.BCs = BCs.center[,i] - min(BCs.center[,i])
      BC_added[3:8,i-1] = -min(BCs.center[,i])
    } else temp.BCs = BCs.center[,i]
    type.cnt=3; BC_skews[type.cnt,i-1] = skewness(temp.BCs^(3/2))
    BCs.3over2[,i] = temp.BCs^(3/2)
    type.cnt=4; BC_skews[type.cnt,i-1] = skewness(temp.BCs^(2/3))
    BCs.2over3[,i] = temp.BCs^(2/3)
    type.cnt=5; BC_skews[type.cnt,i-1] = skewness(temp.BCs^(1/2))
    BCs.1over2[,i] = temp.BCs^(1/2)
    type.cnt=6; BC_skews[type.cnt,i-1] = skewness(temp.BCs^(1/3))
    BCs.1over3[,i] = temp.BCs^(1/3)
    type.cnt=7; BC_skews[type.cnt,i-1] = skewness(temp.BCs^(1/5))
    BCs.1over5[,i] = temp.BCs^(1/5)
    type.cnt=8; BC_skews[type.cnt,i-1] = skewness(temp.BCs^(1/10))
    BCs.1over10[,i] = temp.BCs^(1/10)

    #log10-transform: To take log10 transform, need to make all values positive
    type.cnt = 9
    if (min(BCs.center[,i])<0) {
      temp.BCs = BCs.center[,i]
      BC_added[type.cnt,i-1] = -min(BCs.center[,i])
      temp.BCs = temp.BCs + BC_added[type.cnt,i-1]
      temp.BCs = temp.BCs + 0.1*min(temp.BCs[temp.BCs>0])
      BC_added[type.cnt,i-1] =
        BC_added[type.cnt,i-1] + 0.1*min(temp.BCs[temp.BCs>0])
    } else if (min(BCs[,i])==0) {
      temp.BCs = BCs.center[,i]
      BC_added[type.cnt,i-1] = 0.1*min(temp.BCs[temp.BCs>0])
      temp.BCs = temp.BCs + BC_added[type.cnt,i-1]
    } else temp.BCs = BCs.center[,i]
    BC_skews[type.cnt,i-1] = skewness(log10(temp.BCs))
    BCs.log10[,i] = log10(temp.BCs)

    #powers 2 and 10: To preserve ordering,
    #need to use negative of even powers of negatives
    temp.BCs = BCs.center[,i]
    neg.subs = which(temp.BCs<0)
    if (length(neg.subs)>0) {
      BCs.power2[neg.subs,i] = -(temp.BCs[neg.subs])^2
      BCs.power10[neg.subs,i] = -(temp.BCs[neg.subs])^10
    }
    noneg.subs = which(temp.BCs>=0)
    if (length(noneg.subs)>0) {
      BCs.power2[noneg.subs,i] = temp.BCs[noneg.subs]^2
      BCs.power10[noneg.subs,i] = temp.BCs[noneg.subs]^10
    }
    type.cnt=10; BC_skews[type.cnt,i-1] = skewness(BCs.power2[,i])
    type.cnt=11; BC_skews[type.cnt,i-1] = skewness(BCs.power10[,i])

    #powers 3, 5, and exp10 values (meaning values used as exponent of 10):
    #no special handling needed
    temp.BCs = BCs.center[,i]
    BCs.power3[,i] = temp.BCs^3
    BCs.power5[,i] = temp.BCs^5
    BCs.exp10[,i] = 10^temp.BCs
    type.cnt=12; BC_skews[type.cnt,i-1] = skewness(BCs.power3[,i])
    type.cnt=13; BC_skews[type.cnt,i-1] = skewness(BCs.power5[,i])
    type.cnt=14; BC_skews[type.cnt,i-1] = skewness(BCs.exp10[,i])

    #Add BC.center BC_added values to subsequent types
    for (type.cnt in 3:14) {
      BC_added[type.cnt,i-1] = BC_added[type.cnt,i-1] + BC_added[2,i-1]
    }
  }

  #Writing these out is only done for initial debugging
  if (debug.flg) {
    write.csv(BC_skews,paste(out_folder,".BC_skews.",BC_sfx,".csv",sep=""))
    write.csv(BC_added,paste(out_folder,".BC_added.",BC_sfx,".csv",sep=""))
    write.csv(BCs.center,paste(out_folder,
      ".BCs.centered.",BC_sfx,".csv",sep=""))
    write.csv(BCs.3over2,paste(out_folder,".BCs.3over2.",BC_sfx,".csv",sep=""))
    write.csv(BCs.2over3,paste(out_folder,".BCs.2over3.",BC_sfx,".csv",sep=""))
    write.csv(BCs.1over2,paste(out_folder,".BCs.1over2.",BC_sfx,".csv",sep=""))
    write.csv(BCs.1over3,paste(out_folder,".BCs.1over3.",BC_sfx,".csv",sep=""))
    write.csv(BCs.1over5,paste(out_folder,".BCs.1over5.",BC_sfx,".csv",sep=""))
    write.csv(BCs.1over10,paste(out_folder,
      ".BCs.1over10.",BC_sfx,".csv",sep=""))
    write.csv(BCs.log10,paste(out_folder,".BCs.log10.",BC_sfx,".csv",sep=""))
    write.csv(BCs.power2,paste(out_folder,".BCs.power2.",BC_sfx,".csv",sep=""))
    write.csv(BCs.power10,paste(out_folder,
      ".BCs.power10.",BC_sfx,".csv",sep=""))
    write.csv(BCs.power3,paste(out_folder,".BCs.power3.",BC_sfx,".csv",sep=""))
    write.csv(BCs.power5,paste(out_folder,".BCs.power5.",BC_sfx,".csv",sep=""))
    write.csv(BCs.exp10,paste(out_folder,".BCs.exp10.",BC_sfx,".csv",sep=""))
  }

  #Select appropriate transformation (if any) dependent on which gives
  #skewness closest to zero
  #tran_names = c("orig","center","3over2","2over3","1over2","1over3",
  #               "1over5","1over10","log10","power2","power10",
  #               "power3",power5","exp10")
  transfBCs = BCs
  transfBCs.tran = character(nBCs)
  transfBCs.added = numeric(nBCs)
  transfBCs.skew = numeric(nBCs)
  for (i in 2:ncol(BCs)) {
    minrow = which.min(abs(BC_skews[,i-1]))
    #Do not consider "original" version, because if not actually centered,
    #"center" is the same as original; otherwise "center" is preferred.
    if (minrow==1) minrow = 2
    transfBCs.tran[i-1] = tran_names[minrow]
    transfBCs.skew[i-1] = BC_skews[minrow,i-1]
    transfBCs.added[i-1] = BC_added[minrow,i-1]
    if (minrow==1) {} else
      if (minrow==2) transfBCs[,i] = BCs.center[,i] else
        if (minrow==3) transfBCs[,i] = BCs.3over2[,i] else
          if (minrow==4) transfBCs[,i] = BCs.2over3[,i] else
            if (minrow==5) transfBCs[,i] = BCs.1over2[,i] else
              if (minrow==6) transfBCs[,i] = BCs.1over3[,i] else
                if (minrow==7) transfBCs[,i] = BCs.1over5[,i] else
                  if (minrow==8) transfBCs[,i] = BCs.1over10[,i] else
                    if (minrow==9) transfBCs[,i] = BCs.log10[,i] else
                      if (minrow==10) transfBCs[,i] = BCs.power2[,i] else
                        if (minrow==11) transfBCs[,i] = BCs.power10[,i] else
                          if (minrow==12) transfBCs[,i] = BCs.power3[,i] else
                            if (minrow==13) transfBCs[,i] = BCs.power5[,i] else
                              if (minrow==14) transfBCs[,i] = BCs.exp10[,i]
  }

  #write out pre-recentered version for debugging
  if (!silent) write.csv(transfBCs,paste(out_folder,".transfBCs.",
    BC_sfx,".csv",sep=""))

  #Compute the usual summary statistics for pre-centered transformed BCs
  transfBC.sum_stats = origBC.sum_stats
  row.names(transfBC.sum_stats) = paste("transf",
    row.names(origBC.sum_stats),sep=".")
  transfCVs = numeric(nBCs)
  for (i in 1:nBCs) {
    transfBC.sum_stats[1:6,i] = as.numeric(summary(transfBCs[,i+1]))
    transfBC.sum_stats[7,i] = sd(transfBCs[,i+1])
    transfBC.sum_stats[9,i] = skewness(transfBCs[,i+1])
    if (mean(transfBCs[,i+1])==0) transfCVs[i] = NA else
      transfCVs[i] = sd(transfBCs[,i+1])/mean(transfBCs[,i+1])
    transfBC.sum_stats[8,i] = transfCVs[i]
  }

  #Do post-transform centering as needed (if transfCV is < center.minCV),
  #creating final version of transformed BCs
  #Revision 1, 7/2/2015: Also divide each BC vector Xi by max(abs(Xi))
  #to make the magnitude of BC values among BCs comparable.
  finalBCs = transfBCs; post.transfBCs.added = numeric(nBCs)
  names(post.transfBCs.added) = names(finalBCs)[-1]
  finalBCs.divisor = numeric(nBCs) #added 7/2/2015, TMO
  names(finalBCs.divisor) = names(finalBCs)[-1] #added 7/2/2015, TMO
  for (i in 2:ncol(finalBCs)) {
    if (ifelse(is.na(abs(transfCVs[i-1])),999,abs(transfCVs[i-1])) < center.minCV) {
      finalBCs[,i] = transfBCs[,i] - mean(transfBCs[,i])
      post.transfBCs.added[i-1] = -mean(transfBCs[,i])
    }
    finalBCs.divisor[i-1] = max(abs(finalBCs[,i])) #added 7/2/2015, TMO
    finalBCs[,i] = finalBCs[,i]/finalBCs.divisor[i-1] #added 7/2/2015, TMO
  }

  #write out final BCs
  if (!silent) write.csv(finalBCs,paste(out_folder,
    ".finalBCs.",BC_sfx,".csv",sep=""))

  #Create transfBC.info: constant added before transformation,
  #transformation used, and constant added after transformation
  #(added 7/2/2015): and final BC divisor
  transfBC.info = data.frame(rbind(transfBCs.added, transfBCs.tran,
    post.transfBCs.added,finalBCs.divisor)) #added 7/2/2015
  names(transfBC.info) = names(transfBCs[-1])

  #Re-compute summary statistics (now final)
  finalBC.sum_stats = transfBC.sum_stats
  row.names(finalBC.sum_stats) = paste("final",
    row.names(origBC.sum_stats),sep=".")
  finalCVs = numeric(nBCs)
  for (i in 1:nBCs) {
    finalBC.sum_stats[1:6,i] = as.numeric(summary(finalBCs[,i+1]))
    finalBC.sum_stats[7,i] = sd(finalBCs[,i+1])
    if (mean(finalBCs[,i+1])==0) finalCVs[i] = NA else
      finalCVs[i] = sd(finalBCs[,i+1])/mean(finalBCs[,i+1])
    finalBC.sum_stats[8,i] = finalCVs[i]
  }

  #Compute VIF values for each variable as an
  #additional aid to finding redundant BCs
  finalBC.VIFs = matrix(nrow=1,ncol=ncol(finalBCs)-1,
    dimnames=list("finalBC.VIF",names(finalBCs[-1])))
  for (i in 2:ncol(finalBCs)) {
    x.data = finalBCs[,seq(2,ncol(finalBCs))[-(i-1)]]
    lm.mod = lm(finalBCs[,i] ~ ., data=x.data)
    finalBC.VIFs[i-1] = 1 / (1 - summary(lm.mod)$r.squared)
  }

  #Append transfBC.sum_stats, transfBC.info, finalBC.sum_stats, and finalBC.VIFs
  #to original BC.sum_stats and write out
  combBC.sum_stats = rbind(origBC.sum_stats, transfBC.sum_stats, transfBC.info,
    finalBC.sum_stats,finalBC.VIFs)
  if (!silent) write.csv(combBC.sum_stats,paste(out_folder,
    ".combBC.sum_stats.",BC_sfx,".csv",sep=""))

  #Plot original and transformed BCs for checking, e.g.,
  #to make sure transformation is monotonic.
  if (!silent) {
    pdf(paste(out_folder,".final_v_orig_BCs.",
      BC_sfx,".pdf",sep=""),height=11,width=8.5)
    par(mfrow=c(3,2))
    for (i in 2:ncol(finalBCs)) {
      plot(BCs[,i],finalBCs[,i],
        xlab=paste("orig",names(finalBCs)[i]),
        ylab=paste("final",names(finalBCs)[i]))
      #Most transformations create a concave plot so put legend at bottom right
      #(Could figure out which is which but
      #seems unnecessary for initial version.)
      legend("bottomright",c(paste("pre-transf const added:",
        format(transfBCs.added[i-1])),
        paste("transform:", transfBCs.tran[i-1]),
        paste("post-transf const added:", format(post.transfBCs.added[i-1]))))
    }
    dev.off()

    pdf(paste(out_folder,".final_and_orig_BCs_hists.",BC_sfx,".pdf",sep=""),
      height=11,width=8.5)
    par(mfrow=c(3,2))
    for (i in 2:ncol(finalBCs)) {
      hist(BCs[,i],breaks="FD",xlab="original",main=names(finalBCs)[i])
      hist(finalBCs[,i],breaks="FD",xlab="final",main=names(finalBCs)[i])
    }
    dev.off()
  }

  #Finally check for highly correlated finalBCs
  #Write out correlation matrix
  finalBC.corr = cor(finalBCs[,-1])
  if (!silent) {
    write.csv(finalBC.corr,paste(out_folder,".finalBC.corr.",
      BC_sfx,".csv",sep=""))
    #Write out, for each variable, which other variables have
    #corr > max.corr or < -max.corr
    max.corr = 0.95
    sink(paste(out_folder,".transfBC.big_corr.",BC_sfx,".txt",sep=""))
    for (i in 2:ncol(finalBCs)) {
      big_corr.subs = which(finalBC.corr[i-1,] > max.corr |
          finalBC.corr[i-1,] < (-max.corr))
      big_corr.names = names(finalBCs)[big_corr.subs+1]
      big_corr.vals = finalBC.corr[i-1,big_corr.subs]
      #sort from largest to smallest so, for one thing,
      #the unit value on the diagonal comes first
      sort.out = sort(big_corr.vals,decreasing=T,index.return=T)
      big_corr.vals.sort = sort.out$x
      big_corr.names.sort = big_corr.names[sort.out$ix]
      for (j in 1:length(big_corr.subs)) {
        cat(big_corr.names.sort[j]," ",big_corr.vals.sort[j])
        if (j<length(big_corr.subs)) cat("; ")
      }
      cat("\n")
    }
    sink()
  }

  return(list(finalBCs=finalBCs, transfBC.info=transfBC.info,
    finalBC.VIFs=finalBC.VIFs, finalBC.corr=finalBC.corr))
}
