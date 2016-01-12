#' Cross-compare models.
#'
#' @description
#' A sub-function used by \link{\code{compile.vars}}.
#'
#' @param out Unclear.
#' @param y Unclear.
#' @param ranks Unclear.
#' @param vars Unclear.
#' @param names Unclear.
#'
#' @details
#' Lorem ipsum...
#'
#' @return Unclear.
#'@export
compare<-function(out, y, ranks, vars, names){
  # Function orginially designed by Thomas M. Over and Mike Olsen, 03 July 2015.
  # Modified by William Farmer, 06 July 2015.
  leny = length(y)
  k = 1
  while(k<leny){
    if(is.null(dim(y[[k]])[1])) {
      l = 1
    } else {
      l = length(y[[k]][1,])
    }
    for(i in 1:l){
      rank = ranks[[k]][i]
      whi = names[k]
      num = 1
      #only look at quantiles after this one
      j = (k+1)
      if(j == (length(y)+1)){
        next;
      }
      if(is.null(dim(y[[k]])[1])){
        tempy = y[[k]]
      } else {
        tempy = y[[k]][,i]
      }
      #look for this subset of EVs in each other quantile
      while(j<=length(y)){
        #if all subsets have been removed from quantile, continue
        if(length(y[[j]])==0) {
          j <- j + 1 # WHF, 07/06/2015: Added to avoid infinite loop. ?
          next
        }
        #if this subset appears in another quantile: increase the rank,
        #remove from that quantile (so it does not appear twice)
        if(
          length(which(apply(t(y[[j]]),1, function(x) all(x %in% tempy))))>0 ||
            (vars == 1 && length(which(y[[j]] %in% tempy))>0)
        ){
          if(vars == 1) {
            rank = rank + ranks[[j]][which(y[[j]] %in% tempy)]
          } else {
            rank = rank+ranks[[j]][c(which(apply(t(y[[j]]),1,
              function(x) all(x %in% tempy))))]
          }
          #add this quantile to list of which quantiles this subeset appears in
          whi = paste(whi, ",", names[j], sep = "")
          num = num+1
          if(vars == 1){
            ranks[[j]] = ranks[[j]][-c(which(y[[j]] %in% tempy))]
            y[[j]] = y[[j]][-c(which(y[[j]] %in% tempy))]
          }
          if((is.null(dim(y[[j]]))&& vars!=1) || (length(y[[j]]) == 0)){
            ranks[[j]] = NULL
            names = names[-j]
            y[[j]] = NULL
            j = j-1
          } else if(vars!=1){
            ranks[[j]] = ranks[[j]][-c(which(apply(t(y[[j]]),1,
              function(x) all(x %in% tempy))))]
            y[[j]] = y[[j]][,-c(which(apply(t(y[[j]]),1,
              function(x) all(x %in% tempy))))]
          }
        }
        j = j+1
      }
      t = array(rep(" ", length(out[,1])))
      out = cbind(out, t)
      out[which(out[,1] %in% tempy),length(out[1,])] = c("*")
      out[(length(out[,1])-2),length(out[1,])] = as.numeric(rank)
      out[(length(out[,1])-1),length(out[1,])] = as.numeric(num)
      out[(length(out[,1])),length(out[1,])] = whi
    }
    k = k+1
    if(length(y)<leny){
      leny = length(y)
      k = k-1
    }
    return(out)
  }
}
