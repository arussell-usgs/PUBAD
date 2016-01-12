# winnow_BCs ----
#' Winnow Basin Characteristics
#'
#' @description
#' The function \code{winnow_BCs} removes correlated basin characteristics
#' and aggregate monthly variables to seasonal variables.
#'
#' @param BCs A data frame of basin characteristics from
#' \link{\code{clean_BCs.R}}.
#' @param BC.code6.remflg (optional) A logical value indicating if the code-6
#' variables, those derived from an uncalibrated hydrologic model,
#' should be removed.  The defautl is \code{FALSE}.
#'
#' @details
#' Reads in the "clean" basin characteristics (BCs) produced by
#' \link{\code{clean_BCs.R}} and then "winnows" them by:
#' \item{1.}{"pre-removing" various highly correlated or unlikely to be useful
#' BCs, including "SITE", "STD", "CDL", forest NLCD variables other than total
#' forest (FORESTNLCD06).}
#' \item{2.}{Add a few BCs are combinations of existing BCs, including:
#' \item{(a)}{July - Jan monthly temperatures (then all other monthly
#' temperatures are removed; annual average is also retained.)}
#' \item{(b)}{Phase of annual cycle of monthly water balance WB (WB5100...)
#' variables, then delete monthly WB variables but retain annual WB variable
#' WB5100_ANN_MM.}
#' \item{(c)}{Combining monthly PPT variables into quarterly
#' (DJF, MAM, JJA, SON)}
#' \item{(d)}{Combining woody and emergent NLCD fractions to a single total
#' wetland fraction; also creating a NLCD total wetland plus NLCD water
#' variable}
#' \item{(e)}{Adding ELEV_RANGE_M_BASIN variable as difference of
#' ELEV_MAX_M_BASIN and ELEV_MIN_M_BASIN}}
#'
#' @return A matrix of basin characteristics, with some removed.
#'@export
winnow_BCs <- function(BCs, BC.code6.remflg=F) {
  # Function orginially designed by Thomas M. Over and Mike Olsen, 30 June 2015.
  # Modified by William Farmer, 30 June 2015.

  #Define stations numbers and how many stations based on BC input
  statn.nums = BCs$STAID; nstatns = nrow(BCs)

  #Remove lat-lon at gage, other SITE variables, STD and CDL variables
  BCs$LAT_GAGE = BCs$LNG_GAGE = NULL
  #"SITE" variables are: PPTAVG_SITE, T_AVG_SITE, T_MAX_SITE, T_MIN_SITE,
  #RH_SITE, FST32SITE, LST32SITE, WD_SITE, WDMAX_SITE, WDMIN_SITE, ELEV_SITE_M;
  #there are "BASIN" versions of all which seem physically more meaningful
  BCs = BCs[grep("SITE",names(BCs),invert=T)]
  #"STD" variables are: T_MAXSTD_BASIN, T_MINSTD_BASIN, ELEV_STD_M_BASIN;
  #it is not clear what is physical relevance of these
  BCs = BCs[grep("STD",names(BCs),invert=T)]
  BCs = BCs[grep("CDL",names(BCs),invert=T)]

  #Remove add'l BCs because of generally high correlations with other
  #retained BCs

  #For Upper MO, T_MAX_BASIN and T_MIN_BASIN are highly correlated with
  #T_AVG_BASIN, and that seeems not unlikely to be a general case,
  #so remove them
  BCs$T_MAX_BASIN = BCs$T_MIN_BASIN = NULL

  #For Upper MO, FORESTNLCD06 and EVERGRNLCD06 are very highly correlated.
  #Could imagine that all subtypes of forest (DECID, EVERGR, and MIXEDFOR)
  #should be deleted and FORESTNLCD06 retained, but could also imagine a
  #scenario where say DECID versus EVERGR is indicative of some physically
  #meaningful difference.
  #Nevertheless, for sake of automation, remove all but FOREST.
  BCs$DECIDNLCD06 = BCs$EVERGRNLCD06 = BCs$MIXEDFORNLCD06 = NULL

  #Create combined WOODYWETNLCD06 and EMERGWETNLCD06 variable and
  #remove them separate version;
  if (sum(c('WOODYWETNLCD06','EMERGWETNLCD06') %in% names(BCs)) == 2) {
    # WHF added if statement in case variables are removed.
    BCs$TOTWETNLCD06 = BCs$WOODYWETNLCD06 + BCs$EMERGWETNLCD06
  }
  BCs$WOODYWETNLCD06 = BCs$EMERGWETNLCD06 = NULL

  #Also create a combined wetland and water variable
  #(at risk of raising correlations)
  if (sum(c('TOTWETNLCD06','WATERNLCD06') %in% names(BCs)) == 2) {
    # WHF added if statement in case variables are removed.
    BCs$WATERWETNLCD06 = BCs$TOTWETNLCD06 + BCs$WATERNLCD06
  }
  BCs$TOTWETNLCD06 = BCs$WATERNLCD06 = NULL

  #In Upper MO, soil size variables NO10AVE and NO4AVE are very highly
  #correlated, and this seems likely to be generally true, since they are
  #very similar measures.
  #Since the other soil size variable NO200AVE is more different from NO4AVE,
  #retain that one.
  BCs$NO10AVE = NULL

  #RRMEAN and RRMEDIAN are highly correlated in Upper MO and this seems likely
  #to be generally true because of the similar way they are calculated;
  #removing RRMEAN
  BCs$RRMEAN = NULL

  #ELEV_MEAN_M_BASIN and ELEV_MEDIAN_M_BASIN are highly correlated in Upper MO,
  #and this seems likely to be generally true because of the similar way
  #they are calculated; removing ELEV_MEAN_M_BASIN
  BCs$ELEV_MEAN_M_BASIN = NULL

  # Remove empty BCs
  if (length(which(colSums(is.na(BCs))==nrow(BCs)))>0){
    BCs <- BCs[-which(colSums(is.na(BCs))==nrow(BCs))]
  }

  #Combine monthly ppt to quarterly (DJF, MAM, JJA, SON)
  DJF.ppt.names = paste(c("DEC","JAN","FEB"),"_PPT7100_CM",sep="")
  if (sum(DJF.ppt.names %in% names(BCs)) == 3) {
    #Add sum of monthly ppt to dataframe
    BCs$DJF_PPT7100_CM = as.numeric(unlist(BCs[DJF.ppt.names[1]] +
        BCs[DJF.ppt.names[2]] + BCs[DJF.ppt.names[3]]))
    #Remove separate monthly ppt from dataframe
    for (i in 1:3) BCs[DJF.ppt.names[i]] = NULL
  }
  MAM.ppt.names = paste(c("MAR","APR","MAY"),"_PPT7100_CM",sep="")
  if (sum(MAM.ppt.names %in% names(BCs)) == 3) {
    #Add sum of monthly ppt to dataframe
    BCs$MAM_PPT7100_CM = as.numeric(unlist(BCs[MAM.ppt.names[1]] +
        BCs[MAM.ppt.names[2]] + BCs[MAM.ppt.names[3]]))
    #Remove separate monthly ppt from dataframe
    for (i in 1:3) BCs[MAM.ppt.names[i]] = NULL
  }
  JJA.ppt.names = paste(c("JUN","JUL","AUG"),"_PPT7100_CM",sep="")
  if (sum(JJA.ppt.names %in% names(BCs)) == 3) {
    #Add sum of monthly ppt to dataframe
    BCs$JJA_PPT7100_CM = as.numeric(unlist(BCs[JJA.ppt.names[1]] +
        BCs[JJA.ppt.names[2]] + BCs[JJA.ppt.names[3]]))
    #Remove separate monthly ppt from dataframe
    for (i in 1:3) BCs[JJA.ppt.names[i]] = NULL
  }
  SON.ppt.names = paste(c("SEP","OCT","NOV"),"_PPT7100_CM",sep="")
  if (sum(SON.ppt.names %in% names(BCs)) == 3) {
    #Add sum of monthly ppt to dataframe
    BCs$SON_PPT7100_CM = as.numeric(unlist(BCs[SON.ppt.names[1]] +
        BCs[SON.ppt.names[2]] + BCs[SON.ppt.names[3]]))
    #Remove separate monthly ppt from dataframe
    for (i in 1:3) BCs[SON.ppt.names[i]] = NULL
  }

  #If present, compute amplitude and phase of monthly WB (Water balance model)
  #discharges, following D.S. Wilks (2006), Stat. Methods in the Atmos.
  #Sciences, 2nd ed., Sec. 8.4.3
  #However, as it was found that mon_WB.amp was very strongly correlated with
  #mean annual WB discharge WB WB5100_ANN_MM (for upper MO basins), only the
  #phase (and mean annual) are retained; all the monthly values are removed.
  if (!BC.code6.remflg) {
    month.abbrevs = c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP",
      "OCT","NOV","DEC")
    mon_WB.names = paste("WB5100_",month.abbrevs,"_MM",sep="")
    if (sum(mon_WB.names %in% names(BCs)) == 12) {
      # WHF added if statement in case variables are removed.
      x.cos = cos(2*pi*(1:12)/12); x.sin = sin(2*pi*(1:12)/12)
      mon_WB.phase = mon_WB.amp = numeric(nstatns)
      for (i in 1:nstatns) {
        lm.mod = lm(as.numeric(BCs[mon_WB.names][i,]) ~ x.cos + x.sin)
        cos.coef = coef(lm.mod)[2]; sin.coef = coef(lm.mod)[3]
        mon_WB.amp[i] = sqrt(cos.coef^2 + sin.coef^2)
        if (cos.coef > 0) mon_WB.phase[i] = atan(sin.coef/cos.coef) else
          if (cos.coef < 0) mon_WB.phase[i] = atan(sin.coef/cos.coef) + pi else
            mon_WB.phase[i] = 0.
      }
      BCs$WB5100_MON_PHASE = mon_WB.phase
      for (i in 1:12) BCs[mon_WB.names[i]] = NULL
    }
  }

  #Was going to also compute amplitude and phase of monthly temperatures,
  #but it was found that the phase was almost constant among basins
  #(for Upper MO stations), and the amplitude was very strongly correlated
  #with Jul - Jan temperature, so may as well just use the July - Jan
  #temperature and remove all the others.
  #(Recall mean annual temperature is still present.)
  month.abbrevs = c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP",
    "OCT","NOV","DEC")
  mon_tmp.names = paste(month.abbrevs,"_TMP7100_DEGC",sep="")
  if (sum(c('JUL_TMP7100_DEGC','JAN_TMP7100_DEGC') %in% names(BCs)) == 2) {
    # WHF added if statement in case variables are removed.
    BCs$JUL_JAN_TMP_DIFF_DEGC = BCs$JUL_TMP7100_DEGC - BCs$JAN_TMP7100_DEGC
  }
  for (i in 1:12) BCs[mon_tmp.names[i]] = NULL

  #Create ELEV_RANGE variable
  if (sum(c('ELEV_MAX_M_BASIN','ELEV_MIN_M_BASIN') %in% names(BCs)) == 2) {
    # WHF added if statement in case variables are removed.
    BCs$ELEV_RANGE_M_BASIN = BCs$ELEV_MAX_M_BASIN - BCs$ELEV_MIN_M_BASIN
  }
  BCs$ELEV_MAX_M_BASIN = BCs$ELEV_MIN_M_BASIN = NULL

  return(BCs)
}
