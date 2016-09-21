#getBasinChar----
#' Collect and process basin characteristics for required basins.
#'
#' @description
#' The function \code{getBasinChar} collects all of the basin characteristics
#' in GAGES-II for a specific \code{listofgages}.
#'
#' @param listofgages Dataframe of gage numbers (as a character vector)
#' @param basinChars (optional) A character vector that specifies certain
#' variables to retrieve.  The default is to retrieve all.
#' @param destination (optional) A character string specfying a file name,
#' without an extension, to which to write intermediate results.  The default
#' is not to save results.
#' @param BC.code6.remflg (optional)  A logical flag that indicates whether or
#' not variables derived from a non-calibrated hydrologic model will be removed.
#'   The default, \code{FALSE} requests that these be retained.
#' @param max.BC.oneval.frac (optional)  A constant between 0 and 1, such that
#' if at least that fraction of the BC values has a single value (usually zero)
#' then the BC will be removed.  The default is \code{0.5}.
#' @param debug.flg (optional) A logical indicating if verbose outputs are
#' required from variable transformation.  The default is \code{FALSE}.
#' @param keepWB (optional) A logical passed to
#' \code{\link{clean_BCs}} that indicates whether or not
#' Water Balance model variables are retained. The default is \code{FALSE}.
#'
#' @details
#' Lorem ipsum...
#'
#' @return The function returns a list as output.  This list contains:
#' \item{DescBCs}{A data frame containing the station ID (\code{STAID}),
#' drainage area (\code{DRAIN_SQKM}), and latitude and longitude in Albers
#' projection (\code{LAT_GAGE_UTM} and \code{LNG_GAGE_UTM})
#' and decimal degrees (\code{LAT_GAGE} and \code{LNG_GAGE}) for the gage outlet
#' and basin centroid (\code{LAT_CENT} and \code{LoNG_CENT})
#' for the \code{listofgages}.}
#' \item{AllBCs}{A data frame containing all the variables in the database.}
#' \item{cleanBCs}{A data frame of basin characteristics screened to remove
#' categorical data or variables with a large proportion of
#' identical variables.}
#' \item{winnowedBCs}{A data frame of basin characteristics with correlations
#' removed and climates aggregated to seasonal signals.}
#' \item{transf}{\itemize{A list of the variable transformations:
#'   \item{transfBCs}{A data frame of basin characteristics with
#'   transformations for linearity applied.}
#'   \item{transf.info}{A data frame describing all transformations applied.}
#'   \item{transfBC.VIFs}{A vector of the variance inflation factors (VIFs)
#'   of each bsin characteristic.}
#'   \item{transfBC.corr}{A square matrix of the correlations between basin
#'   characteristics.}}
#' }
#'@export
getBasinChar <- function(listofgages,basinChars=NULL,destination="",
  BC.code6.remflg=F,max.BC.oneval.frac=0.5,debug.flg=F, keepWB=F) {
  # Function orginially designed by Stacey A. Archfield, 02 June 2015.
  # Modified by William Farmer, 03 June 2015.
  # Fused with code by Tom Over and Mike Olson, 30 June 2015.
  #   Includes code to clean, winnow and transform basin characteristic variables.
  # Modified by TMO and Amy Russell, 9/2016, to pass parameter keepWB to
  #   clean_BCs function

  # Cut down to only sites requested.
  # NOTE: BasCharRaw is embedded in sysdata.rda
  ndx <- match(as.double(as.character(listofgages$gages)),BasCharRaw$STAID)
  BasChar <- BasCharRaw[ndx,]

  # Save all variables
  if (destination!="") {
    silent <- TRUE
    destination <- paste0(destination,".AllVar.RData")
    save(list=c("BasChar"),file=destination,compress="bzip2")
  }
  row.names(BasChar) <- NULL

  # Reserve standard variables:
  # Added lat/lon of basin centroid and outlet to list of standard variables
  # TMO, 7/2016
  DescChars <- BasChar[,which(is.element(names(BasChar),
                                         c("STAID","DRAIN_SQKM",
                                           "LAT_GAGE_UTM","LNG_GAGE_UTM",
                                           "LAT_GAGE", "LNG_GAGE",
                                           "LAT_CENT", "LONG_CENT")))]

  row.names(BasChar) <- BasChar$STAID

  BasChar <- BasChar[,-which(is.element(names(BasChar),
    c("LAT_GAGE_UTM","LNG_GAGE_UTM","STAID")))]

  ## CODE DEVELOPED BY TOM AND MIKE

  #Create suffix to be used (when silent = F) as part of output files
  if (BC.code6.remflg) code6_str = "rem_code6" else code6_str = "keep_code6"
  BC_suffix = paste(code6_str,".max_oneval_frac",max.BC.oneval.frac,sep="")

  cleanBCs = clean_BCs(BasChar, BC.code6.remflg, max.BC.oneval.frac,
    BC_suffix, destination, keepWB=keepWB)

  winnowedBCs = winnow_BCs(cleanBCs, BC.code6.remflg)

  transfBCs.list = transform_BCs(winnowedBCs, debug.flg, BC_suffix, destination)
  #Elements of transfBCs.list:
  transfBCs = transfBCs.list[[1]] #Data frame of transformed BCs
  transf.info = transfBCs.list[[2]]
  #Data frame containing information describing transformations,
  #three rows and one column for each BC.
  transfBC.VIFs = transfBCs.list[[3]]
  #Data frame containing VIFs for each BC (one row),
  #computed by predicting each BC from all other BCs
  #(to help determine redundant BCs
  #if too many models kicked out by subsequent leaps processing for high VIF).
  transfBC.corr = transfBCs.list[[4]]
  #Correlation matrix of BCs, also to help determine
  #redundant BCs.

  ## END CODE DEVELOPED BY TOM AND MIKE

  # Cut down to only variables requested.
  if (!is.null(basinChars)) {
    ndx <- which(is.element(names(cleanBCs),basinChars))
    cleanBCs <- cleanBCs[,ndx]
    ndx <- which(is.element(names(winnowedBCs),basinChars))
    winnowedBCs <- winnowedBCs[,ndx]
    ndx <- which(is.element(names(transfBCs),basinChars))
    transfBCs <- transfBCs[,ndx]
    ndx <- ndx-1 # Removes STAID
    transf.info <- transf.info[,ndx]
    transfBC.VIFs <- transfBC.VIFs[ndx]
    transfBC.corr <- transfBC.corr[ndx,ndx]
  }

  result <- list(DescBCs=DescChars,
    AllBCs=BasChar,
    cleanBCs=cleanBCs,
    winnowedBCs=winnowedBCs,
    transf=list(transfBCs=transfBCs,transf.info=transf.info,
      transfBC.VIFs=transfBC.VIFs,transfBC.corr=transfBC.corr))

  return(result)
}
