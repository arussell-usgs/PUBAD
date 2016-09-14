# cleanBCs ----
#' Clean out raw basin characteristics.
#'
#' @description
#' The function \code{clean_BCs} cleans the raw basin characteristics by
#' removing unnecessary or not-useful variables.
#'
#' @param BCs A data frame of basin characteristics from
#' \code{\link{winnow_BCs}}.
#' @param BC.code6.remflg (optional)  A logical flag that indicates whether or
#' not variables derived from a non-calibrated hydrologic model will be removed.
#'   The default, \code{FALSE} requests that these be retained.
#' @param max.BC.oneval.frac (optional)  A constant between 0 and 1, such that
#' if at least that fraction of the BC values has a single value (usually zero)
#' then the BC will be removed.  The default is \code{0.5}.
#' @param BC_suffix (optional) A character string appending certain condition
#' onto saved file names.
#' @param destination (optional) A character string specfying a file name,
#' without an extension, to which to write intermediate results.  The default
#' is not to save results.
#' @param keepWB (optional) A logical flag that indicates whether or not
#'   Water Balance model variables, named "WB5100...", are kept even if a large
#'   fraction of the values are the same in order to later compute their phase.
#'   The default, \code{FALSE} requests that these not be retained.
#'
#' @details
#' "Cleans" basin characteristics (BCs) by:
#' \itemize{
#' \item{1.}{\itemize{Removing BCs that are cannot be used because they are non-numeric
#' or, according to current thinking, should not be used because they:
#' \item{a)}{measure hydrologic modifications such as dams or land development}
#' \item{b)}{use NHD or TIGER roads which are spatially non-uniform
#' in their "resolution"}
#' \item{c)}{use measured streamflow (e.g, kriged streamflow statistics or
#' calibrated model output)}
#' \item{d)}{optionally, BCs based on non-calibrated hydrologic models}}}
#' \item{2.}{Removing BCs that have a large fraction, whose value is given by
#' \code{max.BC.oneval.frac},of values that are the same.}
#' }
#'
#' @return A matrix of basin characteristics, with several removed.
#'@export
clean_BCs <- function(BCs,BC.code6.remflg=F, max.BC.oneval.frac = 0.5,
  BC_suffix, destination="", keepWB=F) {
  # Function originally designed by Tom Over and Mike Olson, 30 June 2015.
  # Modified by William Farmer, 30 June 2015.
  # Revised by TMO, 7/2016 and Modified by Amy Russell, 9/2016
  #   to give option to keep Water Balance model variables (named "WB5100...")
  #   so that their phase could be computed later in winnow_BCs; also fixed a couple cosmetic items

  class_col = which(names(BCs)=="CLASS")
  col_order = c(class_col,(1:ncol(BCs))[-class_col])
  BCs = BCs[,col_order]

  #reading GAGESII.var_code_file commented out: GAGESII.var_codes now embedded in sysdata.rdata
  #GAGESII_BC_var_code_file <- file.path("Data","Raw",
  #  "selectedBCs.var_codes.csv")
  #Read in GAGES II variable code file
  #GAGESII.var_codes = read.csv(GAGESII_BC_var_code_file)
  BC.codes = GAGESII.var_codes$CODE
  BC.code.names = GAGESII.var_codes$VARIABLE_NAME
  nstatns = nrow(BCs)

  #Remove columns based on associated code
  #Remove columns with BC codes 2, 3, 4, and 5 plus 6 if BC.code6.remflg = T
  if (BC.code6.remflg) code.rem_BCs = BC.code.names[which(BC.codes>=2)] else
    code.rem_BCs = BC.code.names[which(BC.codes>=2 & BC.codes<6)]
  keptBCs = BCs[,!(names(BCs) %in% code.rem_BCs)]
  remBC_names = names(BCs)[names(BCs) %in% code.rem_BCs]
  remBC_codes = BC.codes[BC.code.names %in% remBC_names]
  remBC_frame = data.frame(cbind(remBC_names,remBC_codes))

  #Apply one-value fraction max.BC.oneval.frac to
  #remove variables with a large fraction of the same value
  BC.oneval.fracs = numeric(ncol(keptBCs)-1)
  remove.cols = NULL
  for (i in 2:ncol(keptBCs)) {
    #indexing starts at 2 because first column is supposed to be CLASS
    ranks = rank(keptBCs[,i])
    uniq.vals = unique(ranks)
    n.uvals = length(uniq.vals)
    uval.cnts = numeric(n.uvals)
    for (j in 1:n.uvals) uval.cnts[j] = sum(ranks==uniq.vals[j])
    BC.oneval.fracs[i] = max(uval.cnts)/nstatns

    # Added parameter keepWB to keep code6 (WB5100... Water Balance model discharge) variables
    #  even if some have a large fraction of a common value (would be 0 in this case),
    #  because after "winnowing" only the phase and the annual average will be retained -
    #  TMO, 7/2016
    if (keepWB==T){
      if ((BC.oneval.fracs[i] > max.BC.oneval.frac) & (substr(names(keptBCs)[i],1,6)!="WB5100")) {
        remove.cols = c(remove.cols,i)
      }
    } else {
      if (BC.oneval.fracs[i] > max.BC.oneval.frac) remove.cols = c(remove.cols,i)
    }
  }
  clean_BCs = keptBCs[,-remove.cols]

  if (destination!="") {
    if (BC.code6.remflg) code6_str = "rem_code6" else code6_str = "keep_code6"
    write.csv(keptBCs,paste0(destination,".keptBCs.",code6_str,".csv"))
    write.csv(remBC_frame,paste0(destination,".remBC_names.",
                                 code6_str,".csv"),row.names=F)
    write.csv(clean_BCs,paste(destination,".clean_BCs.",BC_suffix,".csv"))
    sink(paste(destination,".cleaned_BCs.",BC_suffix,".txt"))
    for (i in 1:length(remove.cols)) cat(names(keptBCs)[remove.cols[i]],"\n")
    sink()
  }
  return(clean_BCs)
}
