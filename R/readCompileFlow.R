#readCompileFlow----
#' Function to read and compile streamflow data for analysis
#'
#' @description
#' The function \code{readCompileFlow} collects the daily streamflow information
#' available for a specified \code{listofgages}.
#'
#' @param listofgages Dataframe of gage numbers (as a character vector)
#' @param units Specifies the units of the output.  Options include "cms" for
#' cubic meters per second and "cfs" for cubic feet per second
#' @param dataset_name Allows the user to specify a specific output location.
#' Output will be saved as an R workspace.  If no name is given, output data
#' will be returned but not saved.
#' @param startDate Specifies a particular date from which to begin retrieving
#' the streamflow data.  Format should be "YYYY-MM-DD".  If blank, all available
#' data will be retrieved.
#' @param endDate Specifies a particular date up to which to retrieve
#' the streamflow data.  Format should be "YYYY-MM-DD".  If blank, all available
#' data will be retrieved.
#' @param checkmiss Logical.  If \code{TRUE}, if any day is missing a
#' streamflow value for a given site, none of the data for that site will be
#' collected.  If \code{FALSE}, missing days will be allowed.
#' @param keepNAdays Logical.If \code{TRUE}, if there is a day of flow missing,
#'  the date will be included with an \code{NA} streamflow value. If
#'  \code{FALSE}, only dates with streamflow data will be collected.  This input
#'   is only valid if \code{checkmiss=FALSE}.
#'
#' @details
#' Lorem ipsum...
#'
#' @return The function returns a zoo object of observed streamflows.
#'@export
#'@import dataRetrieval
#'@import zoo
readCompileFlow<-function(listofgages,units=c("cfs","cms"),dataset_name="",
  startDate="",endDate="",checkmiss=FALSE,keepNAdays=TRUE){
  # Function orginially designed by Stacey A. Archfield, 02 June 2015.
  # Modified by William Farmer, 03 June 2015.

  # @importFrom dataRetrieval readNWISdv
  # @importFrom  zoo zoo

  # Read in stations from input file; no header
  stations<-listofgages
  colnames(stations)<-c("gage")
  stations$gage<-as.character(stations$gage)

  # Create vector to store revised list of stations which excludes
  #   any stations that are removed
  rev.stations<-NULL

  # Get total number of gages
  M<-length(stations$gage)

  # Parameter code for discharge
  parameterCd<-"00060"

  # If user asks for output to be saved.
  if (dataset_name!="") {
    # Filename to store flow data for all sites
    all_q_filenm<-paste(dataset_name,"_",units,".RData",sep="")
    # Remove any previous versions of this file so it
    #   does NOT append this new data to an existing file
    if (file.exists(all_q_filenm)){
      file.remove(all_q_filenm)
    }
  }

  returnList <- list()
  AllFlow <- list()
  count <- 0

  for (i in 1:M) {

    # Get additional gage information
    #qinfoi<-getSiteFileData(stations$gage[i])

    # Get all available streamflow data for gage i, if streamgage exists
    possibleError <- tryCatch(
      qdatai<-readNWISdv(stations$gage[i],parameterCd,startDate,endDate),
      error=function(e) e
    )

    if(!inherits(possibleError, "error")){

      returnList[[length(returnList)+1]] <- qdatai

      # Check that there are more than 1 measurement in qdatai
      if (nrow(qdatai)!=0) { #There is more than 1 measurement

        #Rename columns
        colnames(qdatai)[1:3]<-c("agency","site_no","datetime")
        colnames(qdatai)[
          which(grepl('00060_00003_cd',colnames(qdatai)))] <-
          "approv" # Correction added by William Farmer, 04 Feb 2016
        colnames(qdatai)[
          which(grepl('00060_00003',colnames(qdatai)))] <-
          "qval" # Correction added by William Farmer, 04 Feb 2016
        # Actual columns names from NWIS:
        #   agency_cd  site_no  Date  X_00060_00003_cd	X_00060_00003

        # Capture only approved data
        qdataapprovi<-subset(qdatai,(qdatai$approv=="A" | qdatai$approv=="A e"))

        # Define start and end to record to check
        #   that record covers entire period
        if (startDate=="") {
          startDate<-qdataapprovi$datetime[1]
        } else {
          startDate<-startDate
        }
        if (endDate=="") {
          #Compute numnber of days spanning beginning and end of record
          numpts<-length(qdataapprovi$datetime)
          endDate<-qdataapprovi$datetime[numpts]
        } else {
          endDate<-endDate
        }

        # Compute number of days of observed
        daysinperiod<-length(qdataapprovi$datetime)

        if (checkmiss==TRUE){
          # Now compute # of expected complete days in the record
          completedays<-as.numeric(as.Date(endDate)-as.Date(startDate)+1)
        } else {
          completedays<-daysinperiod
        }

        # If there are missing days for the period of interest
        #   (specified by startDate and endDate), then remove the
        #   site from the analysis; UNLESS checkmiss=FALSE
        if (completedays==daysinperiod) {

          if (keepNAdays==TRUE) {
            # If desire missing streamflow values to be included as NAs
            seqdate<-data.frame(seq(as.Date(startDate),
              as.Date(endDate), by="day"))
            colnames(seqdate)<-c("datetime")
            # Make sure dates are formatted as dates
            seqdate$datetime<-as.Date(seqdate$datetime, "%Y-%m-%d")
            qdataapprovi$datetime<-as.Date(qdataapprovi$datetime, "%Y-%m-%d")
            # Merge alldates with days of data
            merge_dates<-merge(seqdate,qdataapprovi,by="datetime",
              all.x=T,sort=T)
            # Create vector of site_no and add to merge_dates
            merge_dates$site_no<-qdataapprovi$site_no[1]
            qonlyi<-data.frame(merge_dates$site_no,merge_dates$datetime,
              merge_dates$qval)
          } else {
            # Obtain only fields with site name, date, and discharge
            #   and rename columns
            qonlyi<-data.frame(qdataapprovi$site_no,qdataapprovi$datetime,
              qdataapprovi$qval)
          }
          colnames(qonlyi)<-c("gage","date","flow")

          # Store raw dates without formatting
          rawdates<-qonlyi$date
          # Format date in Q file
          dateaschar<-as.character(rawdates)
          qonlyi$date<-as.Date(dateaschar, "%Y-%m-%d")

          colnames(qonlyi)<-c("gage","date","flow")
          if (units=="cms"){
            # Convert flow from cu.ft/sec to cu.m/sec
            qonlyi$flow<-qonlyi$flow*0.02831685
          }
          count <- count+1

          # Report progress
          cat("\n Successfully processed gage",stations$gage[i],
            "- Done with gage",i," of ",M)
          rev.stations<-c(rev.stations,stations$gage[i])
          # Make zoo object to store all streamflow data
          AllFlow[[count]] <- zoo(qonlyi$flow, qonlyi$date)
          names(AllFlow)[[count]] <- stations$gage[i]

        } else {
          #Write message that the gage was removed
          cat("\n Gage",stations$gage[i],
            "removed due to missing record for the period",
            "of interest. Done with gage",i,"of",M)
        }

      } else {
        # Less than two streamgage measurements exist.
        #   Write message that the gage was removed
        cat("\n Gage",stations$gage[i],
          "removed due to missing record for the",
          "period of interest. Done with gage",i,"of",M)

      }

    } else {
      cat("\n Gage",stations$gage[i],
        "not found. Removed due to gage not being available in NWIS.\n")
    }

  }

  # Output revised gage list with any removed sites
  if (dataset_name!="") {
    rev.stationsf<-data.frame(rev.stations)
    names(rev.stationsf)<-names(listofgages)
    if (!identical(rev.stationsf,listofgages)) {
      names(rev.stationsf)<-c("gage")
      write.table(rev.stationsf,paste(dataset_name,
        "_revlistofgages.txt",sep=""))
    }
    #Save zoo object as compressed .RData
    save(list=c("AllFlow"), file=all_q_filenm, compress="bzip2")
  }

  return(AllFlow)
}
