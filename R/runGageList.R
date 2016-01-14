# #
# # William Farmer
# #
# #
# rm(list=ls())
# HUC02 <- 9
#
# # Define Inputs ####
# GageList.File <-file.path('Data','Raw',paste0("region",HUC02,
#   ".refGages.txt")) # Txt file of gages to be used
# SaveFile <- file.path('Results',paste0('region',HUC02,
#   '.result.RData')) # File for saving final results.
# regSelect <- 'adjR2'
# logFile <- paste0("UserLogs/LogFile.AirDrop.region",HUC02,".txt")
# # Optional
# BasChar.File <-
#   'Data/Raw/BasChars.txt' # Not needed unless you want specific variables.
# # File indicating which variables to pull, if empty, all
# Start.Date <- "1959-10-01"
# Start.Date <- "1980-10-01"
# End.Date <- "2013-09-30"
# WY.Lim <- 10 # minimum number of water years
# CV.fold <- NA # How many folds to randomly cross-validate; NA is LOO
# CV.dense <- T # Which side of the CV folds to use.
# FlowStat <- T # In analysis, do you want flow statistics calculated
# dataset_name <- file.path("Data","Processed",
#   paste0("region",HUC02,".Streamflow"))
# saveName <- file.path("Data","Processed",
#   paste0("region",HUC02,".FDC"))
# useOKLUD <- T # Whether or not to apply OKLUD
# DAR.low <- -Inf
# DAR.high <- Inf
#
# # Get all streamflow data ####
# listofgages<-read.table(GageList.File,header=T,colClasses=c("character"))
# if (ncol(listofgages)>1) {
#   listofgages <- data.frame(listofgages[,2])
# }
# names(listofgages) <- "gages"
# Streamflow <- readCompileFlow(listofgages,units="cfs",
#   dataset_name=dataset_name,
#   startDate=Start.Date,endDate=End.Date,checkmiss=FALSE,keepNAdays=TRUE)
#
# # Get basin characteristics ####
# #basinChars <- as.character(read.table(BasChar.File,sep="\t")$V1)
# expVars <- getBasinChar(listofgages)
#
# # Produce empirical FDC ####
# empFDC.calc <- calcEmpFDCs(flow_data=Streamflow,
#   cens='filled')
# empFDCs.all <- empFDC.calc$empFDC
# p <- as.double(row.names(empFDCs.all))
# FDC.sites <- as.double(colnames(empFDCs.all))
# ndx <- match(as.double(as.character(listofgages$gages)),FDC.sites)
# empFDC.all.filled <- empFDCs.all[,ndx]
#
# empFDC.calc <- calcEmpFDCs(flow_data=Streamflow,
#   saveName=saveName)
# empFDCs.all <- empFDC.calc$empFDC
# p <- as.double(row.names(empFDCs.all))
# FDC.sites <- as.double(colnames(empFDCs.all))
# ndx <- match(as.double(as.character(listofgages$gages)),FDC.sites)
# empFDC.all.unfilled <- empFDCs.all[,ndx]
# FDC.sites.unfilled <- FDC.sites
# empFDC.calc.unfilled <- empFDC.calc
#
# # Screen for complete water years above a limit ####
# listofgages.new <- listofgages
# ndx <- which(empFDC.calc$PORstats.compWyr[,2]>=WY.Lim)
# ndx2 <- which(is.element(as.double(as.character(listofgages$gages)),
#   FDC.sites[ndx]))
# listofgages.new <- data.frame(listofgages[ndx2,])
# names(listofgages.new) <- names(listofgages)
# listofgages.old <- listofgages
# listofgages <- listofgages.new
#
# # create CV scenarios (LOO default) ####
# if (is.na(CV.fold)|CV.fold==1|CV.fold>length(listofgages$gages)) {
#   CV.fold=length(listofgages$gages)
# }
# if (CV.fold == length(listofgages$gages)) {
#   scenNums <- 1:length(listofgages$gages)
#   CV.dense <- TRUE
# } else {
#   sampSize <- round(length(listofgages$gages)/CV.fold)
#   scenNums <- rep(NA,length(listofgages$gages))
#   for (i in 1:(CV.fold-1)) {
#     ndx <- sample(which(is.na(scenNums)),sampSize)
#     scenNums[ndx] <- i
#   }
#   scenNums[is.na(scenNums)] <- CV.fold
# }
#
# #cl=makeCluster(detectCores()-1)
# # cl=makeCluster(2)
# cl=makeCluster(100,outfile="")
# registerDoParallel(cl)
#
# # logFile <- "Test.txt"
#
# # Conduct CV
# cat(paste("Begin CV -",Sys.time(),"\n"),file=logFile)
# # cvList <- list()
# # for (i in 1:CV.fold)
# cvList <-
#   foreach(i = icount(CV.fold),
#     #.export=as.vector(lsf.str()),
#     .packages=c("zoo","geoR","leaps","smwrQW","lubridate","survival")
#   ) %dopar% {
#     stime <- Sys.time()
#     # Setup index set and target set
#     if (CV.dense) {
#       index.gages <- as.double(as.character(listofgages$gages[scenNums!=i]))
#       target.gages <- as.double(as.character(listofgages$gages[scenNums==i]))
#     } else {
#       index.gages <- as.double(as.character(listofgages$gages[scenNums==i]))
#       target.gages <- as.double(as.character(listofgages$gages[scenNums!=i]))
#     }
#     # Get target information
#     ndx <- match(target.gages,expVars$DescBCs$STAID)
#     target.baschar <- expVars$DescBCs[ndx,]
#     target.regvar <- expVars$transf$transfBCs[ndx,]
#     ndx <- match(target.gages,as.double(names(Streamflow)))
#     target.obs <- Streamflow[ndx]
#     ndx <- match(target.gages,as.double(colnames(empFDC.all.unfilled)))
#     target.empFDC.unfilled <- empFDC.all.unfilled[,ndx]
#     ndx <- match(target.gages,as.double(colnames(empFDC.all.filled)))
#     target.empFDC.filled <- empFDC.all.filled[,ndx]
#
#     # Get index information
#     ndx <- match(index.gages,expVars$DescBCs$STAID)
#     index.baschar <- expVars$DescBCs[ndx,]
#     index.regvar <- expVars$transf$transfBCs[ndx,]
#     ndx <- match(index.gages,as.double(names(Streamflow)))
#     index.obs <- Streamflow[ndx]
#     ndx <- match(index.gages,as.double(colnames(empFDC.all.unfilled)))
#     index.empFDC.unfilled <- empFDC.all.unfilled[,ndx]
#     ndx <- match(index.gages,as.double(colnames(empFDC.all.filled)))
#     index.empFDC.filled <- empFDC.all.filled[,ndx]
#     ndx <- match(index.gages,FDC.sites.unfilled)
#     index.compwys <- empFDC.calc.unfilled$PORstats.compWyr[ndx,2]
#
#     # Develop regressions ####
#     leaps.models <- compute.leaps.for(unfilled_FDCs=index.empFDC.unfilled,
#       filled_FDCs=index.empFDC.filled,
#       comp.wys=index.compwys,
#       expl=index.regvar,
#       WY.lim=WY.Lim)
#
#     top.models <- compile.vars(leaps_list=leaps.models)
#
#     best.mods <- best.models(top_n_list=top.models,
#       unfilled_FDCs=index.empFDC.unfilled,
#       expl=index.regvar,comp.wys=index.compwys,WY.lim=WY.Lim)
#
#     # Estimate ungauged FDC ####
#     target.estFDC <- as.matrix(target.empFDC.filled)
#     target.estFDC[] <- NA
#     js <- 0
#     regimes <- best.mods$regimes
#     for (r in 1:length(regimes)) {
#       j0 <- max(js+1)
#       j1 <- j0 + length(regimes[[r]]) - 1
#       js <- j0:j1
#       if (regSelect=="adjR2") {# Get best adjusted R2.
#         m <- 1
#       } else if (regSelect=="AIC") {# Get best AIC.
#         m <- 4
#       }
#       mod.mat <- best.mods$best[[r]][[m]]
#       if (sum(is.na(as.double(mod.mat[nrow(mod.mat)-1,-1])))>0
#         & regSelect=="adjR2") {
#         # If no valid model, switch to AIC.
#         m <- 4
#       } else if (sum(is.na(as.double(mod.mat[nrow(mod.mat),-1])))>0
#         & regSelect=="AIC") {
#         # If no valid model, switch to adjusted R2
#         m <- 1
#       }
#       mod.mat <- best.mods$best[[r]][[m]]
#       ndx <- !is.element(mod.mat[,1],c("SD","VIF","R^2","AIC"))
#       mod.dat <- data.frame(t(apply(mod.mat[ndx,-1],2,as.double)))
#       names(mod.dat) <- mod.mat[ndx,1]
#       ndx <- which(is.element(names(target.regvar),names(mod.dat)))
#       r.target.regvar <- target.regvar[,ndx]
#       target.estFDC[js,] <- mod.dat[,1] +
#         as.matrix(mod.dat[,-1])%*%t(as.matrix(r.target.regvar))
#     }
#     target.estFDC <- as.numeric(10^target.estFDC)
#     dim(target.estFDC) <- dim(target.empFDC.filled)
#     names(target.estFDC) <- names(target.empFDC.filled)
#
#     # Determine Nearest Neighbor ####
#     NN.ranking <- index.NN(index.gages,index.baschar,
#       target.gages,target.baschar)
#
#     # Determine Map Correlation ####
#     MC.ranking <- index.MC(index.gages,index.obs,index.baschar,
#       target.gages,target.obs,target.baschar,
#       method="pearson")
#
#     # Apply DAR ####
#     # NN-DAR
#     #     cat(paste("CV.fold",i,"- NNDAR -",Sys.time(),"\n"),
#     #       file=logFile,append=TRUE)
#     NNDAR <- estDAR(index.network=NN.ranking,index.baschar,
#       target.baschar,index.obs,target.obs,DAR.low=DAR.low,DAR.high=DAR.high)
#
#     # MC-DAR
#     #     cat(paste("CV.fold",i,"- MCDAR -",Sys.time(),"\n"),
#     #       file=logFile,append=TRUE)
#     MCDAR <- estDAR(index.network=MC.ranking,index.baschar,
#       target.baschar,index.obs,target.obs,DAR.low=DAR.low,DAR.high=DAR.high)
#
#     # Apply QPPQ ####
#     # NN-QPPQ
#     #     cat(paste("CV.fold",i,"- NNQPPQ -",Sys.time(),"\n"),
#     #       file=logFile,append=TRUE)
#     NNQPPQ <- estQPPQ(index.network=NN.ranking,index.obs,index.empFDC.filled,
#       zero.val=NA,target.obs,target.empFDC.filled,target.estFDC,pvals=p)
#
#     # MC-DAR
#     #     cat(paste("CV.fold",i,"- MCQPPQ -",Sys.time(),"\n"),
#     #       file=logFile,append=TRUE)
#     MCQPPQ <- estQPPQ(index.network=MC.ranking,index.obs,index.empFDC.filled,
#       zero.val=NA,target.obs,target.empFDC.filled,target.estFDC,pvals=p)
#
#     if (useOKLUD) {
#       # Apply OKDAR ####
#       #       cat(paste("CV.fold",i,"- OKDAR -",Sys.time(),"\n"),
#       #         file=logFile,append=TRUE)
#       OKDAR <- krigeLUD(index.gages,index.baschar,index.obs,
#         target.gages,target.baschar,target.obs)
#       parres <- list(NNDAR=NNDAR,MCDAR=MCDAR,OKDAR=OKDAR,
#         NNQPPQ=NNQPPQ,MCQPPQ=MCQPPQ)
#     } else {
#       parres <- list(NNDAR=NNDAR,MCDAR=MCDAR,
#         NNQPPQ=NNQPPQ,MCQPPQ=MCQPPQ)
#     }
#     ltime <- Sys.time()
#     etime <- ltime - stime
#     cat(paste0("CV.fold ",i,": Started ",stime,"; Completed ",
#       ltime,"(",round(etime,digits=2),")\n"),
#       file=logFile,append=TRUE)
#     #cvList[[i]] <- parres
#     return(parres)
#   }
# for (i in 1:length(cvList)) {
#   names(cvList)[i] <- names(cvList[[i]]$NNDAR)
# }
# cat(paste("All Folds Estimated -",Sys.time(),"\n"),file=logFile,append=TRUE)
#
# # Save CV results
# inputs <- list(GageList.File=GageList.File,BasChar.File=BasChar.File,
#   Start.Date=Start.Date,End.Date=End.Date,WY.Lim=WY.Lim,CV.fold=CV.fold,
#   CV.dense=CV.dense,FlowStat=FlowStat,dataset_name=dataset_name,
#   saveName=saveName,useOKLUD=useOKLUD,DAR.low=DAR.low,DAR.high=DAR.high,
#   regSelect=regSelect)
# save(list=c("inputs","cvList"),file=SaveFile)
#
#
# # Unlist objects
# NNDAR.All <- eval(parse(text=
#     paste0("c(",paste0("cvList[[",c(1:length(cvList)),
#       "]]$NNDAR",collapse=","),")")))
# MCDAR.All <- eval(parse(text=
#     paste0("c(",paste0("cvList[[",c(1:length(cvList)),
#       "]]$MCDAR",collapse=","),")")))
# NNQPPQ.All <- eval(parse(text=
#     paste0("c(",paste0("cvList[[",c(1:length(cvList)),
#       "]]$NNQPPQ",collapse=","),")")))
# MCQPPQ.All <- eval(parse(text=
#     paste0("c(",paste0("cvList[[",c(1:length(cvList)),
#       "]]$MCQPPQ",collapse=","),")")))
# if (useOKLUD) {
#   OKDAR.All <- eval(parse(text=
#       paste0("c(",paste0("cvList[[",c(1:length(cvList)),
#         "]]$OKDAR$OKLUD",collapse=","),")")))
#   OKDAR.Var <- eval(parse(text=
#       paste0("list(",paste0("cvList[[",c(1:length(cvList)),
#         "]]$OKDAR$varPar",collapse=","),")")))
#   sets <- c("NNDAR","MCDAR","OKDAR","NNQPPQ","MCQPPQ")
#   exList <- c("analyzeResult","FlowStat","nse","rmse.like","percent.error",
#     "obs.sim.corr","getFlowStats","calcFlowStats",
#     "NNDAR.All","MCDAR.All","NNQPPQ.All","MCQPPQ.All","OKDAR.All")
# } else {
#   sets <- c("NNDAR","MCDAR","NNQPPQ","MCQPPQ")
#   exList <- c("analyzeResult","FlowStat","nse","rmse.like","percent.error",
#     "obs.sim.corr","getFlowStats","calcFlowStats",
#     "NNDAR.All","MCDAR.All","NNQPPQ.All","MCQPPQ.All")
# }
#
# # Parallel analysis
# cat(paste0("Started Analysis - ",Sys.time(),"\n"),file=logFile,append=TRUE)
# # AnList <- list()
# # for (i in 1:length(sets))
# AnList <-
#   foreach(i = icount(length(sets)),
#     .export = exList,
#     .packages = c("lmomco","zoo")
#   ) %dopar% {
#     stime<-Sys.time()
#     result <- eval(parse(text=
#         paste0("analyzeResult(",sets[i],".All,FlowStat=FALSE)")
#     ))
#     ltime <- Sys.time()
#     etime <- ltime - stime
#     cat(paste0("Completed analysis of ",sets[i],": Started ",stime,
#       "; Completed ",ltime,"(",round(etime,digits=2),")\n"),
#       file=logFile,append=TRUE)
#     #AnList[[i]] <- result
#     return(result)
#   }
# names(AnList) <- sets
# cat(paste0("Completed all analysis - ",Sys.time(),"\n"),
#   file=logFile,append=TRUE)
# stopCluster(cl)
#
# # Save CV results
# save(list=c("inputs","AnList","cvList"),file=SaveFile)
#
# a <- AnList[[1]]$PerfMat$nse
# for (i in 2:length(AnList)) {
#   a <- cbind(a,AnList[[i]]$PerfMat$nse)
# }
