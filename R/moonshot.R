####
#### to build
####
if(FALSE){
# cd /home/ptong1/Backup/Package/; R0
library(roxygen2)
#library(roxygen) # not working
roxygenize("drexplorerExtra")

library(devtools)
build('drexplorerExtra')
install('drexplorerExtra')

devtools::load_all('/home/ptong1/Backup/Package/drexplorerExtra')


##
detach("package:drexplorerExtra", unload=TRUE)
library(drexplorerExtra)

    source(file.path(osbaseZdrive, 'ptong1/Backup/Package/drexplorerExtra/R/moonshot.R'))
}



#' split data from a plate into multiple experiment unit (one cell line one drug).
#'
#' This split function tries to split data from a plate which consists multiple drugs for the same cell line. The control data is shared across
#' multiple drugs. It returns a list of experiement unit which is dose response data for one drug and one cell line including the control data.
#'
#' This is written specific to the lung moonshot screening data where on one plate, only a single cell line is profiled with multiple drugs; The control
#' data is shared across different drugs. Currently it only splits on drug and append control data to each drug. 
#' We may later need to augment this function so that it splits on both cell line and drug.
#' @param dat a data frame
#' @param CL colname for cell line
#' @param Agent colname for drug name
#' @param Dose colname for dose
#' @param Response colname for response
#' @param ControlName name for the control data so that it will not be treated as a drug.
#' @param PlateID colname for PlateID. this is for tracking, not used for split
#' @param PlateName colname for PlateName, usually the run date. this is for tracking, not used for split
#' @param minReplicate if a subset does not have nrow exceeding this, throw an error
#' @return a list of data frame after split with the same columns as dat
#' @export
split_plate <- function(dat, CL='CellLine', Agent='Agent1', Dose='AgentConc1', Response='CellCount', ControlName='Control1',
			PlateID='PlateID', PlateName='PlateName', minReplicate=3){
	# only works on one cell line currently
	if(n_unique(dat[, CL])>1) stop('Multiple cell lines observed, this is not a typical plate data!\n')
	colNames <- colnames(dat)
	if(!CL %in% colNames) stop(sprintf('No column name for %s', CL))
	if(!Agent %in% colNames) stop(sprintf('No column name for %s', Agent))
	if(!Dose %in% colNames) stop(sprintf('No column name for %s', Dose))
	if(!Response %in% colNames) stop(sprintf('No column name for %s', Response))
	if(!ControlName %in% dat[, Agent]) stop(sprintf('The specified control name %s is not observed in the drug column!', ControlName))
	# split by Agent
	datsplit <- dlply(dat, Agent, function(x) x) # a list, each cell line per element
	#dfsplit_labels <- unname(unlist(attributes(datsplit)$split_labels)) # 
	dfsplit_labels <- names(datsplit) # 
	nMeasure <- sapply(datsplit, nrow) # number of rows for each drug-CL
	#browser()
	if(any(nMeasure<minReplicate)){
		indWrong <- which(nMeasure<minReplicate)
		msg <- sprintf('Cell line \n  %s\nhas too few measurements:\n  %s\n 
			Make sure you have the right data!\n', str_c(dfsplit_labels[indWrong], collapse=', '), str_c(nMeasure[indWrong], collapse=', '))
		cat('Error: Data is probably wrong!\n')
		cat('-----------------------\n')
		stop(cat(msg))
		cat('-----------------------\n')
	}
	dat_control <- datsplit[[ControlName]]
	datsplit[[ControlName]] <- NULL
	#browser()
	for(dd in names(datsplit)){
		datsplit[[dd]] <- rbind(dat_control, datsplit[[dd]])
	}
	datsplit
}
#split_plate(dat=dat_plate_list[[1]])
#tt <- split_plate(dat=L_by_CL_plate[[1]], CL='CellLine', Agent='Drug', Dose='Dose', Response='Count', ControlName='DMSO', PlateID='Plate.barcode', PlateName='PlateName', minReplicate=1)

extractResult <- function(ll, plot=FALSE){
	#browser()
	#ICtab <- ldply(ll, function(x) x$ICx)[, -c(1, 5)]
	ICtab <- ldply(ll, function(x) cbind(x$ICx, x$info_experimentSetup))
	if(plot){
		lapply(ll, function(x) plotOneExp(x, show='both'))
	}
	ICtab
}



#' fit model on one experiment and allows to plot a jpeg figure for the curves
#' 
#' this use fitOneExp() from drexplorer to analyze a lung moonshot experiement unit so that plate ID are tracked
#' 
#' @param x a data frame for data from one experiment
#' @param models models to be fitted; default is drexplroer recModels
#' @param CL colname for cell line
#' @param Agent colname for drug name
#' @param Dose colname for dose
#' @param Response colname for response
#' @param ControlName name for the control data so that it will not be treated as a drug.
#' @param PlateID colname for PlateID
#' @param unit unit of drug concentration
#' @param interpolation interpolation passed to fitOneExp
#' @param PlateName colname for PlateName, usually the run date
#' @param dirFigure folder path to store the figures. If set to NA, no figures will be output. With the returned list, it is easy to generate the figures though.
#'        the figures here are one drug/one cell line per file. 
#' @param info information passed from analyze_plate_list() about the plate
#' @param con connection to store log information; default is the console terminal; User can pass a file connection so that log will be saved in a file
#' @return a list of the fits for different models  
#' @export  
fit_one_experiment <- function(x, models=drModels(return=TRUE, verbose=FALSE)[['recommendedModels']], dirFigure=NA, CL='CellLine', Agent='Agent1', Dose='AgentConc1', 
		Response='CellCount', ControlName='Control1', PlateID='PlateID', PlateName='PlateName', debug=FALSE, unit=NA, interpolation=FALSE, con=stdout(), info=''){
	#cat(x[1, 1])
	tpPlateID <- x[1, PlateID]
	tpPlateName <- x[1, PlateName]
	tpCL <- x[1, CL]
	tpDrug <- setdiff(x[, Agent], ControlName)
	if(debug){
		cat(sprintf('Working on plate: %s, plateName: %s, cellLine: %s, drug: %s\n\t', tpPlateID, tpPlateName, tpCL, tpDrug))
	}
	info_experiment <- sprintf('%s; CellLine: %s, Drug: %s', info, tpCL, tpDrug)
	#browser()
	## add experimental info: N dose (non-control), N technical rep (number of obs. at each dose, may be a string), experimental variation including medianSD, minSD, maxSD, meanSD computed from SD across different tech replicates
	resPrepDR <- drexplorer:::prepDRdat(x[, c(Dose, Response)], alpha=0.05, fitCtr=FALSE, standardize=TRUE)
	info_experimentSetup <- drexplorer:::prepDRdat2expInfo(resPrepDR)
	if(!is.null(con))	cat(sprintf("Working on %s...", info_experiment), file = con)
	resL <- fitOneExp(x[, c(Dose, Response)], drug=tpDrug, cellLine=tpCL, models=models, plot=FALSE, transparency=0.95, interpolation=interpolation, unit=unit)
	resL$ICx <- data.frame(PlateID=tpPlateID, PlateName=tpPlateName, resL$ICx)
	resL$info_experimentSetup <- info_experimentSetup
	#browser()
	if(!is.na(dirFigure)){
		#browser()
		# a cell line may named: H69/CR, thus corrupting the name!
		#pathFigure <- file.path(dirFigure, sprintf('plateID%s_plateName%s_%s_%s.jpeg', tpPlateID, tpPlateName, tpCL, tpDrug))
		pathFigure <- file.path(dirFigure, make.names(sprintf('plateID%s_plateName%s_%s_%s.jpeg', tpPlateID, tpPlateName, tpCL, tpDrug)))
		jpeg(file=pathFigure, quality=100, pointsize=38, width = 480*3.5, height = 480*3.5) # increase width and height for better resolution
		plotOneExp(resL, show='both', lwd=5)
		dev.off()
		#browser()
	}
	#browser()
	if(!is.null(con))	cat(sprintf("Done! resL$ICx[1, 'AUCs']=%.3f\n", resL$ICx[1, 'AUCs']), file = con)
	resL
}

#' analyze data from a plate 
#' 
#' this wraps over fit_one_experiment
#' 
#' @param plate a list of data from one plate 
#' @param dirFigure where to save figure; default is no figure to save
#' @param models models to be fitted; default is drexplroer recModels
#' @param unit unit of drug concentration
#' @param info information passed from analyze_plate_list() about the plate
#' @param con connection to store log information; default is the console terminal; User can pass a file connection so that log will be saved in a file
#' @param ... additional parameters passed to split_plate such as CL='CellLine', Agent='Agent1', Dose='AgentConc1', Response='CellCount', ControlName='Control1'
#' @return a list of the fits  
#' @export
analyze_one_plate <- function(plate, dirFigure=NA, models=drModels(return=TRUE, verbose=FALSE)[['recommendedModels']], unit=NA, interpolation=FALSE, con=stdout(), info='', returnData=FALSE, ...){
	dat_experiment_list <- split_plate(dat=plate, ...)
	info_1plate <- sprintf('%s; %d experiments', info, length(dat_experiment_list))
	#browser()
	resL <- llply(dat_experiment_list, fit_one_experiment, dirFigure=dirFigure, models=models, unit=unit, interpolation=interpolation, info=info_1plate, con=con) # figures plotted
	tab <- extractResult(resL) # now a tab of IC for all drugs
	if(returnData){
		return(list(resL=resL, ICtab=tab, dataList=dat_experiment_list))	
	}
	else {
		return(list(resL=resL, ICtab=tab)) # returns the fit as well as IC tab; the fit can be used to generate other figures
	}
}


#' analyze a list of plate data
#'
#' this wraps over analyze_one_plate and use parallel computing within plate
#' paralell is performed per plate as an unit
#'
#' @param  plate_list a list of plate data
#' @param dirFigure where to save figure; default is no figure to save
#' @param parallel logical
#' @param models models to be fitted; default is drexplroer recModels
#' @param core number of cores to register for parallel computing
#' @param unit unit of drug concentration
#' @param con connection to store log information; default is the console terminal; User can pass a file connection so that log will be saved in a file
#'  Notice that con is only possible if not using parallel (foreach does not allow multiple cores to simultaneously writting on the same file, as connection cannot be exported to different node)
#' @param ... additional parameters for colname tracking
#' @return a list
#' @export
analyze_plate_list <- function(plate_list, dirFigure=NA, core=NA, parallel=TRUE, models=drModels(return=TRUE, verbose=FALSE)[['recommendedModels']], unit=NA, interpolation=FALSE, con=stdout(), ...){
	#browser()
	nJobs <- length(plate_list)
	if(parallel){
		con <- NULL # when parallel, multiple cores cannot write on the same file
		if(nJobs<=10) {
			core <- nJobs
		} else {
			core <- 10
		}
		# required variables (including functions) and packages are exported here
		cl <- createCluster(core=core, logfile = "/dev/null", export=c('analyze_one_plate'), lib = c('drexplorer', 'drexplorerExtra', 'plyr'))
		on.exit(stopCluster(cl))
	}
	#browser()
	ll <- foreach(i=1:length(plate_list)) %dopar% {
		plate <- plate_list[[i]]
		pid <- as.character(unique(plate[, 1])[1])
		info_plateL <- sprintf('Plate index: %d; PlateID: %s', i, pid)
		if(!is.null(con))	cat(sprintf(">>>>> Starting new plate: Plate ID=%s\n", pid), file = con)
		#browser()
		res_onePlate <- analyze_one_plate(plate, dirFigure=dirFigure, models=models, unit=unit, interpolation=interpolation, info=info_plateL, con=con, ...)
		if(!is.null(con))	cat(sprintf(">>>>> >>>>> Plate ID=%s succeeded; nrow(res_onePlate$ICtab)=%d; length(res_onePlate$resL)=%d\n", pid, nrow(res_onePlate$ICtab), length(res_onePlate$resL)), file = con)
		res_onePlate	
	}
	#browser()	
	if(!is.null(con))	cat(sprintf(">>>>> >>>>> >>>>>> All plates done for %d plates; foreach i=%d; foreach output length(ll)=%d\n", length(plate_list), i, length(ll)), file = con)
	ll
}
#time2 <- system.time(ll <- analyze_plate_list(dat_plate_list, dirFigure=file.path(WD, 'Figure_DR_curves'), models=moonshotModels, 
#		CL=CL_name, Agent=Agent_name, Dose=Dose_name, Response=Response_name, ControlName=Control_name, PlateID=PlateID_name, PlateName=Plate_name, unit=UNIT, parallel=TRUE))

#ll <- analyze_plate_list(dat_plate_list, dirFigure=file.path(WD, 'Figure_DR_curves'), models=drModels(return=TRUE, verbose=FALSE)[['recommendedModels']], 
#		CL=CL_name, Agent=Agent_name, Dose=Dose_name, Response=Response_name, ControlName=Control_name, PlateID=PlateID_name, PlateName=Plate_name, unit=unit, parallel=FALSE)


# extract SE of say IC50 from results of fit_one_experiment
extractICse <- function(resL1, whichSE=c('IC50')) {
	#resL1 <- fit$resL
	RSEs <- resL1$RSEs
	ICmat_variance <- resL1$ICmat_variance
	#browser()
	seVec <- sqrt(ICmat_variance[, whichSE])
	# extract feasible SE: not NA but with least RSE
	isNotNA <- !is.na(seVec)
	td <- data.frame(isNotNA=isNotNA, SE=seVec, RSE=RSEs)
	td <- td[order(td$isNotNA*(-1), td$RSE, decreasing=FALSE), ]
	td1 <- subset(td, isNotNA==TRUE)
	if(nrow(td1)>0){
		res <- td1[1, 'SE']
	} else {
		res <- NA
	}
	# if(sum(isNotNA)>0){
	# 	res <- min(seVec, na.rm=TRUE)
	# } else {
	# 	res <- NA
	# }
	names(res) <- str_c(whichSE, '_se', sep='')
	res
}
# extract multiple SEs
extractICse_vec <- function(fit, whichSEs='IC50'){
	res <- sapply(whichSEs, function(x) extractICse(fit, whichSE=x))
	names(res) <- str_c(whichSEs, '_se', sep='')
	res
}

#' parse the result from analyze_plate_list
#'
#' the result is a nested list; on the top level, we have resL and ICtab; resL is another list each element is the fit for one drug; ICtab is a table for 
#' the IC values for all experiments in a plate.
#'
#' @param ll list as returned by analyze_plate_list
#' @param type whether to extract IC table or a list of fit for plotting
#' @param whichSE a vector of IC names for SE values, i.e. s.e. for IC50 by whichSE=c('IC50')
#' @return a table or list of fit
#' @export
parse_res_from_plate_list <- function(ll, type=c('ICtable', 'fitList'), whichSE='IC50'){
	if(type == 'ICtable') {
		ICtable <- foreach(l=ll, .combine='rbind') %do% {
			#browser()
			#extractICse(l$resL[[1]], whichSE='IC50')	
			if(!is.null(whichSE)){
				dfSE <- data.frame(sapply(l$resL, extractICse, whichSE=whichSE))
				colnames(dfSE) <- str_c(whichSE, '_se', sep='')
				rownames(dfSE) <- NULL
				ICtable <- data.frame(l$ICtab, dfSE)
			} else {
				ICtable <- l$ICtab
			}
			
		}
		ICtable$fitIndex <- 1:nrow(ICtable) # this index can be used to retrieve the fits
		return (ICtable)
	}	
	if(type == 'fitList') {
		fitList <- vector('list')
		for(l in ll) {
			list_all_drugs <- l$resL
			for(fit in list_all_drugs){
				fitList <- lappend(fitList, fit)
			}	
		}
		return (fitList)
	}
	#browser()
}


#' parse the result from fit_one_experiment
#'
#' the result is a list
#' @param whichSE a vector of IC names for SE values, i.e. s.e. for IC50 by whichSE=c('IC50')
#' @return a table or list of fit
#' @export
parse_res_from_exp_fit <- function(fit, whichSE='IC50'){
	ICtab <- cbind(fit$ICx, fit$info_experimentSetup)
	dfSE <- extractICse(fit, whichSE=whichSE)
	#browser()
	if(!is.null(whichSE)){
		ICtable <- data.frame(ICtab, dfSE)
		colnames(ICtable)[(ncol(ICtab)+1):ncol(ICtable)] <- names(dfSE)
	} else {
		ICtable <- data.frame(ICtab)
	}
	ICtable
}
#parse_res_from_exp_fit(fit)

#' prune a vector of IC values according to minDose and maxDose
#'
#' This is useful if we want to truncate non-achievable IC values (i.e. -30, 30) with observed minimum and maximum dose
#'
#' @param vec a vector of IC values
#' @param minD minimum dose adminised
#' @param maxD maximum dose adminised
#' @return a vector of the pruned IC values
#' @export
pruneICvec <- function(vec, minD, maxD){
	truncByLimit(vec, Lower=minD, Upper=maxD)
}
#' a wrapper of pruneICvec to prune a IC table
#' @param tab a table containing IC values and possibly other info
#' @param IC_names colnames for IC values
#' @param minD_name colname for minimum administered dose
#' @param maxD_name colname for maximum administered dose
#' @return a pruned table
#' @export
pruneICtable <- function(tab, IC_names, minD_name, maxD_name){
	res <- foreach(i=1:nrow(tab), .combine='rbind') %do% {
		tab[i, IC_names] <- pruneICvec(tab[i, IC_names], minD=tab[i, minD_name], maxD=tab[i, maxD_name])
		tab[i, ]
	}
	#browser()
	# IC50_se: if IC50 not achieved, SE is problematic, so mask it with NA
	col_se <- str_detect(colnames(tab), '_se')
	if(sum(col_se)>0){
		colname_se <- colnames(tab)[col_se]
		icname_se <- str_split(colname_se, '_')[[1]][1]
		#ind_infeasible <- tab[, icname_se] - tab[, 'maxLog10Dose'] >=0 | tab[, icname_se] - tab[, 'minLog10Dose'] <=0
		ind_infeasible <- tab[, icname_se] - tab[, 'maxLog10Dose'] >=0  # this maybe enough
		IC_se_mask <- tab[, colname_se]
		IC_se_mask[ind_infeasible] <- NA
		res$IC_se_mask <- IC_se_mask
		colnames(res)[ncol(res)] <- str_c(colname_se, '_pruned', sep='')
	}
	rownames(res) <- rownames(tab)
	res
}	