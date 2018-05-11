####
#### to build
####
if(FALSE){
# cd /home/ptong1/Backup/Package/
# Z:\Backup\Package\drexplorerExtra\R\base.R

#
# for drexplorer: cd /data/bioinfo2/ptong1/Projects/Coombes/IC50Package/Package/
#
library(roxygen2)
#library(roxygen) # not working
roxygenize("drexplorerExtra")

library(devtools)
build('drexplorerExtra')
install('drexplorerExtra')

load_all('drexplorerExtra')


##
detach("package:drexplorerExtra", unload=TRUE)
library(drexplorerExtra)

osbaseZdrive <- ifelse(.Platform$OS.type=="windows", "//mymdafiles/usersdqs1", "/home") 

    source(file.path(osbaseZdrive, 'ptong1/Backup/Package/drexplorerExtra/R/moonshot.R'))
    source(file.path(osbaseZdrive, 'ptong1/Backup/Package/drexplorerExtra/R/base.R'))

	
}


# maxum response value to control ylim
getMaxResponse <- function(fits){
	max(sapply(fits, function(x) max(x@fitDat[, 2]))) # maxum response value to control ylim
}

#' shadow plot for multiple fits
#' 
#' this plots multiple curves from a set of experiments
#'
#' when passing a list with one fit object, this can be used to plot the best model for a experiment
#'
#' @param fitResList a list of fit result from multiple experiments 
#' @param indFocus which experiment to highlight and other are treated as background shadow
#' @param h horizontal line added to the figure, i.e. indicating IC50, IC70
#' @param main main
#' @param cex.main cex.main to adjust main title size
#' @param style either 'full' or 'simple' indicating if data points to be added. 
#' @param xlab xlab
#' @param ylab ylab 
#' @param lwd line width for the curves
#' @param fgcol col for focused line
#' @param bgcol col for background lines; can be a vector so that each line have different color
#' @param alpha alpha blending for background lines
#' @export
shadowPlot <- function(fitResList, indFocus=1, h=c(0.3, 0.5), main=NA, cex.main=1, type=c('plot', 'line'), tag='', 
	style='simple', xlab=NA, ylab=NA, ylim=NA, xlim=NA, lwd=2, lty=1, fgcol=NA, bgcol=80, alpha=NA){
	#browser()
	drugs <- sapply(fitResList, function(x) x$drug)
	if(length(unique(drugs))!=1) warning(sprintf('There are more than 1 drugs!'))
	if(is.na(main) & tag ==''){
		main <- str_c(unique(drugs), sep=', ')
	}
	indBg <- setdiff(1:length(fitResList), indFocus)
	if(is.na(alpha)){
		alpha <- guessAlpha_line(length(fitResList))
	}
	if(length(bgcol)==1) {
		# this allows to specify bgcol for each curve
		bgcol <- rep(bgcol, length(indBg))
	}
	if(is.na(fgcol)){
		# no focusing
		fgcol <- bgcol[1]
	}
	#browser() 
	# need to add plot to set up the figure
	plotOneExp(fitResList[[indFocus]], ind2plot='best', col=fgcol, type=c('plot', type), show='None', main=main, cex.main=cex.main, tag=tag, xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, lwd=lwd, lty=lty)
	# focused line: in case not plotted or color needs to be highlighted
	#plotOneExp(fitResList[[indFocus]], ind2plot='best', type=type, col=fgcol, show='None', main=main, cex.main=cex.main, tag=tag)
	#browser()
	j <- 1
	for(i in indBg) {
		#browser()
		plotOneExp(fitResList[[i]], ind2plot='best', type=setdiff(type, 'plot'), col=scales::alpha(bgcol[j], alpha), style=style, show='None', main=main, cex.main=cex.main, tag=tag)
		j <- j+1
	}
}
#shadowPlot(fitResList=fitList[ind], indFocus=1, fgcol='gray', alpha=0.3)

#shadowPlot(fitResList=fitList[1:10], indFocus=1)


#' guess alpha in a line plot
#' @param G number of genes or lines
#' @return a scalar between 0 and 1
guessAlpha_line <- function(G) {
	if(G<20) {
		return(0.9)
	} else if(G<100){
		return(0.5)
	} else if(G<200) {	
		return(0.3)
	} else if(G<500) {	
		return(0.2)
	} else if (G<1000) {	
		return(0.1)
	} else if (G<10000) {
		return(0.05)
	} else {
		return(0.01)
	}
}
