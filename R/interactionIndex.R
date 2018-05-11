if(FALSE){
# this is put into drexplorerExtra as the site for development; we can copy the source code and Rd file to add into drexplorer
#
## cd /home/ptong1/Backup/Package/

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
	source(file.path(osbaseZdrive, 'ptong1/Backup/Package/drexplorerExtra/R/base.R'))
    source(file.path(osbaseZdrive, 'ptong1/Backup/Package/drexplorerExtra/R/moonshot.R'))
    source(file.path(osbaseZdrive, 'ptong1/Backup/Package/drexplorerExtra/R/interactionIndex.R'))

}

#'  The drug combination data between SCH66336 and 4-HPR. This is an all-possible combination design.
#'
#' This dataset is published in Lee et al, 2007, Journal of Biopharmaceutical Statistics, Table 1
#'
#' \itemize{
#'   \item schd SCH66336 dose
#'   \item hpr 4-HPR dose
#'   \item resp response
#' }
#'
#' @format A data frame with 30 rows and 3 columns
#' @source \url{https://biostatistics.mdanderson.org/SoftwareDownload/}
#' @references Lee, J. Jack, et al. Interaction index and different 
#'    methods for determining drug interaction in combination therapy. Journal of biopharmaceutical statistics 17.3 (2007): 461-480.
#' @name nl22B2
NULL


#' The drug combination data of squamous cell carcinoma cells (UMSCC22B) surviving after 72 hours of treatment by single and
#' combination dose levels of SCH66336 and 4-HPR. This is a fixed-ratio design.
#'
#' This dataset is published in Lee et al, 2009, Statistics in biopharmaceutical research, Table 2
#'
#' \itemize{
#'   \item SCH66336 SCH66336 dose
#'   \item 4HPR 4-HPR dose
#'   \item resp response
#' }
#'
#' @format A data frame with 13 rows and 3 columns
#' @source \url{https://biostatistics.mdanderson.org/SoftwareDownload/}
#' @references Lee, J. J., & Kong, M. (2009). Confidence intervals 
#'   of interaction index for assessing multiple drug interaction. Statistics in biopharmaceutical research, 1(1), 4-17.
#' @name UMSCC22B
NULL








### downloaded from: https://biostatistics.mdanderson.org/SoftwareDownload/
##################################################################################################################
###                                                                                                            ###
###   Confidence Bound of Interaction Indices vs Effects at a fixed ray based on Lee and Kong (2006)           ###
###           Section 3 in this paper                                                                          ###
###                                                                                                            ###
###  INPUT:                                                                                                    ###
###  d1 and e1:     observed doses and effects for drug 1                                                      ###
###  d2 and e2:     observed doses and effects for drug 2                                                      ###
###  d12 and e12:   observed doses and effects for mixture at the fixed ratio d2/d1=d2.d1                      ###
###  E:             fixed effects, their corresponding interaction indices and confidence intervals are estimated.
###  alpha:         1-alpha is the size of the confidence intervals, alpha has the default value of 0.05.      ###
###                                                                                                            ###
###  OUTPUT:                                                                                                   ###
###  ii:            the estimated interaction indices corresponding to the input effects E                     ###
###  ii.low:        the estimated lower confidence intervals for ii                                            ###                        
###  ii.up:         the estimated upper confidence intervals for ii                                            ###                        
###                                                                                                            ###
##################################################################################################################  
#
# originally called: CI.delta from CI_IIV2.SSC                          
#
 
#' Estimate the interaction index as well as its confidence interval using delta method for fixed ratio design
#'
#' this code is extracted from the source code distributed at https://biostatistics.mdanderson.org/SoftwareDownload/
#' Two versions are included in the source code: CI_IIV2 2008.SSC and CI_IIV2.SSC
#' At first we implement CI_IIV2.SSC; this gives wider CI band; so we decide to use  CI_IIV2 2008.SSC. Further, CI_IIV2 2008.SSC has
#' a more recent date and thus most updated.
#'
#' @param d1 dose for drug 1 where drug 2 has dose equal to 0                                         
#' @param e1 corresponding scaled response for drug 1, must between 0 and 1                                          
#' @param d2 dose for drug 2 where drug 1 has dose equal to 0                                               
#' @param e2 corresponding scaled response for drug 2, must between 0 and 1                                
#' @param d12 the sum of doses from drug 1 and drug 2 when both drugs are administered                                          
#' @param e12 corresponding scaled response when both drugs are administered                                          
#' @param d2.d1 the ratio of the ray design. This is the fixed ratio of dose 2 divided by dose 1
#' @param E a vector of responses (between 0 and 1) where IAI and confidence interval are to be computed from.                                           
#' @param alpha significance level of confidence interval             
#' @return a data frame with columns IAI, IAI.low, IAI.up, E, dx1 (corresponding dose of drug 1), dx2 (corresponding dose of drug 1), 
#'  dx12 (corresponding dose of combined drug, same as definition of d12)
CI.delta <- function(d1, e1, d2, e2, d12, e12, d2.d1, E, alpha=0.05)
{
     lm1 <- lm(log(e1/(1-e1))~log(d1))
     dm1 <- exp(-summary(lm1)$coef[1,1]/summary(lm1)$coef[2,1])
     lm2 <- lm(log(e2/(1-e2))~log(d2))
     dm2 <- exp(-summary(lm2)$coef[1,1]/summary(lm2)$coef[2,1])
     lmcomb <- lm(log(e12/(1-e12))~log (d12))
     dm12 <- exp(-summary(lmcomb)$coef[1,1]/summary(lmcomb)$coef[2,1]) 
     Dx1 <- dm1*(E/(1-E))^(1/summary(lm1)$coef[2,1])
     Dx2 <- dm2*(E/(1-E))^(1/summary(lm2)$coef[2,1])
     dx12 <- dm12*(E/(1-E))^(1/summary(lmcomb)$coef[2,1])
     iix <- (dx12/(1+d2.d1))/Dx1+(dx12*d2.d1/(1+d2.d1))/Dx2
     lm1.s <-summary(lm1)
     lm2.s <-summary(lm2)
     lm12.s <-summary(lmcomb)
     c1 <- 1.0/lm1.s$coef[2,1]^2*lm1.s$coef[1,2]^2
     temp <- - mean(log(d1))*lm1.s$coef[2,2]^2
    ### temp <- lm1.s$coef[1,2]*lm1.s$coef[2,2]*lm1.s$cor[1,2]   ### covariance of b0 and b1
     c1 <- c1+2.0*(log(E/(1-E))-lm1.s$coef[1,1])/lm1.s$coef[2,1]^3*temp
     c1 <- c1+(log(E/(1-E))-lm1.s$coef[1,1])^2/lm1.s$coef[2,1]^4*lm1.s$coef[2,2]^2
     c2 <- 1.0/lm2.s$coef[2,1]^2*lm2.s$coef[1,2]^2
     temp <- - mean(log(d2))*lm2.s$coef[2,2]^2
    ### temp <- lm2.s$coef[1,2]*lm2.s$coef[2,2]*lm2.s$cor[1,2]   ### covariance of b0 and b1
     c2 <- c2+2.0*(log(E/(1-E))-lm2.s$coef[1,1])/lm2.s$coef[2,1]^3*temp
     c2 <- c2+(log(E/(1-E))-lm2.s$coef[1,1])^2/lm2.s$coef[2,1]^4*lm2.s$coef[2,2]^2
     c12 <- 1.0/lm12.s$coef[2,1]^2*lm12.s$coef[1,2]^2
     temp <- - mean(log(d12))*lm12.s$coef[2,2]^2
   ### temp <- lm12.s$coef[1,2]*lm12.s$coef[2,2]*lm12.s$cor[1,2]   ### covariance of b0 and b1
     c12 <- c12+2.0*(log(E/(1-E))-lm12.s$coef[1,1])/lm12.s$coef[2,1]^3*temp
     c12 <- c12+(log(E/(1-E))-lm12.s$coef[1,1])^2/lm12.s$coef[2,1]^4*lm12.s$coef[2,2]^2
     var.ii <-((dx12/Dx1)^2*c1+(dx12*d2.d1/Dx2)^2*c2+(1.0/Dx1+d2.d1/Dx2)^2*dx12^2*c12)/(1+d2.d1)^2 
     t975 <- qt(1-alpha/2,length(d1)+length(d2)+length(d12)-6)
     iix.low1 <- iix*exp(-t975*var.ii^0.5/iix)
     iix.up1 <- iix*exp(t975*var.ii^0.5/iix)
	 #browser()
     #return(list(ii=iix, ii.low=iix.low1, ii.up=iix.up1))
     # add the corresponding dose for visualization
	 IAI <- data.frame(IAI=iix, IAI.low=iix.low1, IAI.up=iix.up1, E=E, dx1=Dx1, dx2=Dx2, dx12=dx12)
	 #browser()
	 return(IAI)
}



if(FALSE){
	# example
	#par(mfrow=c(1,2),mai=c(0.4,0.4,0.4,0.0), mgp=c(1.0, 0.2, 0))
# data published in Table 2 of: Confidence Intervals of Interaction Index for Assessing Multiple Drug Interaction
dose1 <-c(0.1,0.5,1,2,4,0,0,0,0, 0.1, 0.5, 1, 2)
dose2 <-c(0,0,0,0,0,0.1,0.5,1,2, 0.1, 0.5, 1, 2)
fa <-c(0.6701,0.6289,0.5577, 0.455,0.3755,0.7666, 0.5833,0.5706,0.4934,0.6539,0.4919,0.3551,0.2341)
UMSCC22B <- data.frame(dose1=dose1, dose2=dose2, e=fa)
#save(UMSCC22B, file='UMSCC22B.RData')
d2.d1 <-1
#####  median effect plots (Figure 3, Panel A )               #########
median.effect.plots(dose1, dose2, fa, name1="SCH66336", name2="4HPR", d2.d1=1/1)
####  Construct confidence intervals and confidence bounds  (Figure 3, Panel B )  ###
ind1 <- dose1!=0 & dose2==0
ind2 <- dose1==0 & dose2!=0
ind12 <- dose1!=0 & abs(dose2-d2.d1*dose1)<0.00001
E<- seq(0.05, 0.95, 0.005)
Er <-c(E,fa[ind12])
#CId.out <- fitIAI_ray(d1=dose1[ind1], e1=fa[ind1], d2=dose2[ind2], e2=fa[ind2], d12=dose1[ind12]+dose2[ind12], e12=fa[ind12], d2.d1, E, alpha=0.05)

    osbaseZdrive <- ifelse(.Platform$OS.type=="windows", "//mymdafiles/usersdqs1", "/home") 
	source(file.path(osbaseZdrive, 'ptong1/Backup/Package/drexplorerExtra/R/base.R'))
    source(file.path(osbaseZdrive, 'ptong1/Backup/Package/drexplorerExtra/R/moonshot.R'))
    source(file.path(osbaseZdrive, 'ptong1/Backup/Package/drexplorerExtra/R/interactionIndex.R'))

library(drexplorerExtra)
data(nl22B2)	
data(UMSCC22B)	
fit_fixedRay <- fitIAI(d1=UMSCC22B$dose1, d2=UMSCC22B$dose2, e=UMSCC22B$e, name1='SCH66336', name2='4HPR')
plotIAI(fit_fixedRay, type='IAI', mode='both') # this reproduce the example but not the paper at first; then we shift to CI_IIV2 2008.SSC and now it works!
plotIAI(fit_fixedRay, type='medianEffect', mode='both') # this reproduce the CI paper
plotIAI(fit_fixedRay, type='doseResponseCurve', mode='both')

#### this reproduces the example figures exactly in the source code: SYNERGY_V3_original
fit_allPoss <- fitIAI(d1=nl22B2$schd, d2=nl22B2$hpr, e=nl22B2$y1, name1='SCH66336', name2='4HPR')
#medianEffect <- fit_median_efect(d1=nl22B2$schd, d2=nl22B2$hpr, e=nl22B2$y1, d2.d1=res_design$d2.d1)
res_design <- detect_ray_design(d1=nl22B2$schd, d2=nl22B2$hpr, e=nl22B2$y1)
plotIAI(fit_allPoss, type='IAI', mode='both') # this use the CI method; approximately
plotIAI(fit_allPoss, type='contour', mode='both') # this is slightly different from the paper, due to the author's inconsistency
plotIAI(fit_allPoss, type='medianEffect', mode='both')
plotIAI(fit_allPoss, type='doseResponseCurve', mode='both')

fit <- fitIAI(d1=dose1, d2=dose2, e=fa)
plotIAI(fit, type='IAI', mode='both')
x11()
plotIAI(fit, type='IAI', mode='dose')
x11()
plotIAI(fit, type='IAI', mode='response')
plotIAI(fit, type='medianEffect')
plotIAI(fit, type='doseResponseCurve')
plotIAI(fit, type='contour')
}