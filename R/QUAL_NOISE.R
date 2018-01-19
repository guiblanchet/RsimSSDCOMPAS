#' @title Add qualitative noise to a species simulated in Rcompas
#'
#' @description Add qualitative noise to a species simulated in Rcompas, based on its beta parameters.
#'
#' @param SPi Piece of a BETAA list that contains a number of parameters that will be used later on to evaluate species abundances from the Rcompas function.
#' @param Mi m parameter for the specie.
#' @param Ri r parameter for the specie.
#' @param ALPHAi alpha parameter for the specie.
#' @param GAMMAi gamma parameter for the specie.
#' @param X the gradients on which the specie are simulated.
#'
#' @return
#'
#' the abundances of a species
#'
#' @author Marie-Hélène Ouellette
#'
#' @keywords datagen
QUAL_NOISE<-function(SPi,Mi,Ri,ALPHAi,GAMMAi,X)
{
	P0_EVALi<-runif(1, min=0, max=1)
	# Qualitative response function
	PROBSi<-BETA(P0_EVALi,Mi,Ri,ALPHAi,GAMMAi,X)
	UNIF_P<-runif(nrow(X), min=0, max=1)
	# If the random value is higher than the probability value, the species does not appear in the sample.

	if(any(UNIF_P>PROBSi,na.rm = TRUE)) SPi[which(UNIF_P>PROBSi)]=0
	return(SPi)
}
