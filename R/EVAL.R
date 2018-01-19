#' @title Build a matrix of parameters
#'
#' @description Build a matrix of parameters from the list BETAA for the main function Rcompas
#'
#' @param BETAAjx Piece of a BETAA list that contains a number of parameters that will be used later on to evaluate species abundances from the Rcompas function.
#' @param MATi Matrix of parameters that will be modified, build by this functions.
#' @param i Number corresponding the gradient.
#' @param j Number corresponding to the profile.
#' @param LINE Number corresponding to the species, the line number of Mati to fill.
#' @param col The column of Mati to fill.
#'
#'
#' @return
#'
#' The matrix Mati filled
#'
#' @author Marie-Hélène Ouellette
#'
#' @keywords datagen
EVAL<-function(BETAAjx,MATi,i,j,LINE,col)
{
	if(BETAAjx[[1]][i]=='value') # It is set
	{
		# Than for all species, values are set and are in the second element of the list
		MATi[LINE:(LINE+length(unlist(BETAAjx[[2]][i]))-1),col]<-unlist(BETAAjx[[2]][i])

	}

	if(BETAAjx[[1]][i]=='random') # It is not set
	{

		# Evaluate function DIST at parameters specified in the list
		MATi[LINE:(LINE+as.numeric(substr(as.character(BETAAjx[[3]][i]),3,3))-1),col]<-eval(parse(text=noquote(paste('DIST("',BETAAjx[[2]][i],'",',BETAAjx[[3]][i],')', sep=''))))
	}

	return(MATi)
}
