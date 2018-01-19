#' @title Evaluating generalized Beta distribution for use in the Rcompas
#'
#' @description Evaluates values of beta distribution with parameters provided by the user.
#'
#' @param A0 A value.
#' @param M A vector of m values, as many as their is gradients.
#' @param R A vector of r values, as many as their is gradients.
#' @param ALPHA A vector of 'alpha' values, as many as their is gradients.
#' @param GAMMA A vector of 'gamma' values, as many as their is gradients.
#' @param X A matrix of gradients.
#'
#' @return
#'
#' A vector of values
#'
#' @author Marie-Hélène Ouellette
#'
#' @keywords datagen
BETA<-function(A0,M,R,ALPHA,GAMMA,X)
{
	b<-abs(ALPHA)/(abs(ALPHA)+abs(GAMMA))
	d<-sum((b^abs(ALPHA))*(1-b)^abs(GAMMA))
	# For all gradients
	for(i in 1:length(M))
	{
		A<-(A0/d)*(((X[,i]-M[i])/R[i]+b[i])^abs(ALPHA[i]))*(1-((X[,i]-M[i])/R[i]+b[i]))^abs(GAMMA[i])
	}
	return(A)
}
