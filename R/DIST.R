#' @title Generating a value from a specified distribution for use in the Rcompas
#'
#' @description From which distribution should we take numbers (name), and how many (n)?
#'
#' @param name One of 'unifrandom', 'lognormal', 'lograndom', 'exponential', 'gamma' and 'Poisson' distributions.
#' @param n The number of values to be generated.
#' @param \dots other arguments to be passed to generating functions.
#'
#'
#' @return
#'
#' A vector of \code{n} values generated from the distribution specified by the user.
#'
#' @author Marie-Hélène Ouellette
#'
#' @keywords datagen
DIST<-function(name,n,...)
{
	# From which distribution should we take numbers (name), and how many (n).
	# Implemeted are :

	# 'unifrandom'
	if(name=='unifrandom') num<-runif(n,...)
	#'normal'
	if(name=='normal') num<-rnorm(n,...)
	#'lognormal'
	if(name=='lognormal') num<-rlnorm(n,...)
	#'lograndom'
	if(name=='lograndom') num<-rlog(n,...)
	#'exponential'
	if(name=='exponential') num<-rexp(n,...)
	#'gamma'
	if(name=='gamma') num<-rgamma(n,...)
	#'Poisson'
	if(name=='Poisson') num<-rpois(n,...)

	return(num)
}
