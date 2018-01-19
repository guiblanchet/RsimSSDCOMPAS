#' @title Simulation of a covariate on top of a deterninistic effect
#'
#' @description Simulation of a covariate W that has an effect on the species distribution on top of a deterninistic effect.
#'
#' @param SP Species matrix to be modified.
#' @param ENV Environmental variable matrix to which W is linked.
#' @param type Type of covariate effect, either 'LIN' for linear, 'FACT' for factor (can be used for the binary case) or 'NORM' for normal density.
#' \itemize{
#'	\item{LIN}{List with the following arguments in order : a vector containing the mean and sd of W, a vector giving the beta parameter of the regression : one for each species.}
#'	\item{FACT}{In order the number of repetition of each states, and their effect (for example, 3 states, rep 10 times each on 30 sites, and their respective effect 3,5,8 would give the vector : c(3,10,10.10,3,5,8)}
#'	\item{NORM}{List with the following arguments in order : a vector containing the mean and sd of W, a vector containing the maximum height of the effect of W : one for each species.}
#' }
#' @param corr Correlation matrix between E and W compeled by user. Must be square and symetric. Will change for FACT with the transformation of W.
#' @param parms Vector of parameters for the relationship between species and W. The number of elements and their meaning changes in function of the relationship between W and species.
#'
#' @details
#'
#' The relationship between E and W is simulated using the Cholesky decomposition of the correlation matrix and the scaled E matrix both specified by the user. First, we use the function solve to get the inverse of the Cholesky decomposition of the variance matrix of the combined matrix with E and W (a random standardize normal vector). We multiply the combined matrix E, W with the resulting matrix. This insures we get EXACTLY the correlation intended. To conclude the simulation, we multiply the previoulsy obtained matrix with the Cholesky decomposition of the correlation matrix. The modified W vector now has the requested correlation with the original E matrix. Once this is establish, the effect on the species is simulated in different ways depending on the user's specification of type.
#'
#' @return
#'
#' A list containing the new species data (SPW) and the resulting covariate (W).
#'
#' @author Marie-Hélène Ouellette
#'
#' @keywords datagen
#' @export

# Wrapper of SIMSSDR data to add linear or multinormal gradient effect #

SimPartial<-function(SP,ENV,type='LIN',parms,corr)
{
	# Initialisation
	SPW<-SP
	# Step 1 : scale E
	E_scale<-scale(ENV)
	cor.mat<-corr
	W_rand<-matrix(rnorm(length(E_scale)))
	EW_rand<-cbind(E_scale,W_rand)
	EW_solve<-EW_rand%*%solve(chol(var(EW_rand)))
	EW_new<-EW_solve %*% chol(cor.mat)


	W<-EW_new[,2]

	if(type=='LIN')
	{
		# Unscale W (preserving correlation)
		W<-W*parms[[1]][2] + parms[[1]][1]

		# Step 2 : Add W effect to species
		for(i in 1:ncol(SP))
		{
			# First simulate effect
			EffSPi_W<-parms[[2]][i]*W
			# Add effect to species in order of W # ADD THIS ORDER
			SPW[,i]<-SP[,i]+EffSPi_W
		}
	}
	if(type=='FACT')
	{
		# For'FACT', in order the number of repetition of each states, and their effect (for example, 3 states, rep 10 times each on 30 sites, and their respective effect 3,5,8 would give the vector : c(10,10.10,3,5,8))

		# Step 2 : modify the W vector into a factor
		W_ord<-order(W)
		W_fact<-rep(parms[[1]],times=parms[[2]])
		W<-W_fact
		# Effect in order of W
		EffSPi_W<-W_fact[W_ord]
		# Add effect to species
		SPW[,i]<-SP[,i]+EffSPi_W
	}
	if(type=='NORM')
	{
		# Unscale W (preserving correlation)
		W<-W*parms[[2]][1]+parms[[2]][2]
		# Step 2 : Add W effect to species
		for(i in 1:ncol(SP))
		{
			# First simulate effect (at height desired)
			EffSPi_W<-dnorm(scale(W))* parms[[3]][i]/dnorm(0)
			# Make distribution at maximum height
			# Add effect to species
			SPW[,i]<-SP[i]+EffSPi_W
		}
	}

	if(any(SPW<0)) warning('SPW contains negative value. Try different parameters for the W distribution.')
	res<-list(SPW,W)
	names(res)<-c('SPW','W')

	return(res)
}
