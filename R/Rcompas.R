#' @title R version of simulation programs COMPAS and JCOMPAS
#'
#' @description Generates species data that follows specific concepts and
#' hypotheses about the properties of the community patterns. Particularly,
#' concepts of asymetric physiological response (generalized beta distribution),
#' species guilds, interspecific interactions, carrying capacity, qualitative
#' noise (presence or absence of species at a particular site) and quantitative
#' noise (scatter around the fitted response) are integrated in the simulation
#' framework. It follows very closely the code from COMPAS Minchin(1987) and
#' JCOMPAS from Miquel de Caceres, January 2003. The main difference between
#' Rcompas and its two sister packages is that the explanatory variables must be
#' provided. This permits implemetation into \code{SimSSDR} for example, in
#' which autocorrelation can then be added to the species data. There is no
#' sampling design implemented.
#'
#' @param BETAA A list. For each profile PX (each profile corresponds to a number of species chosen by the user), list of objects defining the parameter of the beta distribution for the species on the gradients. The names of the list objects contained in this list are PX$A,PX$m,PX$r,PX$alpha and PX$gamma. See the Details section for more precision on the usage.
#' @param X Matrix of environmental gradient(s).
#' @param G Guilds or groups. G is a logical vector: \code{TRUE} if species are prganized in guilds.
#' @param MAJOR The percentage of major species (species with higher modal abundances) which will be regularly spread out on the gradient. If \code{FALSE}, are scattered as other 'non major' species. Maximum area under the curve at which a species is not considered major (Gauch et Whittaker 1972).}
#' @param HISI Include or omit (0) the hypothesis of Inter-specific interaction. If included, HISI is a list with the following levels:
#' \itemize{
#'	 \item \code{n} Number of interactions per species permitted (maximum, could be less because of subsequent rules). 999 means no limit.
#'	 \item \code{g} Logical. Species must be from the same guild to interact.
#'	 \item \code{mc} Must have similar modal coordinates, the proximity given by the user in percentage of the range of the gradient. If 999, than this rule does not apply.
#'	 \item \code{mv} Must have similar modal value, the disimilarity given by the minimum gower distance. If 999, than this rule does not apply.
#' }
#' @param HQUAL Logical. Include or omit the hypothesis of qualitative noise, aka appearance probability. Applies to all species.
#' @param HQUANT Logical.
#' @param HCCA Logical. Include or omit ($type=0) the carrying capacity adjustment hypothesis. Other than 0 (no limit in carying capacity), this argument can take 3 other values : ($type=1) the stand abundance is constant over all the gradient (must be given in argument $parm) ($type=2) it is in a linear relationship with the gradients (intercept and slopes (in that order) must be provided in $parm as a vector) ($type=3) the carrying capacity is defined by a beta function with specified parameters ($A,$m,$r,$alpha,$gamma; provided in $parm). For this third option, the dipersion parameter of the distribution (if it exists), can be set to be fixed or proportional to a simple function of A by the user.Either way, the species abundances are rescaled accordingly.
#' @param HSYM Logical. Hypothesis of symmetry of the distribution for each profile. Provided by a vector of \code{TRUE} (symmetry) or \code{FALSE} (non symmetry, thus alpha and gamma are different) with one value per profile. If \code{TRUE}, alpha is simulated are the user has asked, and gamma is set to be equal to alpha.
#' @param SAR Logical. Add autocorrelation to species data or not.
#' @param \dots Further arguments for the autocorrelation application.
#'
#' @details
#'
#' As stated above, BETAA is a list of object describing the profile PX, where each profile corresponds to a number of species chosen by the user. Each profile is in turn a list of objects defining the parameter of the beta distribution for the species on the gradients. The names of the list objects contained in this list are PX$A,PX$m,PX$r,PX$alpha and PX$gamma. Two types of structure can be used for each of them, depending on if the user wants to specify the parameters himself or want them to be picked at random from a specific distribution. To simplify the following instructions, lets suppose that the object we're working on is 'A' in the list 'PARM' for a certain profile 'P' :
#'
#' \itemize{
#'   \item For the first type, the value(s) is (are) set by the user. The first element of the list object states that the values are readily given by the user by stating the character string 'value'. The second and last element is the value(s) coerced in a vector. For example, we would have BETAA$P1$A<-list('value',c(1,2,3)) for 3 species having the A values 1, 2 and 3 respectivlely on one gradient. The extention to more than one gradient is given in the form of a list, see \code{examples}.
#'   \item For the second type, the values are chosen to be picked at random out of a specific distribution. The first element of the list object is then the character string 'random'. The user has to state two more elements to specify exactly from which distribution the parameters should be simulted. Thus the second element of the list is the name of the distribution (one of 'unifrandom', 'normal', 'lognormal', 'lograndom', 'exponential', 'gamma' or 'Poisson'). The third element contains the extra arguments needed to specify the parameters of the distribution, coerced in a vector. For example, BETAA$P1$A<-list('random','normal',c(n=2,mean=1,sd=0)). Note that the argument 'n' gives the number of species. See the 'BETA' function for more details on those (and other) parameters.
#' }
#' The user may define as many profiles as needed : one just has to provide them in the form of a list P1 - PX... . The number of species is stated in the arguments : for set values it's in the second element (ex : number of elements in the vector PARM$P$A[2]), and for the simulated ones, it's in the third (ex : n in PARM$P$A[3] = c(n=3,mean=0,sd=1)).
#'
#' @note
#'
#' Guilds in Rcompas have certain particularities. Guilds can be simulated by expressively setting A0 the m parameters in the BETAA list. Also, the argument HSYM should be set to \code{FALSE} so major species in the guilds STAY in the guilds, and are not relocated on the gradient.
#'
#' @return
#'
#' A matrix giving the simulated species abundances
#'
#' @author Marie-Hélène Ouellette
#'
#' @references
#'Minchin, P. R. 1987. Simulation of multidimensional community patterns: towards a comprehensive model. \emph{Plant Ecology} \strong{71}: 145-156.
#'
#'de Cáceres, M. 2003. JCOMPAS. Dept. Biologia Vegetal, Universitat de Barcelona: JCOMPAS 1.0 user's manual.
#'
#' @examples
#'
#' #---------------------------------------------#
#' ### Examples of how to build the BETAA list ###
#' #---------------------------------------------#
#'
#' # *** In the case where we have many gradients, we need to specify the parameters for each gradients. This is done by extending the elements stated above into lists. Here are three examples. *** #
#'
#' #-----------#
#' ### Example 1
#' ### Values specified by user for three gradients (2 species, 3 gradients)
#' #-----------#
#' BETAA<-list()
#' BETAA$P1$A <- list(list('value'), list(c(2,3)),list(n=2)) # On all 3 gradients, A0 is the same for the species.
#' BETAA$P1$m <- list(list('value','value','value'), list(c(2,4),c(3,3),c(6,4),list(n=2))) # On each gradient, the values are set. Again, for each gradient, we give the values for both species.
#'
#' #-----------#
#' ### Example 2
#' ### Values picked at random for three gradients (2 species, 3 gradients)
#' #-----------#
#' BETAA<-list()
#' BETAA$P1$A <- list(list('random'), list('normal'), list(list(n=2,mean=1,sd=0))) # On all 3 gradients, A0 are the same for the species.
#' BETAA$P1$m <- list(list('random','random','random'), list('normal','lognormal','lognormal'), list(c(n=2,mean=1,sd=0),c(n=2, meanlog = 0, sdlog = 1),c(n=2, meanlog = 0, sdlog = 3))) # Two values are picked at random on the three gradients at random form different distributions. The distributions are the same for both species.
#'
#' #-----------#
#' ### Example 3
#' ### A mixture of values set and picked at random from different distributions (2 species, 3 gradients)
#' #-----------#
#' BETAA<-list()
#' BETAA$P1$A <- list(list('random'), list('lognormal'), list(c(n=2, meanlog = 0, sdlog = 1)))
#' BETAA$P1$m <- list(list('value','random','random'), list(c(2,3),'lognormal','lognormal'), list(n=2,c(n=2, meanlog = 0, sdlog = 1)),list(n=2,c(n=2, meanlog = 1, sdlog = 3)))
#'
#' #-----------#
#' ### Example 4
#' ### A complete example of how the list should be build
#' #-----------#
#'	BETAA<-list() # 2 gradients
		# P1 has 3 species
#'	BETAA$P1$A<-list(list('value'), list(c(20,30,50)),list(c('n=3')))
#'	BETAA$P1$m<-list(list('random','random'), list('normal','lognormal'), list(c('n=3,mean=50,sd=10'),c('n=3, meanlog = 0, sdlog = 1')))
#'	BETAA$P1$r<-list(list('random','random'), list('normal','lognormal'), list(c('n=3,mean=20,sd=3'),c('n=3, meanlog = 4, sdlog = 1')))
#'	BETAA$P1$alpha<-list(list('value','random'), list(c(2,3,4),'lognormal'), list(c('n=3'),c('n=3, meanlog = 0, sdlog = 1')))
#'	BETAA$P1$gamma<-list(list('random','random'), list('normal','lognormal'), list(c('n=3,mean=0,sd=1'),c('n=3, meanlog = 0, sdlog = 1')))
		# P2 has 2 species
#'	BETAA$P2$A <- list(list('random'), list('normal'), list(c('n=2, mean=20, sd=1')))
#'	BETAA$P2$m <- list(list('random','random'), list('normal','lognormal'), list(c('n=2,mean=30,sd=1'),c('n=2, meanlog = -1, sdlog = 1')))
#'	BETAA$P2$r <- list(list('random','random'), list('normal','lognormal'), list(c('n=2,mean=20,sd=1'),c('n=2, meanlog = 2, sdlog = 1')))
#'	BETAA$P2$alpha <- list(list('random','random'), list('normal','lognormal'), list(c('n=2,mean=0,sd=1'),c('n=2, meanlog = 0, sdlog = 1')))
#'	BETAA$P2$gamma <- list(list('random','value'), list('normal',c(4,7)), list(c('n=2,mean=0,sd=1'),c('n=2')))
#'
#' #-----------#
#' ### Example 5
#' ### Full examples of the main function usage using the last BETAA constructed
#' #-----------#
#'
#' X<-cbind(c(1:100),c(rnorm(100)))
#' G<-c(TRUE,TRUE) # The two profiles are two guilds
#' MAJOR<-FALSE
#' HISI<-list(n=2,g=TRUE,mc=999,mv=999)
#' HQUAL<-TRUE
#' HQUANT<-TRUE
#' HCCA<-list(type=1,parm=10)
#' HSYM<-c(TRUE,TRUE) # For the two profiles
#' SAR=FALSE
#'
#' Rcompas(BETAA=BETAA,X=X,G=G,MAJOR=MAJOR,HISI=HISI,HQUAL=HQUAL,HQUANT=HQUANT,HCCA=HCCA,HSYM=HSYM,SAR=SAR)
#'
#' @keywords datagen
#' @export
Rcompas<-function(BETAA,X,G,MAJOR=FALSE,HISI,HQUAL,HQUANT,HCCA,HSYM,SAR,...)
{
	require(VGAM)
	require(FD)


### Main data simulation ###
# Here we sought to construct the response parameter matrix. For each profile, simulate the species's parameters on all gradients. We simulate one matrix per gradient. Assembling comes later, when the abundances are simulated.

# -- Simulation of Physiological Response -- #

# All parameters (A,m,r,alpha and gamma) are stoked in a list PARM, abundances will be estimated later. EAch element of PARM is a matrix that contains the parameters per gradient.

	PARM<-list()

	if(length(BETAA)>0) # length(BETAA) gives the number of profiles
	{
		# Get total number of species N
		N=0
		for(j in 1:length(BETAA))
		{
			N=N+as.numeric(substr(as.character(BETAA[[j]]$A[[3]]),3,3))
		}

		A_EVAL=0
			for(j in 1:length(BETAA))
			{
				# A : the same for all gradients.
				# Is it a set value or a value picked at random
				if(BETAA[[j]]$A[[1]]=='value') # It is set
				{
					# Than for all species, values are set and are in the second element of the list
					A_EVALi<-unlist(BETAA[[j]]$A[[2]])
				}

				if(BETAA[[j]]$A[[1]]=='random') # It is not set
				{
					# Evaluate function DIST at parameters specified in the list

					A_EVALi<-eval(parse(text=noquote(paste('DIST("',BETAA[[j]]$A[[2]],'",',BETAA[[j]]$A[[3]],')', sep=''))))
				}
				A_EVAL<-c(A_EVAL,A_EVALi)
			}
			A_EVAL<-A_EVAL[-1]

		for(i in 1:ncol(X)) # One matrix per gradient
		{
			MATi<-mat.or.vec(N,5)
			LINE<-1

			# A
			MATi[1:N,1]<-A_EVAL
			NSP<-mat.or.vec(1,length(BETAA))

			for(j in 1:length(BETAA)) # For each profile (with a certain number of species associated), get the parameters, and place them in a matrix PARM. The column are in order : A, m, r, alpha, gamma. Individual species are lines.
			{

				NSP[,j]<-as.numeric(substr(as.character(BETAA[[j]]$A[[3]]),3,3))
				# m

				MATi<-EVAL(BETAA[[j]]$m,MATi,i,j,LINE,2)

				# r

				MATi<-EVAL(BETAA[[j]]$r,MATi,i,j,LINE,3)

				# alpha

				MATi<-EVAL(BETAA[[j]]$alpha,MATi,i,j,LINE,4)

				# gamma

				if(HSYM[j]) MATi[LINE:(LINE + NSP[,j]-1),5]<-MATi[LINE:(LINE + NSP[,j]-1),4]
				if(!HSYM[j]) MATi<-EVAL(BETAA[[j]]$gamma,MATi,i,j,LINE,5)

				# Number of species per profile. Usefull for the section 'Modification of the community data' (guild).


				LINE<-LINE + NSP[,j]


			}
			PARM[[i]]<-MATi
		}
	}

	# If 'major' species option is TRUE: regularly distribute the major species, simulating the other parameters accordingly. For the others, simulated other parameters of the distribution according to the user's demands. #
	if(MAJOR!=FALSE)
	{
		# The percentage of species to be considered as 'major' is MAJOR. Order all values of A0, adjust m accordingly on the 'major' species.

		# Species to be spaced evenly : major species
		MAJOR_SP<-order(PARM[[1]][,1],decreasing = TRUE)[1:round(MAJOR/100*length(PARM[[1]][,1]))]

		# Change their m values on all gradients : evenly spaced
		for(i in 1:ncol(X)) # One matrix per gradient
		{
			SPACEi<-(range(X[,i])[2]-range(X[,i])[1])/(length(MAJOR_SP)+1)
			Mi<-SPACEi*c(1:(length(MAJOR_SP)))

			# Change m in the matrix
			PARM[[i]][,2][MAJOR_SP]<-Mi
		}
	}


	### SIMULATING SPECIES DATA ###
	# For each species, we simulate their abundances as stated by the parameters here above chosen.
	SP<-as.matrix(mat.or.vec(nrow(X),nrow(PARM[[1]])))
	A0<-as.matrix(mat.or.vec(nrow(PARM[[1]]),1))
	M<-as.matrix(mat.or.vec(nrow(PARM[[1]]),ncol(X)))
	R<-as.matrix(mat.or.vec(nrow(PARM[[1]]),ncol(X)))
	ALPHA<-as.matrix(mat.or.vec(nrow(PARM[[1]]),ncol(X)))
	GAMMA<-as.matrix(mat.or.vec(nrow(PARM[[1]]),ncol(X)))

	for(i in 1:nrow(PARM[[1]]))
	{


		A0[i]<-PARM[[1]][,1][i]
		for(j in 1:ncol(X)) # gather parameters for all gradients
		{

			M[i,j]<-PARM[[j]][i,2]
			R[i,j]<-PARM[[j]][i,3]
			ALPHA[i,j]<-PARM[[j]][i,4]
			GAMMA[i,j]<-PARM[[j]][i,5]
		}
		SP[,i]<-BETA(PARM[[1]][,1][i],M[i,],R[i,],ALPHA[i,],GAMMA[i,],X)
	}

	# Set NAN to 0
	SP[is.na(SP)]=0


	# -- Modification of the community data -- 4 hypothesis + 1 : autocorrelation #

	#0) autocorrelation
	if(SAR)
	{
		for(j in 1:ncol(SP))
		{
			resultj<-vector(mode = "numeric", length = nx*ny)
			SAR_vecj<-.Fortran("inpsgsm",as.integer(nx*ny),as.integer(nx),as.integer(ny),as.numeric(nug2),as.double(range21),as.double(range22),as.double(resultj),PACKAGE="RsimSSDCOMPAS")
			names(SAR_vecj)<-c('nxymax','nx','ny','nugget','range1','range2','result')

			SP[,j]=as.vector(SP[,j])+as.vector(SAR_vecj$result)
		}
	}



	#1)	Inter-specific interaction
	if(all(HISI!=0))
	{
		# We must build a list of candidate for interactions for each species. This list will be symetric, in the sens that if SPj can interact with SPi, than SPi can interact with SPj.

		# GG, MC and MV are square matrices that define the candidate for interactions for each species. By multiplying them (Adamart multiplication), we will get the possible candidate respecting all rules.
		GG<-mat.or.vec(ncol(SP),ncol(SP))
		MC<-mat.or.vec(ncol(SP),ncol(SP))
		MV<-mat.or.vec(ncol(SP),ncol(SP))

		# RULE 1 : g=TRUE or g=FALSE
		if(!HISI$g) GG<-matrix(1,ncol(SP),ncol(SP))
		if(HISI$g)
		{
			# G : TRUE or FALSE for guilds for each profile
			# NSP : Number of species per profile
			GUILD_TF<-rep(G,NSP)
			PROFILE_NUM<-rep(c(1:length(NSP)),NSP)

			# Species belonging to a group are numbered according to their affiliation. Others are 0.
			GUILD_NUM<-GUILD_TF*PROFILE_NUM

			WHICH_SP_GUILD<-which(GUILD_NUM!=0)

			D<-dist(GUILD_NUM,method='maximum',diag = TRUE, upper = TRUE)
			D_MAT<-as.matrix(D)
			S_MAT<-D_MAT+1
			S_MAT[which(S_MAT>1,arr.ind = TRUE)]<-0

			GG[,WHICH_SP_GUILD]<-S_MAT[,WHICH_SP_GUILD]
		}

		# RULE 2 : mc=dist_max or mc=999
		if(HISI$mc==999) MC<-matrix(1,ncol(SP),ncol(SP))
		if(HISI$mc!=999)
		{
			# The first step is to calculate the proximity of modal coordinates of species and establish if we pass a certain threshold set by the user (HISI$mc). We use the Gower distance for quantitative variables (gowdis of the FD library). The range of each variable is taken as : [MIN(min(X),min(m)):MAX(max(X),max(m))].

			MIN<-apply(rbind(X,M),2,min)
			MAX<-apply(rbind(X,M),2,max)
			RANGES<-MAX-MIN
			DISTA<-gowdis(M,Rj=RANGES)
			SIM_MAT<-1-as.matrix(DISTA)
			SIM_MAT[which(SIM_MAT>HISI$mc/100,arr.ind = TRUE)]=1
			SIM_MAT[which(SIM_MAT<=HISI$mc/100,arr.ind = TRUE)]=0
			MC<-SIM_MAT
		}

		# RULE 3 : mv=gower dist_max or mv=999
		if(HISI$mv==999) MV<-matrix(1,ncol(SP),ncol(SP))
		if(HISI$mv!=999)
		{
			DISTA<-gowdis(as.matrix(A0))
			SIM_MAT<-1-as.matrix(DISTA)
			SIM_MAT[which(SIM_MAT>HISI$mc/100,arr.ind = TRUE)]=1
			SIM_MAT[which(SIM_MAT<=HISI$mc/100,arr.ind = TRUE)]=0
			MV<-SIM_MAT
		}
		# Product of all rules minus the maximum number of interactions
		ISI<-GG*MC*MV
		diag(ISI)=0

		# Apply last rule : maximum number of interactions
		# RULE 4 : n=# or n=999
		if(HISI$n==999) C1<-ISI
		if(HISI$n!=999)
		{
			# For all columns of ISI (same as rows as the matrix is symetric), count the ones and if we are over, pick 'n' at random to stay 1, others go to 0.
			n<-HISI$n
			if(any(colSums(ISI)>n))
			{
				# Columns (species) in which there is more than n interactions
				ID_COL_N<-which(colSums(ISI)>n)
				# Coordinates to where the ones are
				WHICH_ISI_1<-which(ISI==1,arr.ind=TRUE)
				for(i in 1:length(ID_COL_N)) # For each species where there is more than the maximum number of interactions
				{
					# Pick, at random, amongst the interactions for this species, n interactions for that species (i).
					N_inter<-sample(which(ISI[,ID_COL_N[i]]==1, arr.ind=TRUE),n)
					# Others are set to 0
					NON_INTER<-which(ISI[,ID_COL_N[i]]==1)[-match(N_inter,which(ISI[,ID_COL_N[i]]==1))]
					ISI[NON_INTER,ID_COL_N[i]]<-0
				}
			}
		}
	}

	#2)	Qualitative noise : simulated from the generalized BETA distribution, all parameters are the same as with the species simulation, except for one : A0 is replaced by P0, the probability of occurence which is chosen at random from a random dsitribution with values between 0 and 1.
	if(HQUAL)
	{
		# Simulate a probability of occurency for each species
		for(i in 1:sum(NSP))
		{
			SP[,i]<-QUAL_NOISE(SP[,i],M[i,],R[i,],ALPHA[i,],GAMMA[i,],X)
		}
	}

	#3)	Carrying capacity adjustments
	if(HCCA$type!=0)
	{
		if(HCCA$type==1)
		{
			CCA<-matrix(HCCA$parm,1,nrow(X))
		}

		if(HCCA$type==2)
		{
			CCA<-HCCA$parm[1]+HCCA$parm[-1]*X
		}

		if(HCCA$type==3)
		{
			CCA<-BETA(HCCA$parm$A,HCCA$parm$m,HCCA$parm$r,HCCA$parm$alpha,HCCA$parm$gamma,X)
		}

		# Ponderate lines by CCA

		POND_LINES<-which(rowSums(SP)>CCA)
		SP[POND_LINES,]<-SP[POND_LINES,]/rowSums(SP[POND_LINES,])*CCA[POND_LINES]
	}

	#4)	Quantitative noise
	if(HQUANT)
	{
		# Poisson noise is added to the original abundance. The higher the abundance, the more variance there is around the value, as the Poisson distribution states.

		RPOIS<-function(SPij) rpois(1,SPij)
		for(i in 1:sum(NSP))
		{
			SP<-apply(SP,c(1,2),RPOIS)
		}
	}
	# Set NAN to 0
	SP[is.na(SP)]=0

	SP[which(SP<0,arr.ind=TRUE)]<-0

	return(list(SP=apply(SP,c(1,2),round), PARM=PARM))
	# THE END #
}
