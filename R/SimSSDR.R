#' @title Simulating the environment and species composition in a deterministique environment
#'
#' @description Generates species and environment data similarly to what had been done for the Ecology Note simulations on beta diversity in Legendre etal (2002) and (2004). It is derived from the programs SimSSDx written in Fortran from 1999 to 2004 by Pierre Legendre. R code from Marie-Helene Ouellette, calling the Fortran SIMSG routine of (Deutsch and Journel 1992) for autocorrelation structure.
#'
#' @param nx Number of pixels on the x axis of the simulated surface. Default = 100.
#' @param ny Number of pixels on the y axis of the simulated surface. Default = 100.
#' @param deter Vector (length = number of i gradients to simulate) that determines the type of deterministic structure in the resources variables. 0 means no deterministic structure (default), only N(0,\code{varnor}). 1 means gradient from top (low values) to bottom (high values) plus N(0,\code{varnor}). 2 means gradients in two directions plus N(0,\code{varnor}). 3 means one big patch in center of the field plus N(0,/code{varnor}). 4 means 4 waves plus N(0,\code{varnor}). 5 means two zones separated by a sharp step plus N(0,\code{varnor}).
#' @param height Vector (same length as deter) of maximum value of the deterministic structures.
#' @param beta List (size = length(deter)) of vectors of regression parameters between environment and species. Default=1 for all species. Each listed vector has nsp2[i] values.
#' @param nug1 Vector (length = length(deter)) of nugget values (default 0) for E (environment)
#' @param range11 Vector (length = length(deter)) of range value (0 if none) in the horizontal direction for Ei. Default is 10.
#' @param range12 Vector (length = length(deter)) of range in vertical direction (different from range11 if anisotropy) for Ei. Default is 10.
#' @param nug2 Vector (length = length(deter)) of nugget value (default 0) for R (response, i.e. species). This value is the same for all species on one gradient, but can change from one gradient to another.
#' @param range21 Vector (length = length(deter)) of range value (0 if none) in the horizontal direction for the species. Default is 10. This value is the same for all species on one gradient, but can change from one gradient to another.
#' @param range22 Vector (length = length(deter)) of range in vertical direction (different from range11 if anisotropy) for the species. Default is 10. This value is the same for all species on one gradient, but can change from one gradient to another.
#' @param varnor List (size = length(deter)) of vectors of variance for normal errors. The first value of each vector is for the environment, the following are for the species. There is 1 + nsp1 + nsp2 values given in that order. Default is 1 for all.
#' @param splike List (size = length(deter)) of vectors of transformation for species data: (-1) None, (0) Standardize only; no further transformation, (1) y = int(exp(y)), (2) y = int(sqrt(exp(y))), (3) y = 0 or 1 (binary). Default is (2) for all species.
#' @param cste List (size = length(deter)) of vectors of constants to Multiply each species' RND by, (1, default) : keep the species' RNDs as they are.
#' @param nsp1 Vector (length = length(deter)) of the number of species related to an environmental variable. Default is 5.
#' @param nsp2 Vector (length = length(deter)) of the number of species not related to an environmental variable. Default is 5.
#' @param SAE Vector (length = length(deter)) of logicals : \code{TRUE} if a spatial autocorrelation structure should be added to the environment. If not, \code{FALSE}.
#' @param SAR Vector (length = length(deter)) of logicals : \code{TRUE} if a spatial autocorrelation structure should be added to each species. If not, \code{FALSE}.
#' @param centroide Vector of size 2 giving the position of the center of the multinormal distribution when using \code{deter}=3. Default is 0, putting the center of the multinormal distribution at the center of the grid.
#' @param BETAA If NA (default), no COMPAS species added, if a list, use \code{Rcompas} to simulated species data. See code/{Rcompas} for details on the structure of this list.
#' @param \dots Arguments to be passed when autocorrelation is used for the species with niche and generalized beta distribution : see \code{Rcompas} for more details.
#'
#' @return
#'
#' \itemize{
#'   \item{E}{Simulated environmental data}
#'   \item{R}{Simulated response data, i.e. species (R1,R2,\code{dot}}
#'   \item{x}{spatial coordinates x}
#'   \item{y}{spatial coordinate y}
#' }
#'
#' @author Marie-Hélène Ouellette
#'
#' @references
#'
#' Legendre, P., M. R. T. Dale, M.-J. Fortin, J. Gurevitch, M. Hohn and D. Myers. 2002. The consequences of spatial structure for the design and analysis of ecological field surveys. \emph{Ecography} \strong{25}: 601-615.
#'
#' Legendre, P., M. R. T. Dale, M.-J. Fortin, P. Casgrain and J. Gurevitch.  2004. Effects of spatial structures on the results of field experiments. \emph{Ecology} \strong{85}: 3202-3214.
#'
#' Deutsch, C. V., and A. G. Journel. 1992. GSLIB: Geostatistical software library and user s guide. Oxford University Press, new York, New Yorkm USA.
#'
#' @examples
#'
#'# First example :
#'
#' res1<-SimSSDR() # All default settings
#'
#' # Complete example using both SIMSSD and COMPAS features
#'
#' # -- EX2 : complete example of how the list should be built
#' BETAA<-list() # 2 gradients
#' # P1 has 5 species
#' BETAA$P1$A<-list(list('value'), list(c(20,30,50,30,40)),list(c('n=5')))
#' BETAA$P1$m<-list(list('value','value'), list(c(0,1,2,3,4),c(2,3,4,5,6)), list(c('n=5'),c('n=5')))
#' BETAA$P1$r<-list(list('value','value'), list(c(2,3,2,3,2),c(2,3,2,3,2)), list(c('n=5'),c('n=5')))
#' BETAA$P1$alpha<-list(list('value','random'), list(c(2,3,4),'lognormal'), list(c('n=5'),c('n=5, meanlog = 0, sdlog = 1')))
#' BETAA$P1$gamma<-list(list('random','random'), list('normal','lognormal'), list(c('n=5,mean=0,sd=1'),c('n=5, meanlog = 0, sdlog = 1')))
#' # P2 has 6 species
#' BETAA$P2$A <- list(list('random'), list('normal'), list(c('n=6, mean=20, sd=1')))
#' BETAA$P2$m <- list(list('value','value'), list(c(1,2,3,4,5,6),c(1,2,3,4,5,6)), list(c('n=6'),c('n=6')))
#' BETAA$P2$r <- list(list('value','value'), list(c(1,2,3,3,2,1),c(2,2,2,2,2,2)), list(c('n=6'),c('n=6')))
#' BETAA$P2$alpha <- list(list('random','random'), list('normal','lognormal'), list(c('n=6,mean=0.5,sd=1'),c('n=6, meanlog = 0.01, sdlog = 1')))
#' BETAA$P2$gamma <- list(list('random','value'), list('normal',c(0.5,0.6,0.7,0.8,0.9,1)), list(c('n=6,mean=3,sd=1'),c('n=6')))
#'
#' ## --- Example with plotted results --- ##
#' nx=10
#' ny=10
#' deter=c(1,2)
#' height=c(5,6)
#' beta=list(c(5,5,5),c(5,5))
#' nug1=c(0,0)
#' range11=c(10,10)
#' range12=c(10,10)
#' nug2=c(0,0)
#' range21=c(10,10)
#' range22=c(10,10)
#' varnor=list(c(0,0,0,0,0,0,0),c(0,0,0,0,0))
#' nsp1=c(3,2)
#' nsp2=c(3,2)
#' SAE=c(TRUE,TRUE)
#' SAR=c(TRUE,TRUE)
#' splike=list(c(-1,0,1,2,3,3),c(1,2,3,3))
#' cste=list(c(1,2,3,4,5,6),c(1,2,3,4))
#' centroide=list(0,0)
#'
#' G<-c(TRUE,TRUE) # The two profiles are two guilds
#' MAJOR<-FALSE
#' HISI<-list(n=2,g=TRUE,mc=999,mv=999)
#' HQUAL<-FALSE
#' HQUANT<-FALSE
#' HCCA<-list(type=0)
#' HSYM<-c(TRUE,FALSE) # For the two profiles
#'
#' res<-SimSSDR(nx=nx,ny=ny,deter=deter,height=height,beta=beta,nug1=nug1,range11=range11,range12=range12,nug2=nug2,range21=range21,range22=range22,varnor=varnor,nsp1=nsp1,nsp2=nsp2,splike=splike,cste=cste,SAE=SAE,SAR=SAR,centroide=centroide,BETAA=BETAA,G=G,MAJOR=MAJOR,HISI=HISI,HQUAL=HQUAL,HQUANT=HQUANT,HCCA=HCCA,HSYM=HSYM)
#'
#' # Plot to test
#'
#'
#' # for all deterministic structure
#' for(i in 1:(length(deter))){
#'   k=1
#'   ##  Environment  ##
#'   name=paste('Deterministic structure ',i)
#'   pos<-res$E[,i]>=0
#'   symbols(res$x[pos],res$y[pos],circles=res$E[,i][pos]/max(abs(res$E[,i]))/2,fg='grey',inches=FALSE)
#'   neg<-res$E[,i]<0
#'
#'   if(any(neg)){
#'     symbols(res$x[neg],res$y[neg],circles=abs(res$E[,i][neg]/max(abs(res$E[,i])))/2,fg='black',add=TRUE,inches=FALSE)
#'   }
#'
#'   ##  Species  ##
#'   for(j in 1:ncol(res$R)){
#'     name=paste('Species ',j, ' : Dtmc. sctr. ',i)
#'	   pos<-res$R[,k]>=0
#' 	   symbols(res$x[pos],res$y[pos],circles=res$R[pos,k]/max(abs(res$R[,k]))/2,fg='grey',inches=FALSE,main=name)
#'	   neg<-res$R[,k]<0
#'	   if(any(neg)){
#'	     symbols(res$x[neg],res$y[neg],circles=abs(res$R[neg,k]/max(abs(res$R[,k])))/2,fg='black',add=TRUE,inches=FALSE)
#'	   }
#'       k=k+1
#'   }
#' }
#'
#' @keywords datagen
#' @export
SimSSDR<-function(nx=100,ny=100,deter=0,height=5,nug1=0,range11=10,range12=10,nug2=0,range21=10,range22=10,nsp1=5,nsp2=5,cste=list(mat.or.vec(1,nsp1+nsp2)+1),varnor=list(mat.or.vec(1,nsp1+nsp2+1)+1),beta=list(mat.or.vec(1,nsp1)+1),splike=list(mat.or.vec(1,nsp1+nsp2)+2),SAE=TRUE,SAR=TRUE,centroide=0,BETAA=NA,...){
# For RCOMPAS : BETAA,X,G,MAJOR=FALSE,HISI,HQUAL,HQUANT,HCCA,HSYM,SAR

# BEFORE RCOMPAS
# nsp3=0,discrete=TRUE,nsp3CENT=NULL,nsp3RANGE=NULL,nsp3HEIGHT=NULL,nsp3SD=NULL
	if(nx!=ny) stop('nx!=ny : Implemented only for square surfaces')


	for(i in 1:length(deter)) # For each gradient
	{
		# Vectors

		if(!any(deter[i]==c(0:5))) stop(paste('deter arguments takes only values from 0 to 5. Check value for gradient',i))

		if(nug1[i]<0) stop(paste('nug1 must be equal or greater than 0. Check value for gradient',i))

		if(nug2[i]<0) stop(paste('nug2 must be equal or greater than 0. Check value for gradient',i))

		if(range11[i]<0) stop(paste('range11 must be equal or greater than 0. Check value for gradient',i))

		if(range12[i]<0) stop(paste('range12 must be equal or greater than 0. Check value for gradient',i))

		if(range21[i]<0) stop(paste('range21 must be equal or greater than 0. Check value for gradient',i))

		if(range22[i]<0) stop(paste('range22 must be equal or greater than 0. Check value for gradient',i))

		if(height[i]<=0) stop(paste('height argument must be positive. Check value for gradient',i))

		if(nsp1[i]<0) stop(paste('nsp1 must be equal or greater than 0. Check value for gradient',i))

		if(nsp2[i]<0) stop(paste('nsp2 must be equal or greater than 0. Check value for gradient',i))

		# Lists

		if(length(beta[[i]])!=nsp1[i]) stop(paste('length(beta) must be equal to nsp1. Check value for gradient',i))

		if(!(length(centroide[[i]])==2||centroide[[i]]==0)) stop(paste('centroide was not properly defined. Check value for gradient',i))

		if(length(cste[[i]])!=(nsp1[i]+nsp2[i])) stop(paste('size of cste is not proper for gradient',i))

		if(length(splike[[i]])!=(nsp1[i]+nsp2[i])) stop(paste('size of splike is not properly defined for gradient',i))

		if(length(varnor[[i]])!=(nsp1[i]+nsp2[i]+1)) stop(paste('size of varnor is not proper for gradient',i))

			for(j in 1:length(deter))
			{
				if(!any(splike[[i]][j]==c(-1,0,1,2,3))) stop(paste('a value for splike was not properly defined. Check values for gradient',i))

				if(any(varnor[[i]][j]<0)) stop(paste('varnor should be positive for all. Check values for gradient',i))
			}

	}

	E<-as.matrix(mat.or.vec(ny*nx,length(deter)))
	x<-rep(c(1:nx),ny)
	y<-rep(c(1:ny),each=nx)

	SP1<-mat.or.vec(ny*nx,1)


	# FOR ALL DETERMINISTIC STRUCTURE
	for(d in 1:length(deter))
	{
		varnord_count<-1
		# --- Simluating Evironmental data : E = D(deter) + SA + error --- #

		# -- First take care of deterministic structure of the environment

		E[,d]<-DETER(deter[d],height[d],varnor[[d]][varnord_count],nx,ny,centroide[[d]])
		varnord_count=varnord_count+1

		# -- Add autocorrelation if requested
		if(SAE[d])
		{
			result<-vector(mode = "numeric", length = nx*ny)
			SAE_vec<-.Fortran("inpsgsm",as.integer(nx*ny),as.integer(nx),as.integer(ny),as.numeric(nug1),as.double(range11),as.double(range12),as.double(result),PACKAGE="RsimSSDCOMPAS")
			names(SAE_vec)<-c('nxymax','nx','ny','nugget','range1','range2','result')

			E[,d]=E[,d]+SAE_vec$result
		}

		# Simulation Species response (linked to environment) : R = E*beta + SA + error
		R<-mat.or.vec(ny*nx,1)

		if(nsp1[d]>0)
		{
			for(i in 1:nsp1[d])
			{

				rand_norm<-rnorm(ny*nx,0,varnor[[d]][varnord_count])
				varnord_count=varnord_count+1

				Ri<-E[,d]*beta[[d]][i]+rand_norm
				R<-cbind(R,Ri)
			}
		}

		# Simulation Species response (not linked to environment)
		if(nsp2[d]>0)
		{
			for(k in 1:nsp2[d])
			{

				rand_norm<-rnorm(ny*nx,0,varnor[[d]][varnord_count])
				varnord_count=varnord_count+1

				Ri<-rand_norm
				R<-cbind(R,Ri)

			}
		}

		# BEFORE RCOMPAS
		# Simulation Species response (linked to environment with niche)
		#if(nsp3>0)
		#{
		#	for(l in 1:nsp3)
		#	{
		#		rand_norm<-rnorm(ny*nx,0,varnor)
		#		Centi<-nsp3CENT[l]*nx*ny/100
		#		nN=nsp3RANGE[l]*nx*ny/100
		#		N_gr<-bkde(rnorm(nN,mean = 50, sd = nsp3SD[l]),gridsize=nN)$y

		#							N<-c(vector(length=(Centi-nN/2)),N_gr,vector(length=(nx*ny-length(N_gr)-length(vector(length=(Centi-nN/2))))))*100
			# Make distribution at maximum height
			# Rl<-(N+rand_norm)*nsp3HEIGHT[l]/max(N+rand_norm)
		#		Rl<-N+rand_norm
		#		R<-cbind(R,Rl)

		#	}
		#}

		# -- Add autocorrelation if requested
		if(SAR[d] & (nsp1[d] + nsp2[d])>0)
		{
			for(j in 1:ncol(R))
			{
				resultj<-vector(mode = "numeric", length = nx*ny)
				SAR_vecj<-.Fortran("inpsgsm",as.integer(nx*ny),as.integer(nx),as.integer(ny),as.numeric(nug2),as.double(range21),as.double(range22),as.double(resultj),PACKAGE="RsimSSDCOMPAS")
				names(SAR_vecj)<-c('nxymax','nx','ny','nugget','range1','range2','result')

				R[,j]=as.vector(R[,j])+as.vector(SAR_vecj$result)
			}
		}

		if((nsp1[d] + nsp2[d])>0)
		{
			R<-R[,-1]
			# Make abundances have the characteristics defined by the user, and at proper height for niche species

				for(m in 1:(nsp1[d]+nsp2[d]))
				{


					#(-1) None
					#(0) Standardize only; no further transformation'
					#(1) y = int(exp(y))
					#(2) y = int(sqrt(exp(y)))
					#(3) y = 0 or 1 (binary)

					# First, standardize the species
					if(splike[[d]][m]>=0)
					{
						R[,m]<-scale(R[,m])

						if(splike[[d]][m]>0)
						{
							# Multiply it by a random normal vector and the constant
							tmp_R=R[,m]*rnorm(1)*cste[[d]][m]

							# Make the form of distribution the one the user chose
							if(splike[[d]][m]==1) R[,m]=as.integer(exp(tmp_R))
            				if(splike[[d]][m]==2) R[,m]=as.integer(sqrt(exp(tmp_R)))
							if(splike[[d]][m]==3)
							{
								R[which(as.integer(exp(tmp_R))!=0),m]=1
								R[which(as.integer(exp(tmp_R))==0),m]=0
							}
						}
					}
				}
			}

# 			BEFORE RCOMPAS
#			z=1
#			if(nsp3>0)
#			{
#				for(n in (nsp1+nsp2+2):(nsp1+nsp2+nsp3+1))
#				{
#					# Make maximum value the heigth assigned by user
#					Rn_heigth<-R[,n]*nsp3HEIGHT[z]/max(R[,n])
#
#					# Choice #1 : add the absolute minimum value to 0, bringing all abundances up.
#					#if(min(Rn_heigth)<0)
#					#{
#					#	R[,n]<-round(Rn_heigth+abs(min(Rn_heigth)))
#					#}
#					#if(min(Rn_heigth)>=0)
#					#{
#					#	R[,n]<-round(Rn_heigth)
#					#}
#					z=z+1
#					# Choice #2 : bring all negative values to 0.
#					R[,n]<-round(Rn_heigth)
#					R[which(R[,n]<0),n]=0
#
#				}
#			}

		SP1<-cbind(SP1,R)

	}
	if((nsp1[1] + nsp2[1])>0) SP<-SP1[,-1]
	if((nsp1[1] + nsp2[1])==0) SP<-mat.or.vec(ny*nx,1)


	if(!is.vector(SP))
	{
		colnames(SP)<-paste('SP',1:ncol(SP),sep='')
	}


	if(!is.vector(E))
	{
		colnames(E)<-paste('ENV',1:ncol(E),sep='')
	}
	res<-list(E,SP,x,y)
	names(res)<-c('E','R','x','y')
	# Add RCOMPAS species : flexible species with niches

	if(!is.na(BETAA))
	{
		SP_compas<-Rcompas(BETAA=BETAA,X=E,G=G,MAJOR=MAJOR,HISI=HISI,HQUAL=HQUAL,HQUANT=HQUANT,HCCA=HCCA,HSYM=HSYM,SAR=SAR[1])
		SP<-cbind(SP,SP_compas$SP)
		if((nsp1[1] + nsp2[1])==0) SP<-SP[,-1]
		if(!is.vector(SP))
		{
			colnames(SP)<-paste('SP',1:ncol(SP),sep='')
		}
		res<-list(E,SP,x,y,SP_compas$PARM)
		names(res)<-c('E','R','x','y','PARM')
	}

	return(res)
}

##################################################
# ---  Simulation the deterninistic effect --- ###
##################################################

DETER<-function(deter,height,varnor,nx,ny,centroide)
{
	if(deter==0)
	{
		rand_norm<-rnorm(ny*nx,0,varnor)
		E<-rand_norm
		# leveling to height
		E<-E*(height/max(E))
	}
	if(deter==1)
	{
		D_NtoS<-c(1:ny)
		# leveling to height
		D_NtoS_rep<-rep(D_NtoS,nx)
		# leveling to height
		D<-matrix(D_NtoS_rep,ny,nx)*(height/max(D_NtoS))
		rand_norm<-rnorm(ny*nx,0,varnor)
		rand_norm_mat<-matrix(rand_norm,ny,nx)
		E<-D+rand_norm_mat
		E<-as.vector(t(E))
	}
	if(deter==2)
	{

		D_NtoS<-c(1:ny)
		D_NtoS_rep<-rep(D_NtoS,nx)
		D_NtoS_mat<-matrix(D_NtoS_rep,ny,nx)
		D_WtoE<-c(1:nx)
		D_WtoE_rep<-rep(D_WtoE,ny)
		D_WtoE_mat<-matrix(D_WtoE_rep,ny,nx,byrow=TRUE)
		D=D_NtoS_mat+D_WtoE_mat
		# leveling to height
		D<-D*(height/max(D))
		rand_norm<-rnorm(ny*nx,0,varnor)
		rand_norm_mat<-matrix(rand_norm,ny,nx)
		E<-D+rand_norm_mat
		E<-as.vector(t(E))
	}
	if(deter==3)
	{
		if(centroide==c(0,0))
		{
			bivn <- mvrnorm(nx*ny, mu = c(100, 100), Sigma = matrix(c(1, 0, 0, 1), 2))
			bivn.kde <- kde2d(bivn[,1], bivn[,2], n = nx)
			D=bivn.kde$z*100
			# leveling to height
			D<-D*(height/max(D))
			# Plotting to verify
			#x<-rep(c(1:nx),ny)
			#y<-rep(c(1:ny),each=nx)
			#symbols(x, y, circles=as.vector(bivn.kde$z))
			rand_norm<-rnorm(ny*nx,0,varnor)
			rand_norm_mat<-matrix(rand_norm,ny,nx)
			E<-D+rand_norm_mat
			# Plotting to verify
			#x<-rep(c(1:nx),ny)
			#y<-rep(c(1:ny),each=nx)
			#symbols(x, y, circles=as.vector((E+abs(min(E))+0.1)/10), inches=FALSE)
			E<-as.vector(t(E))
		}
		if(centroide!=0)
		{
			# Enlarging grid to reposition center later on
			grid_size<-max(c(abs(centroide-nx),abs(nx-((-centroide)%%nx+1))))*2+1
			bivn <- mvrnorm(grid_size*grid_size, mu = c(100, 100), Sigma = matrix(c(1, 0, 0, 1), 2))
			bivn.kde <- kde2d(bivn[,1], bivn[,2], n = grid_size)
			# Resizing to original grid size, leaving the center of the distribution at the centroide chosen by the user
			UP<-centroide[2]<=nx/2
			DOWN<-centroide[2]>nx/2
			RIGHT<-centroide[1]<=nx/2
			LEFT<-centroide[1]>nx/2

			if(UP&&RIGHT)
			{
				D=bivn.kde$z[c((1:nx)+grid_size-nx),c((1:nx)+grid_size-nx)]*100
			}
			if(UP&&LEFT)
			{
				D=bivn.kde$z[c((1:nx)+grid_size-nx),1:nx]*100
			}
			if(DOWN&&RIGHT)
			{
				D=bivn.kde$z[1:nx,c((1:nx)+grid_size-nx)]*100
			}
			if(DOWN&&LEFT)
			{
				D=bivn.kde$z[1:nx,1:nx]*100
			}

			# leveling to height
			D<-D*(height/max(D))
			#Plotting to verify
			#x<-rep(c(1:nx),ny)
			#y<-rep(c(1:ny),each=nx)
			#symbols(x, y, circles=as.vector(D))
			rand_norm<-rnorm(ny*nx,0,varnor)
			rand_norm_mat<-matrix(rand_norm,ny,nx)
			E<-D+rand_norm_mat
			# Plotting to verify
			x<-rep(c(1:nx),ny)
			y<-rep(c(1:ny),each=nx)
			symbols(x, y, circles=as.vector((t(E)+abs(min(E))+0.1)/10), inches=FALSE)
			E<-as.vector(t(E))
		}

	}
	if(deter==4)
	{
		#Initialisation of one column of D
		D_1<-mat.or.vec(ny,1)
		# Top of waves
		top_pos<-ceiling(c(ny/4-ny/8,2*ny/4-ny/8,3*ny/4-ny/8,ny-ny/8))
		# Top of waves height
		top_height<-ceiling(ny/4)
		# Filling in between waves
		for(i in 1:length(top_pos))
		{
			# -- First, lets put the top heights where they belong
			D_1[top_pos]=top_height

			# -- Then fill from wave down
			# counter for position in D
			k=top_pos[i]
			# counter for height value
			j=1
			while( k+1<=ny && D_1[k+1]==0)
			{
				D_1[k+1]=top_height-j
				j=j+1
				k=k+1
			}

			# -- Then fill from wave up, stopping in the crest of the wave
			# counter for position in D
			k=top_pos[i]
			# counter for height value
			j=1
			while( k-1>0 && D_1[k-1]<top_height/2)
			{
				D_1[k-1]=top_height-j
				j=j+1
				k=k-1
			}
		}

		# Now build the matrix D
		D<-matrix(rep(D_1,times=nx),ny,nx)
		# leveling to height
		D<-D*(height/max(D))
		rand_norm<-rnorm(ny*nx,0,varnor)
		rand_norm_mat<-matrix(rand_norm,ny,nx)
		E<-D+rand_norm_mat
		# Plotting to verify
		#x<-rep(c(1:nx),ny)
		#y<-rep(c(1:ny),each=nx)
		#symbols(x, y, circles=as.vector((t(E)+abs(min(E))+0.1)/10), inches=FALSE)
		E<-as.vector(t(E))
	}
	if(deter==5)
	{
		#Initialisation of one column of D
		D_1<-mat.or.vec(ny,1)
		# Discontinuity position : one fifth of points
		disc_pos<-ceiling(ny/2)
		wide<-ceiling(ny/5)

		D_1[c(1:(disc_pos-ceiling(ny/5/2)),c(disc_pos+ceiling(ny/5/2)):ny)]<-3
		# Now build the matrix D
		D<-matrix(rep(D_1,times=nx),ny,nx)
		# leveling to height
		D<-D*(height/max(D))
		rand_norm<-rnorm(ny*nx,0,varnor)
		rand_norm_mat<-matrix(rand_norm,ny,nx)
		E<-D+rand_norm_mat
		# Plotting to verify
		#x<-rep(c(1:nx),ny)
		#y<-rep(c(1:ny),each=nx)
		#symbols(x, y, circles=as.vector((t(E)+abs(min(E))+0.1)/10), inches=FALSE)
		E<-as.vector(t(E))
	}
return(E)
}
