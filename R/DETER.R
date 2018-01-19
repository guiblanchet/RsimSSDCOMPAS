#' @title Simulation of the deterministic effect
#'
#' @description Simulation of the deterministic effect that structures species distribution
#'
#' @param nx Number of pixels on the x axis of the simulated surface.
#' @param ny Number of pixels on the y axis of the simulated surface.
#' @param deter Type of deterministic structure in resource variable ;(0), no deterministic structure, only N(0,varnor) is considered; (1), gradient from top (low values) to bottom (high values) + N(0,varnor); (2), gradients in two directions + N(0,varnor); (3), one big patch in center of the field + N(0,varnor); (4), 4 waves + N(0,varnor); (5), two zones separated by a sharp step + N(0,varnor).
#' @param height maximum value of the deterministic structure.
#' @param varnor variance of normal error.
#' @param centroide Vector of length 2 giving the position of the center of the multinormal distribution when using \code{deter}=3. Default is 0, putting the center of the multinormal distribution in the middle of the grid.
#'
#'
#' @return
#'
#' A vector E containing the environment with deterministic structure
#'
#' @author Marie-Hélène Ouellette
#'
#' @references
#' 	Legendre, P., M. R. T. Dale, M.-J. Fortin, J. Gurevitch, M. Hohn and D. Myers. 2002. The consequences of spatial structure for the design and analysis of ecological field surveys. \emph{Ecography} \strong{25}: 601-615.
#'
#'	Legendre, P., M. R. T. Dale, M.-J. Fortin, P. Casgrain and J. Gurevitch.  2004. Effects of spatial structures on the results of field experiments. \emph{Ecology} \strong{85} : 3202-3214.
#'
#' @keywords datagen
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
