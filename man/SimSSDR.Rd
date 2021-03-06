% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimSSDR.R
\name{SimSSDR}
\alias{SimSSDR}
\title{Simulating the environment and species composition in a deterministique environment}
\usage{
SimSSDR(nx = 100, ny = 100, deter = 0, height = 5, nug1 = 0,
  range11 = 10, range12 = 10, nug2 = 0, range21 = 10, range22 = 10,
  nsp1 = 5, nsp2 = 5, cste = list(mat.or.vec(1, nsp1 + nsp2) + 1),
  varnor = list(mat.or.vec(1, nsp1 + nsp2 + 1) + 1),
  beta = list(mat.or.vec(1, nsp1) + 1), splike = list(mat.or.vec(1, nsp1 +
  nsp2) + 2), SAE = TRUE, SAR = TRUE, centroide = 0, BETAA = NA, ...)
}
\arguments{
\item{nx}{Number of pixels on the x axis of the simulated surface. Default = 100.}

\item{ny}{Number of pixels on the y axis of the simulated surface. Default = 100.}

\item{deter}{Vector (length = number of i gradients to simulate) that determines the type of deterministic structure in the resources variables. 0 means no deterministic structure (default), only N(0,\code{varnor}). 1 means gradient from top (low values) to bottom (high values) plus N(0,\code{varnor}). 2 means gradients in two directions plus N(0,\code{varnor}). 3 means one big patch in center of the field plus N(0,/code{varnor}). 4 means 4 waves plus N(0,\code{varnor}). 5 means two zones separated by a sharp step plus N(0,\code{varnor}).}

\item{height}{Vector (same length as deter) of maximum value of the deterministic structures.}

\item{nug1}{Vector (length = length(deter)) of nugget values (default 0) for E (environment)}

\item{range11}{Vector (length = length(deter)) of range value (0 if none) in the horizontal direction for Ei. Default is 10.}

\item{range12}{Vector (length = length(deter)) of range in vertical direction (different from range11 if anisotropy) for Ei. Default is 10.}

\item{nug2}{Vector (length = length(deter)) of nugget value (default 0) for R (response, i.e. species). This value is the same for all species on one gradient, but can change from one gradient to another.}

\item{range21}{Vector (length = length(deter)) of range value (0 if none) in the horizontal direction for the species. Default is 10. This value is the same for all species on one gradient, but can change from one gradient to another.}

\item{range22}{Vector (length = length(deter)) of range in vertical direction (different from range11 if anisotropy) for the species. Default is 10. This value is the same for all species on one gradient, but can change from one gradient to another.}

\item{nsp1}{Vector (length = length(deter)) of the number of species related to an environmental variable. Default is 5.}

\item{nsp2}{Vector (length = length(deter)) of the number of species not related to an environmental variable. Default is 5.}

\item{cste}{List (size = length(deter)) of vectors of constants to Multiply each species' RND by, (1, default) : keep the species' RNDs as they are.}

\item{varnor}{List (size = length(deter)) of vectors of variance for normal errors. The first value of each vector is for the environment, the following are for the species. There is 1 + nsp1 + nsp2 values given in that order. Default is 1 for all.}

\item{beta}{List (size = length(deter)) of vectors of regression parameters between environment and species. Default=1 for all species. Each listed vector has nsp2[i] values.}

\item{splike}{List (size = length(deter)) of vectors of transformation for species data: (-1) None, (0) Standardize only; no further transformation, (1) y = int(exp(y)), (2) y = int(sqrt(exp(y))), (3) y = 0 or 1 (binary). Default is (2) for all species.}

\item{SAE}{Vector (length = length(deter)) of logicals : \code{TRUE} if a spatial autocorrelation structure should be added to the environment. If not, \code{FALSE}.}

\item{SAR}{Vector (length = length(deter)) of logicals : \code{TRUE} if a spatial autocorrelation structure should be added to each species. If not, \code{FALSE}.}

\item{centroide}{Vector of size 2 giving the position of the center of the multinormal distribution when using \code{deter}=3. Default is 0, putting the center of the multinormal distribution at the center of the grid.}

\item{BETAA}{If NA (default), no COMPAS species added, if a list, use \code{Rcompas} to simulated species data. See code/{Rcompas} for details on the structure of this list.}

\item{\dots}{Arguments to be passed when autocorrelation is used for the species with niche and generalized beta distribution : see \code{Rcompas} for more details.}
}
\value{
\itemize{
  \item{E}{Simulated environmental data}
  \item{R}{Simulated response data, i.e. species (R1,R2,\code{dot}}
  \item{x}{spatial coordinates x}
  \item{y}{spatial coordinate y}
}
}
\description{
Generates species and environment data similarly to what had been done for the Ecology Note simulations on beta diversity in Legendre etal (2002) and (2004). It is derived from the programs SimSSDx written in Fortran from 1999 to 2004 by Pierre Legendre. R code from Marie-Helene Ouellette, calling the Fortran SIMSG routine of (Deutsch and Journel 1992) for autocorrelation structure.
}
\examples{

# First example :

res1<-SimSSDR() # All default settings

# Complete example using both SIMSSD and COMPAS features

# -- EX2 : complete example of how the list should be built
BETAA<-list() # 2 gradients
# P1 has 5 species
BETAA$P1$A<-list(list('value'), list(c(20,30,50,30,40)),list(c('n=5')))
BETAA$P1$m<-list(list('value','value'), list(c(0,1,2,3,4),c(2,3,4,5,6)), list(c('n=5'),c('n=5')))
BETAA$P1$r<-list(list('value','value'), list(c(2,3,2,3,2),c(2,3,2,3,2)), list(c('n=5'),c('n=5')))
BETAA$P1$alpha<-list(list('value','random'), list(c(2,3,4),'lognormal'), list(c('n=5'),c('n=5, meanlog = 0, sdlog = 1')))
BETAA$P1$gamma<-list(list('random','random'), list('normal','lognormal'), list(c('n=5,mean=0,sd=1'),c('n=5, meanlog = 0, sdlog = 1')))
# P2 has 6 species
BETAA$P2$A <- list(list('random'), list('normal'), list(c('n=6, mean=20, sd=1')))
BETAA$P2$m <- list(list('value','value'), list(c(1,2,3,4,5,6),c(1,2,3,4,5,6)), list(c('n=6'),c('n=6')))
BETAA$P2$r <- list(list('value','value'), list(c(1,2,3,3,2,1),c(2,2,2,2,2,2)), list(c('n=6'),c('n=6')))
BETAA$P2$alpha <- list(list('random','random'), list('normal','lognormal'), list(c('n=6,mean=0.5,sd=1'),c('n=6, meanlog = 0.01, sdlog = 1')))
BETAA$P2$gamma <- list(list('random','value'), list('normal',c(0.5,0.6,0.7,0.8,0.9,1)), list(c('n=6,mean=3,sd=1'),c('n=6')))

## --- Example with plotted results --- ##
nx=10
ny=10
deter=c(1,2)
height=c(5,6)
beta=list(c(5,5,5),c(5,5))
nug1=c(0,0)
range11=c(10,10)
range12=c(10,10)
nug2=c(0,0)
range21=c(10,10)
range22=c(10,10)
varnor=list(c(0,0,0,0,0,0,0),c(0,0,0,0,0))
nsp1=c(3,2)
nsp2=c(3,2)
SAE=c(TRUE,TRUE)
SAR=c(TRUE,TRUE)
splike=list(c(-1,0,1,2,3,3),c(1,2,3,3))
cste=list(c(1,2,3,4,5,6),c(1,2,3,4))
centroide=list(0,0)

G<-c(TRUE,TRUE) # The two profiles are two guilds
MAJOR<-FALSE
HISI<-list(n=2,g=TRUE,mc=999,mv=999)
HQUAL<-FALSE
HQUANT<-FALSE
HCCA<-list(type=0)
HSYM<-c(TRUE,FALSE) # For the two profiles

res<-SimSSDR(nx=nx,ny=ny,deter=deter,height=height,beta=beta,nug1=nug1,range11=range11,range12=range12,nug2=nug2,range21=range21,range22=range22,varnor=varnor,nsp1=nsp1,nsp2=nsp2,splike=splike,cste=cste,SAE=SAE,SAR=SAR,centroide=centroide,BETAA=BETAA,G=G,MAJOR=MAJOR,HISI=HISI,HQUAL=HQUAL,HQUANT=HQUANT,HCCA=HCCA,HSYM=HSYM)

# Plot to test


# for all deterministic structure
for(i in 1:(length(deter))){
  k=1
  ##  Environment  ##
  name=paste('Deterministic structure ',i)
  pos<-res$E[,i]>=0
  symbols(res$x[pos],res$y[pos],circles=res$E[,i][pos]/max(abs(res$E[,i]))/2,fg='grey',inches=FALSE)
  neg<-res$E[,i]<0

  if(any(neg)){
    symbols(res$x[neg],res$y[neg],circles=abs(res$E[,i][neg]/max(abs(res$E[,i])))/2,fg='black',add=TRUE,inches=FALSE)
  }

  ##  Species  ##
  for(j in 1:ncol(res$R)){
    name=paste('Species ',j, ' : Dtmc. sctr. ',i)
   pos<-res$R[,k]>=0
	   symbols(res$x[pos],res$y[pos],circles=res$R[pos,k]/max(abs(res$R[,k]))/2,fg='grey',inches=FALSE,main=name)
   neg<-res$R[,k]<0
   if(any(neg)){
     symbols(res$x[neg],res$y[neg],circles=abs(res$R[neg,k]/max(abs(res$R[,k])))/2,fg='black',add=TRUE,inches=FALSE)
   }
      k=k+1
  }
}

}
\references{
Legendre, P., M. R. T. Dale, M.-J. Fortin, J. Gurevitch, M. Hohn and D. Myers. 2002. The consequences of spatial structure for the design and analysis of ecological field surveys. \emph{Ecography} \strong{25}: 601-615.

Legendre, P., M. R. T. Dale, M.-J. Fortin, P. Casgrain and J. Gurevitch.  2004. Effects of spatial structures on the results of field experiments. \emph{Ecology} \strong{85}: 3202-3214.

Deutsch, C. V., and A. G. Journel. 1992. GSLIB: Geostatistical software library and user s guide. Oxford University Press, new York, New Yorkm USA.
}
\author{
Marie-Hélène Ouellette
}
\keyword{datagen}
