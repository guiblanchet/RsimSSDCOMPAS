% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EVAL.R
\name{EVAL}
\alias{EVAL}
\title{Build a matrix of parameters}
\usage{
EVAL(BETAAjx, MATi, i, j, LINE, col)
}
\arguments{
\item{BETAAjx}{Piece of a BETAA list that contains a number of parameters that will be used later on to evaluate species abundances from the Rcompas function.}

\item{MATi}{Matrix of parameters that will be modified, build by this functions.}

\item{i}{Number corresponding the gradient.}

\item{j}{Number corresponding to the profile.}

\item{LINE}{Number corresponding to the species, the line number of Mati to fill.}

\item{col}{The column of Mati to fill.}
}
\value{
The matrix Mati filled
}
\description{
Build a matrix of parameters from the list BETAA for the main function Rcompas
}
\author{
Marie-Hélène Ouellette
}
\keyword{datagen}
