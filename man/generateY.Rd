% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateModel.R
\name{generateY}
\alias{generateY}
\title{Generates the response variable}
\usage{
generateY(X, b, rand, ...)
}
\arguments{
\item{X}{a design matrix, possibly prouced by \code{\link{generateX}}}

\item{b}{a coefficient vector, possibly produced by \code{\link{generateb}}}

\item{rand}{character string giving the name of a function generating random numbers}

\item{...}{optional named arguments passed to the random number generating function}
}
\value{
Returns a vector of length \code{nrow(X)}.
}
\description{
Generates the response variable
}
\examples{
X = generateX(100, diag(1,10), 'rnorm')
b = generateb(10, 'constNorm')
Y = generateY(X, b, 'rnorm', sd=.1)
}
