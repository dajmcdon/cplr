#' Generates the design matrix
#'
#' @param n number of observations
#' @param covMatSqrt square root of the design covariance, possibly produced by \code{\link{generateCov}}
#' @param rand character string giving the name of a function generating random numbers
#' @param m optional mean vector for the design matrix
#' @param ... optional named arguments passed to the random number generating function
#'
#' @return Returns a \code{n x ncol(covMatSqrt)} matrix.
#' @export
#'
#' @examples
#' generateX(100, diag(1,10), 'rnorm')
generateX <- function(n, covMatSqrt, rand='rnorm', m = NULL,...){
  rand = get(rand)
  X = matrix( rand(n * ncol(covMatSqrt),...), nrow=n) %*% covMatSqrt
  if(!is.null(m)) X = sweep(X, 2, m, '+')
  return(X)
}

#' Generates a vector of population regression coefficients
#'
#' @param p the number of coefficients
#' @param method one of eight named methods of generation, see Details below
#' @param b1 for \code{method='expdecay'}, see Details below
#' @param s for \code{method='sparse'}, see Details below
#' @param spar the l2 norm of the resulting coefficient vector, see Details below
#' @param ... optional arguments to random number generators, see Details below
#'
#' @return a vector of length p
#'
#' @details This function generates a vector of coefficients using one of eight methods.
#' \describe{
#'  \item{expdecay}{The result is \code{b1^(1:p)}}
#'  \item{lindecay}{p:1/p}
#'  \item{sparse}{Produces \code{s} coefficents drawn independently from the Laplace distribution, followed by \code{p-s} zeros}
#'  \item{randU}{uniform random coefficents, takes arguments in \code{...}}
#'  \item{randN}{Random Normal coefficients, takes arguments in \code{...}}
#'  \item{constNorm}{all coefficients are equal and the vector has l2 norm \code{spar}}
#'  \item{constLinDecay}{Same as \code{lindecay} except that the resulting vector has norm \code{spar}}
#'  \item{constant}{produces a p-vector of 1's}
#'  }
#' @export
#'
#' @examples
#' generateb(10, 'constNorm')
#' generateb(10, 'randN', mean=1, sd=.1)
#' generateb(10, 'expdecay', b1=.9)
#' generateb(10, 'expdecay', b1=1) # same as 'constant'
#' generateb(10, 'constant')
generateb <- function(p, method = 'randN', b1 = NULL, s = 10, spar = 10,...) {
  ## default method is 'sparse' b, s, is the number of nonzero entries
  ## b1 is for exponential decay, with |b1|<1
  ## other methods are 'lindecay', 'randN', and 'randU' with the later 2
  ## accepting ... arguments
  if(method=='sparse' && s>p){
    warning('The number of non-zero elements (s) exceeds the desired length (p). The result has s=p)')
    s=p
  }
  if(method=='expdecay' && is.null(b1)){
    stop('For expdecay, user must provide an initial coefficient (b1!=NULL)')
  }
  b = switch(method,
             expdecay = b1^(1:p),
             lindecay = p:1 / p,
             sparse = as.vector(Matrix::rsparsematrix(p,1,nnz=s,rand.x = rnorm)/sqrt(s)),
             randU = runif(p,...),
             randN = rnorm(p,...),
             constNorm = rep_len(1 / sqrt(p), p) * spar,
             constLinDecay = (p:1) / sqrt(sum((1:p)^2)) * spar,
             constant = rep(1, p),
             ones = rep(c(1,-1), times=ceiling(p/2))[1:p],
             stop('Invalid b type')
  )
  return(b)
}
#' Generates the response variable
#'
#' @param X a design matrix, possibly prouced by \code{\link{generateX}}
#' @param b a coefficient vector, possibly produced by \code{\link{generateb}}
#' @param rand character string giving the name of a function generating random numbers
#' @param ... optional named arguments passed to the random number generating function
#'
#' @return Returns a vector of length \code{nrow(X)}.
#' @export
#'
#' @examples
#' X = generateX(100, diag(1,10), 'rnorm')
#' b = generateb(10, 'constNorm')
#' Y = generateY(X, b, 'rnorm', sd=.1)
generateY <- function(X, b, rand, ...){
  if(p<-ncol(X) != length(b)){
    warning('Coefficient vector has been extended to match the columns of X')
    b = rep_len(b, p)
  }
  rand = get(rand)
  noise = rand(nrow(X),...)
  Y = X %*% b + noise
  return(Y)
}

#' Generates the square root of the covariance matrix
#'
#' @param p the number of predictors
#' @param type \code{'diag'} (the default) gives the identity matrix, \code{'corrUniform'} gives a matrix with 1 on the diagonal and \code{rho} off the diagonal, and \code{'corrDecay'} gives the AR covariance matrix (1 on the diagonal, \eqn{\rho^|i-j|} off the diagonal)
#' @param rho the correlation, ignored for \code{type='diag'}
#'
#' @return Returns the square root of the covariance
#' @export
#'
#' @examples
#' generateCov(10,'diag')
#' generateCov(10, 'corrDecay')
#' generateCov(10, 'corrDecay',.9)
generateCov <- function(p, type=c('diag','corrUniform','corrDecay'), rho=.8){
  ## returns the squareroot of the covariance matrix
  ## current methods are 'diag','corrUniform', and 'corrDecay'
  covMat = switch(type,
                  diag = diag(p),
                  corrUniform = toeplitz(c(1, rep(rho, p-1))),
                  corrDecay = toeplitz(c(1, rho^(1:(p-1)))),
                  stop('Invalid covariance type')
  )
  ev = eigen(covMat, symmetric = TRUE)
  covMatSqrt = t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
  return(covMatSqrt)
}

#' Generates standard Lapace random variates
#'
#' @param n the number of desired random numbers
#'
#' @return an n-vector of independent Laplace random variates
rlaplace <- function(n){
  u = runif(n)
  return(-sign(u-0.5) * (log(2) + ifelse(u < 0.5, log(u), log1p(-u))))
}
