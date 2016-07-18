#' Generates data as in the paper
#'
#' @param ntrain number of training samples
#' @param ntest number of testing samples
#' @param p number of predictors
#' @param covType for use with \code{\link{generateCov}}
#' @param rho for use with \code{\link{generateCov}}
#' @param bmethod for use with \code{\link{generateb}}
#' @param SNR if not na, sets the standard deviation of the noise to \eqn{\sqrt{b'Cov b} / SNR}
#' @param noisesd the standard deviation of the noise if \code{SNR=NA}
#' @param ... optional arguments passed to \code{generateb}
#'
#' @return A list with components \code{Xtrain}, \code{Ytrain}, \code{Xtest}, \code{Ytest}, \code{bstar}, and \code{sig} giving respectively the training data, testing data, population regression coefficients, and noise standard deviation. The matrix \code{Xtest} is the same as the training data adjusted to the correct size (either by truncating or recycling) so that predictions are at the same design points. Note that this function is simply a wrapper of the various \code{generate*} functions for ease of simulation. If more flexibility is desired, use those directly or similar.
#' @export
#'
#' @examples
#' createData(2500, 2, 250, 'corrUniform', .8, 'constant', NA, 2500-250-1)
createData <- function(ntrain, ntest, p, covType, rho, bmethod,
                       SNR, noisesd, ...){
  cms = generateCov(p, covType, rho)
  Xtrain = generateX(ntrain, cms)
  b = generateb(p, bmethod,...)
  if(!is.na(SNR)) noisesd = sqrt(crossprod(cms  %*% b) / SNR)
  Ytrain = generateY(Xtrain, b, 'rnorm', sd=noisesd)
  if(ntest > ntrain){
    Xtest = cbind(Xtrain, generateX(ntest-ntrain, cms))
  } else {
    Xtest = Xtrain[1:ntest,,drop=FALSE]
  }
  Ytest = generateY(Xtest, b, 'rnorm', sd=noisesd)
  return(list(Xtrain=Xtrain, Ytrain=Ytrain, Xtest=Xtest, Ytest=Ytest,
              bstar=b, sig=noisesd ))
}
