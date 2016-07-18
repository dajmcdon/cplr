
#' Calculates estimation error for PLS
#'
#' This function is used by \code{\link{BatchExperiments}} to run the simulations.
#' See that package for documentation.
#'
#' @param job see \code{\link{BatchExperiments}}
#' @param static see \code{\link{BatchExperiments}}
#' @param dynamic see \code{\link{BatchExperiments}}
#' @param ... see \code{\link{BatchExperiments}}
#'
#' @return A list with components \code{estimRisk} and \code{ptm}. The first gives the estimation risk of the particular method while the second calculates how long the solution took to compute.
#' @seealso empEstimRisk
#' @export
plsEstimAlgo <- function(job, static, dynamic, ...){
  ptm = proc.time()
  bhat = pls(dynamic$Xtrain, dynamic$Ytrain, ...)
  ptm = proc.time()-ptm
  estimRisk = empEstimRisk(bhat$bhat, dynamic$bstar)
  return(list(estimRisk = estimRisk, ptm=ptm[3]))
}

#' Calculates estimation error for OLS
#'
#' This function is used by \code{\link{BatchExperiments}} to run the
#' simulations. See that package for documentation.
#'
#' @inheritParams plsEstimAlgo
#'
#' @return A list with components \code{estimRisk} and \code{ptm}. The first gives the estimation risk of the particular method while the second calculates how long the solution took to compute.
#' @seealso empEstimRisk
#' @export
olsEstimAlgo <- function(job, static, dynamic, ...){
  ptm = proc.time()
  bhat = qr.solve(dynamic$Xtrain, dynamic$Ytrain)
  ptm = proc.time()-ptm
  estimRisk = empEstimRisk(bhat, dynamic$bstar)
  return(list(estimRisk = estimRisk, ptm=ptm[3]))
}

#' Calculates estimation error for ridge regression (optimal lambda)
#'
#' This function is used by \code{\link{BatchExperiments}} to run the simulations.
#' See that package for documentation.
#'
#' @inheritParams plsEstimAlgo
#'
#' @return A list with components \code{estimRisk} and \code{ptm}. The first gives the estimation risk of the particular method while the second calculates how long the solution took to compute.
#' @seealso empEstimRisk
#' @export
ridgeEstimAlgo <- function(job, static, dynamic,...){
  ptm = proc.time()
  ridgeOut = ridgeRegression(dynamic$Xtrain, dynamic$Ytrain,...)
  ptm = proc.time()-ptm
  estimRisk = empEstimRisk(ridgeOut$bhat, dynamic$bstar)
  return(list(estimRisk = estimRisk, lamStar=ridgeOut$lamStar, ptm=ptm[3]))
}

#' Calculates estimation error for ridge regression (optimal fixed lambda)
#'
#' This function is used by \code{\link{BatchExperiments}} to run the simulations.
#' See that package for documentation.
#'
#' @inheritParams plsEstimAlgo
#'
#' @return A list with components \code{estimRisk} and \code{ptm}. The first gives the estimation risk of the particular method while the second calculates how long the solution took to compute.
#' @seealso empEstimRisk
#' @export
ridgeEstimAlgoLam <- function(job, static, dynamic,...){
  ptm = proc.time()
  ridgeOut = ridgeRegressionLam(dynamic$Xtrain, dynamic$Ytrain,...)
  ptm = proc.time()-ptm
  estimRisk = empEstimRisk(ridgeOut, dynamic$bstar)
  return(list(estimRisk = estimRisk, ptm=ptm[3]))
}


#' Calculates prediction error for PLS
#'
#' This function is used by \code{\link{BatchExperiments}} to run the simulations.
#' See that package for documentation.
#'
#' @inheritParams plsEstimAlgo
#'
#' @return A list with components \code{risk}, \code{bias}, \code{var}, \code{noiseVar}, and \code{ptm}. These give, respectively, the prediction risk, the bias, the variance, the noise variance, and the computation time.
#' @seealso empPredRisk, predBias, predVar
plsAlgorithm <- function(job, static, dynamic, ...){
  ptm = proc.time()
  plsOut = pls(dynamic$Xtrain, dynamic$Ytrain, ...)
  ptm = proc.time()-ptm
  risk = empPredRisk(plsOut$bhat, dynamic$Xtest, dynamic$Ytest)
  sqbiasPred = predBias(dynamic$bstar, dynamic$Xtest, plsOut$XstuffTrain, plsOut$bhat)
  varPred = predVar(dynamic$sig, dynamic$Xtest, plsOut$XstuffTrain, plsOut$bhat)
  return(list(risk=risk, bias=sqbiasPred, var=varPred,  noiseVar = dynamic$sig^2, ptm=ptm[3]))
}

#' Calculates prediction error for OLS
#'
#' This function is used by \code{\link{BatchExperiments}} to run the simulations.
#' See that package for documentation.
#'
#' @inheritParams plsEstimAlgo
#'
#' @return A list with components \code{risk}, \code{bias}, \code{var}, \code{noiseVar}, and \code{ptm}. These give, respectively, the prediction risk, the bias, the variance, the noise variance, and the computation time.
#' @seealso empPredRisk, predBias, predVar
olsAlgorithm <- function(job, static, dynamic){
  ptm = proc.time()
  bhat = qr.solve(dynamic$Xtrain, dynamic$Ytrain)
  ptm = proc.time()-ptm
  Xstufffun <- function(X) solve(t(X) %*% X) %*% t(X)
  Xstuff = Xstufffun(dynamic$Xtrain)
  risk = empPredRisk(bhat, dynamic$Xtest, dynamic$Ytest)
  sqbiasPred = predBias(dynamic$bstar, dynamic$Xtest, Xstuff)
  varPred = predVar(dynamic$sig, dynamic$Xtest, Xstuff)
  return(list(risk=risk, bias=sqbiasPred, var=varPred, noiseVar = dynamic$sig^2, ptm=ptm[3]))
}
