
#' Calculates estimation error for PLS
#'
#' This function is used by \code{\link[=https://cran.r-project.org/web/packages/BatchExperiments/index.html]{BatchExperiments}} to run the simulations.
#' See that package for documentation.
#'
#' @param job see \code{BatchExperiments}
#' @param static see \code{BatchExperiments}
#' @param dynamic see \code{BatchExperiments}
#' @param ... see \code{BatchExperiments}
#'
#' @return A list with components \code{estimRisk}, \code{predRisk}, \code{baseline}, and \code{ptm}. The first gives the estimation risk of the particular method, the second is the MSE on the test set, the third is the MSE of the true coefficients, while the last calculates how long the solution took to compute.
#' @seealso empEstimRisk
#' @export
plsEstimAlgo <- function(job, static, dynamic, ...){
  ptm = proc.time()
  plsout = pls(dynamic$Xtrain, dynamic$Ytrain, ...)
  ptm = proc.time()-ptm
  estimRisk = empEstimRisk(plsout$bhat, dynamic$bstar)
  predRisk = mean((dynamic$Ytest - dynamic$Xtest %*% plsout$bhat)^2)
  baseline = mean((dynamic$Ytest - dynamic$Xtest %*% dynamic$bstar)^2)
  return(list(estimRisk = estimRisk, predRisk=predRisk,
              baseline=baseline, ptm=ptm[3]))
}

#' Calculates estimation error for OLS
#'
#' This function is used by \code{\link[=https://cran.r-project.org/web/packages/BatchExperiments/index.html]{BatchExperiments}} to run the
#' simulations. See that package for documentation.
#'
#' @inheritParams plsEstimAlgo
#'
#' @return A list with components \code{estimRisk}, \code{predRisk}, \code{baseline}, and \code{ptm}. The first gives the estimation risk of the particular method, the second is the MSE on the test set, the third is the MSE of the true coefficients, while the last calculates how long the solution took to compute.
#' @seealso empEstimRisk
#' @export
olsEstimAlgo <- function(job, static, dynamic, ...){
  ptm = proc.time()
  bhat = qr.solve(dynamic$Xtrain, dynamic$Ytrain)
  ptm = proc.time()-ptm
  estimRisk = empEstimRisk(bhat, dynamic$bstar)
  predRisk = mean((dynamic$Ytest - dynamic$Xtest %*% bhat)^2)
  baseline = mean((dynamic$Ytest - dynamic$Xtest %*% dynamic$bstar)^2)
  return(list(estimRisk = estimRisk, predRisk=predRisk,
              baseline=baseline,ptm=ptm[3]))
}

#' Calculates estimation error for ridge regression (optimal lambda)
#'
#' This function is used by \code{\link[=https://cran.r-project.org/web/packages/BatchExperiments/index.html]{BatchExperiments}} to run the simulations.
#' See that package for documentation.
#'
#' @inheritParams plsEstimAlgo
#'
#' @return A list with components \code{estimRisk}, \code{predRisk}, \code{baseline}, and \code{ptm}. The first gives the estimation risk of the particular method, the second is the MSE on the test set, the third is the MSE of the true coefficients, while the last calculates how long the solution took to compute.
#' @seealso empEstimRisk
#' @export
ridgeEstimAlgo <- function(job, static, dynamic,...){
  ptm = proc.time()
  ridgeOut = ridgeRegression(dynamic$Xtrain, dynamic$Ytrain,...)
  ptm = proc.time()-ptm
  estimRisk = empEstimRisk(ridgeOut$bhat, dynamic$bstar)
  predRisk = mean((dynamic$Ytest - dynamic$Xtest %*% ridgeOut$bhat)^2)
  baseline = mean((dynamic$Ytest - dynamic$Xtest %*% dynamic$bstar)^2)
  return(list(estimRisk = estimRisk, predRisk=predRisk,
              baseline=baseline, lamStar=ridgeOut$lamStar, ptm=ptm[3]))
}



#' Calculates estimation error for ridge regression (GCV-minimizing
#'  lambda)
#'
#' This function is used by \code{\link[=https://cran.r-project.org/web/packages/BatchExperiments/index.html]{BatchExperiments}} to run the simulations.
#' See that package for documentation.
#'
#' @inheritParams plsEstimAlgo
#'
#' @return A list with components \code{estimRisk}, \code{predRisk}, \code{baseline}, and \code{ptm}. The first gives the estimation risk of the particular method, the second is the MSE on the test set, the third is the MSE of the true coefficients, while the last calculates how long the solution took to compute.
#' @seealso empEstimRisk
#' @export
ridgeEstimAlgoLam <- function(job, static, dynamic,...){
  ptm = proc.time()
  ridgeOut = ridgeRegressionLam(dynamic$Xtrain, dynamic$Ytrain,...)
  ptm = proc.time()-ptm
  estimRisk = empEstimRisk(ridgeOut, dynamic$bstar)
  predRisk = mean((dynamic$Ytest - dynamic$Xtest %*% ridgeOut)^2)
  baseline = mean((dynamic$Ytest - dynamic$Xtest %*% dynamic$bstar)^2)
  return(list(estimRisk = estimRisk, predRisk=predRisk,
              baseline=baseline, ptm=ptm[3]))
}

#' Calculates degrees of freedom for compressed methods
#'
#' @inheritParams plsEstimAlgo
#' @return A list with components:
#' \describe{
#'  \item{\code{df}}{the (exact or approximate) degrees of freedom using the
#'  smoothing matrix}
#'  \item{\code{steinR}}{the SURE estimate of the risk using \code{df}}
#'  \item{\code{gcvR}}{the GCV estimate of the risk using \code{df}}
#'  \item{\code{divdf}}{the divergence-based degrees-of-freedom estimator for
#'  the combined methods}
#'  \item{\code{steindivR}}{the SURE estimate of the risk using \code{divdf}}
#'  \item{\code{gcvdivR}}{the GCV estimate of the risk using \code{divdf}}
#'  \item{\code{train}}{the training error}
#'  \item{\code{test}}{the testing error}
#'  }
#'
#' @export
plsRiskEstimAlgo <- function(job, static, dynamic,...){
  plsout = pls(dynamic$Xtrain, dynamic$Ytrain, ...)
  steinR = steinRisk(dynamic$Xtrain, dynamic$Ytrain, plsout$bhat, plsout$df)
  gcvR = gcvRisk(dynamic$Xtrain, dynamic$Ytrain, plsout$bhat, plsout$df)
  train = mean((dynamic$Ytrain - dynamic$Xtrain %*% plsout$bhat)^2)
  test = mean((dynamic$Ytest - dynamic$Xtest %*% plsout$bhat)^2)
  if(is.null(plsout$divdf)){
    steindivR = NA
    gcvdivR = NA
  }else{
    steindivR = steinRisk(dynamic$Xtrain, dynamic$Ytrain, plsout$bhat,
                           plsout$divdf)
    gcvdivR = gcvRisk(dynamic$Xtrain, dynamic$Ytrain, plsout$bhat,
                       plsout$divdf)
  }
  return(list(df=plsout$df, steinR=steinR, gcvR=gcvR,
              divdf=plsout$divdf, steindivR=steindivR, gcvdivR=gcvdivR,
              train=train, test=test))
}

#' Calculates degrees of freedom for ridge regression
#'
#' @inheritParams plsEstimAlgo
#' @return A list with components:
#' \describe{
#'  \item{\code{test}}{the testing error}
#'  \item{\code{train}}{the training error}
#'  \item{\code{df}}{the degrees of freedom using the
#'  smoothing matrix}
#'  \item{\code{gcvR}}{the GCV estimate of the risk using \code{df}}
#'  \item{\code{steinR}}{the SURE estimate of the risk using \code{df}}
#'  \item{\code{lamStar}}{the lambda value used, chosen via GCV}
#'  }
#'
#' @export
ridgePredAlgo <- function(job, static, dynamic,...){
  ridgeOut = ridgeRegression(dynamic$Xtrain, dynamic$Ytrain,...)
  test = mean((dynamic$Ytest - dynamic$Xtest %*% ridgeOut$bhat)^2)
  train = mean((dynamic$Ytrain - dynamic$Xtrain %*% ridgeOut$bhat)^2)
  df = sum(colSums(crossprod(dynamic$Xtrain) *
                  solve(crossprod(dynamic$Xtrain) + ridgeOut$lamStar *
                          diag(ncol(dynamic$Xtrain)))))
  gcvR = gcvRisk(dynamic$Xtrain, dynamic$Ytrain, ridgeOut$bhat, df)
  steinR = steinRisk(dynamic$Xtrain, dynamic$Ytrain, ridgeOut$bhat, df)
  return(list(test = test, train=train, df=df, gcvR=gcvR,
              steinR=steinR, lamStar=ridgeOut$lamStar))
}


#' Calculates degrees of freedom for OLS
#'
#' @inheritParams plsEstimAlgo
#' @return A list with components:
#' \describe{
#'  \item{\code{train}}{the training error}
#'  \item{\code{test}}{the testing error}
#'  \item{\code{df}}{the degrees of freedom using the
#'  smoothing matrix}
#'  \item{\code{steinR}}{the SURE estimate of the risk using \code{df}}
#'  \item{\code{gcvR}}{the GCV estimate of the risk using \code{df}}
#'  }
#'
#' @export
olsPredAlgo <- function(job, static, dynamic,...){
  ols = qr.solve(dynamic$Xtrain, dynamic$Ytrain)
  test = mean((dynamic$Ytest - dynamic$Xtest %*% ols)^2)
  train = mean((dynamic$Ytrain - dynamic$Xtrain %*% ols)^2)
  df = ncol(dynamic$Xtrain)
  gcvR = gcvRisk(dynamic$Xtrain, dynamic$Ytrain, ols, df)
  steinR = steinRisk(dynamic$Xtrain, dynamic$Ytrain, ols, df)
  return(list(test = test, train=train, df=df, gcvR=gcvR,
              steinR=steinR))
}
