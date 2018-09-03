#' Calculate the MSE
#'
#' @param bhat estimated coefficient vector
#' @param bstar true coefficient vector
#'
#' @return The average squared estimation error
#' \deqn{mean((bhat - bstar)^2).}
empEstimRisk <- function(bhat, bstar) mean((bhat-bstar)^2)

#' Calculate the MSPE on new data
#'
#' @param bhat coefficient vector
#' @param newX new design matrix
#' @param newY new response vector
#'
#' @return The average squared prediction error
#' \deqn{mean((newY-newX bhat)^2).}
empPredRisk <- function(bhat, newX, newY){
  return(mean((newY - newX %*% bhat)^2))
}

#' Calculate the average squared bias (per prediction)
#'
#' @param bstar true coefficient vector
#' @param newX new design matrix
#' @param XstuffTrain the portion of the estimated coefficient vector that depends
#' on the design matrix from the training set or \code{NULL}.
#' @param bhat estimated coefficient vector
#'
#' @details If \code{XstuffTrain} is given, the linear smoother gives the
#'  \eqn{\E[\hat{Y}]}. Otherwise (NULL), evaluate the predictions on the
#'  test data to get monte carlo bias
#'
#' @return The average squared difference between the predictions and their
#' expectation
predBias <- function(bstar, newX, XstuffTrain, bhat=NULL){
  expectY = newX %*% bstar
  if(!is.null(XstuffTrain)){
    ## this is a linear predictor, calculate exactly using hat matrix
    predY = newX %*% XstuffTrain %*% expectY
  } else {
    ## not a linear predictor, evaluate on test data to get monte carlo bias
    predY = newX %*% bhat
  }
  return(mean((predY-expectY)^2))
}

#' Calculate the average variance (per prediction)
#'
#' @param sigma the standard deviation of the noise
#' @param newX new design matrix
#' @param XstuffTrain the portion of the estimated coefficient vector that depends
#' on the design matrix from the training set or \code{NULL}.
#' @param bhat estimated coefficient vector or \code{NULL}
#'
#' @details If \code{XstuffTrain} is given, the linear smoother gives the
#'  variance exactly as the trace of the squared "hat" matrix (not necessarily
#'  idempotent). Otherwise (NULL), evaluate the predictions on the
#'  test data to get monte carlo variance
#'
#' @return The average variance per prediction
predVar <- function(sigma, newX, XstuffTrain, bhat=NULL){
  if(!is.null(XstuffTrain)){
    ## this is a linear predictor, calculate exactly using hat matrix
    Hmat = newX %*% XstuffTrain
    outVar = sigma^2 * mean(rowSums(H^2))
  } else {
    ## not a linear predictor, evaluate on test data to get monte carlo variance estimate
    ## easier, but not necessarily the best (could use all replications to get est. of the mean
    preds = newX %*% bhat
    outVar = mean((preds-mean(preds))^2)
  }
  return(outVar)
}
