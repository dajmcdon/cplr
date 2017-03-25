#' Generalized cross validation risk estimate
#'
#' @inheritParams steinRisk
#'
#' @details The calculation is given by
#'  \deqn{\sum (X \hat{beta} - Y)^2 / (1 - df/n)^2}
#'
#' @return The risk estimate.
#' @export
gcvRisk <- function(X, Y, bhat, df){
  n = length(Y)
  mse = mean((X %*% bhat - Y)^2)
  denom = (1 - df/n)^2
  return(mse/denom)
}

#' Computes Stein's Risk Estimate
#'
#' @param X the design matrix
#' @param Y the response vector
#' @param bhat the estimate of the regression coefficients
#' @param df the estimated (or exact) degrees of freedom
#'
#' @return The risk estimate. Note that this is not necessarily unbiased as the variance is unknown
#'  and the estimator used here is not necessarily unbiased.
#' @export
steinRisk <- function(X, Y, bhat, df){
  n = length(Y)
  sse = sum((X %*% bhat - Y)^2)
  sig2hat = sse / (n - df)
  risk = sse - n*sig2hat + 2*sig2hat*df
  return(risk)
}
