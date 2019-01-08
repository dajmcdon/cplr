#' Create a sequence of numbers equally spaced on the log scale
#'
#' @param from the left endpoint
#' @param to the right endpoint
#' @param len the length of the resulting sequence, a non-negative number
#' @param integers if TRUE, rounds the result is rounded to the nearest integer
#'
#' @return a vector of length \code{len}
#'
#' @export
#' @examples
#' logScaleSeq(1,1000,20)
#' logScaleSeq(1,1000,20, FALSE)
logScaleSeq <- function(from=1, to=1, len=1, integers = TRUE){
  x = exp(seq(log(from), log(to), length.out=len))
  if(integers) x = round(x)
  return(x)
}

#' No intercept ridge regression, optimal lambda
#'
#' Computes the ridge regression estimate, optimizing over the tuning parameter
#' lambda using generalized cross validation (GCV)
#'
#' @param X design matrix
#' @param Y response vector
#' @param lam.upper a vector with the smallest and largest lambda values to try
#'
#' @return A list containing th vector of coefficients at the GCV minimizing value of lambda (\code{bhat}) and the best value of lambda (\code{lamStar}).
#' @export
ridgeRegression <- function(X, Y, lam.upper = min(dim(X))){
  p = ncol(X)
  n = nrow(X)
  XX = crossprod(X) / n
  XY = crossprod(X, Y) / n

  lamStar = optimize(GCVridge, interval=c(0,lam.upper), XX=XX, XY=XY, X=X, Y=Y)$minimum
  bhat = solve(XX + lamStar*diag(p), XY)

  return(list(bhat=bhat, lamStar = lamStar))
}



#' The generalized cross validation (GCV) score for ridge regression
#'
#' This function is designed to be optimized over using already computed
#' information. It is not intended to stand alone.
#'
#' @param lam The tuning parameter
#' @param XX The inner product of the design with itself (a pxp matrix)
#' @param XY The inner product of the design with the response (a px1 vector)
#' @param X The design matrix
#' @param Y The response vector
#'
#' @return The generalized cross validation score.
GCVridge <- function(lam, XX, XY, X, Y){
  p = ncol(X)
  n = nrow(X)
  hmat = X %*% solve(XX + lam* diag(p))
  num = sum((Y-hmat %*% XY)^2)
  denom = sum(1-rowSums(hmat*X)/n)^2
  return(num/denom)
}
