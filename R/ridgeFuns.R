pls_scale <- function(X, y, n=nrow(X), p=ncol(X)){
  Xm = colMeans(X)
  X = scale(X, center=Xm, scale=FALSE)
  Xscale = sqrt(colSums(X^2))
  X = scale(X,center=FALSE,scale=Xscale)
  ym = mean(y)
  y = matrix(y-ym,ncol=1)
  return(list(Xs=X, Xm=Xm, Xscale=Xscale, ys=y, ym=ym))
}

ridge_svd <- function(S, scaled,
                      type = c('xy','qxqy','qxy'), comp,
                      lam=NULL, lam.max=NULL, lam.min=1e-6, nlam=100,
                      tol.lam0=1e-10){
  n = nrow(scaled$Xs)
  p = ncol(scaled$Xs)
  type = match.arg(type)
  Xs = scaled$Xs
  ys = scaled$ys
  Xm = scaled$Xm
  Xscale = scaled$Xscale
  ym = scaled$ym

  U = S$u
  d = S$d # / sqrt(n)
  dx = length(d)
  V = S$v
  if(type=='qxy') XV = Xs %*% V
  rhs = switch(type,
               'xy' = crossprod(U, ys), #/ n,
               'qxqy' = crossprod(U, comp$QY), #/ n,
               'qxy' = crossprod(XV, ys), # / n,
               stop('Incorrect type argument')
  )
  lam = fix_lam(d, p, lam, lam.max, lam.min, nlam, rhs, tol.lam0)

  k = length(lam)
  div = d^2 + rep(lam, rep(dx, k))
  if(type=='qxy') a = drop(rhs)/div
  else a = drop(d * rhs)/div
  dim(a) = c(dx, k)
  bhatsc = V %*% a
  mat = matrix(d^2/div, dx)
  df = switch(type,
              'xy' = colSums(mat) + 1,
              'qxqy' = colSums(mat) + 1,
              'qxy' = colSums(XV^2) %*% matrix(1/div, dx) + 1
  )
  fitted = scaled$Xs %*% bhatsc
  residuals = c(ys) - scaled$Xs %*% bhatsc
  train = colMeans( (residuals)^2 )
  GCV = train / (1-df/n)^2
  bhat = bhatsc/Xscale
  bhat0 = ym - colSums(bhat*Xm)
  out = list(intercept=bhat0, bhat=bhat, bhatsc=bhatsc,
             fitted=fitted, residuals=residuals, GCV=GCV, lam=lam, df=df,
             train=train)
  class(out) = 'cplr'
  return(out)
}



#' Perform compressed (penalized) linear regression
#'
#' @param X the design matrix \code{(n x p)}
#' @param Y the response vector
#' @param compression either none ('xy'), full compression ('qxqy'), partial
#' compression ('qxy'), linear combination ('linComb') or  the convex combination ('convexComb')
#' @param q columns in the compression matrix
#' @param lam optional values of the tuning parameter. Default is \code{[lam.min, lam.max]} with
#' \code{n}lam entries equally spaced on the log scale
#' @param lam.max defaults to the maximum L1-norm of the covariates divided by 1e-3
#' @param lam.min defaults to 1e-6, unless the smallest singular value is larger than \code{tol},
#' in which case 0 is included
#' @param nlam number of lambda values, default is 100, 101 if 0 is included
#' @param s \code{1/s} is the probability of non-zero entries in the compression matrix, default
#' is 3
#' @param tol.lam0 determines how close to singular the design can be and still try to
#' use \code{lam=0}. Ignored if \code{lam} is given. Compares with the smallest singular value of
#' the design.
#' @param tol.lc determines how close to singular the matrix of fitted values can be. Compares with
#' the smallest singular value. Default is then a 50-50 combination of \code{qxy} and \code{qxqy}.
#'
#' @return A list with components of class 'cplr':
#'  \describe{
#'  \item{\code{intercept}}{a vector of length \code{nlam} containing the intercept}
#'    \item{\code{bhat}}{a matrix of size \code{p x nlam} containing the estimated coefficients}
#'    \item{\code{bhatsc}}{a matrix of size \code{p x nlam} containing the scaled,
#'    estimated coefficients. For use with the \code{plot} method.}
#'    \item{\code{fitted}}{a matrix of size \code{n x nlam} containing fitted values}
#'    \item{\code{residuals}}{a matrix of size \code{n x nlam} containing residuals}
#'    \item{\code{GCV}}{The generalized cross validation score for model selection.}
#'    \item{\code{lam}}{The sequence of \code{lam} values used. Either generated or user supplied.}
#'    \item{\code{df}}{The degrees of freedom of the procedure. If full or partial
#'     compression, this is the trace of the smoothing matrix. For the other cases,
#'     the procedure is not a linear smoother, but rather a weighted sum of
#'     two linear smoothers. In that case, this value is simply the weighted sum
#'     of the two linear procedures. Note that this likely underestimates the
#'     true degrees of freedom.}
#'     \item{\code{train}}{The training error for each value of \code{lam}}
#' }
#' @export
#'
#' @examples
#' n = 100
#' p = 5
#' q = 50
#' X = generateX(n, diag(1,p), 'rnorm')
#' Y = generateY(X, p:1, 'rnorm')
#' out = compressedRidge(X, Y, 'linComb', q=q, lam.max=1)
compressedRidge <- function(X, Y, compression = c('xy','qxqy','qxy','linComb','convexComb'),
                            q, lam=NULL, lam.max = NULL, lam.min=1e-6, nlam=100, s=3,
                            tol.lam0=1e-8, tol.lc=1e-10){
  type = match.arg(compression)
  p = ncol(X)
  n = length(Y)
  scaled = pls_scale(X, Y, n, p) # scale before compression
  if(type != 'xy'){
    comp = compressCpp(scaled$Xs, q, scaled$ys, s)
    S = svd(comp$QX)
  }
  switch(type,
         xy = {
           S = svd(scaled$Xs)
           return(ridge_svd(S, scaled, type='xy', NULL, lam, lam.max, lam.min, nlam, tol.lam0))
         },
         qxqy = {
           return(ridge_svd(S, scaled, type='qxqy',comp, lam, lam.max, lam.min, nlam, tol.lam0))
         },
         qxy = {
           return(ridge_svd(S, scaled, type='qxy',comp, lam, lam.max, lam.min, nlam, tol.lam0))
         },
         linComb = {
           return(compress_Comb(S, scaled, comp, lam, lam.max, lam.min, nlam, ahat_linComb,
                                tol.lam0, tol.lc))
         },
         convexComb = {
           return(compress_Comb(S, scaled, comp, lam, lam.max, lam.min, nlam, ahat_convexComb,
                                tol.lam0, tol.lc))
         }
           )
}

fix_lam <- function(d, p, lam, lam.max, lam.min, nlam, rhs, tol.lam0 = 1e-10){
  if(is.null(lam)){
    dx = length(d)
    if(is.null(lam.max)) lam.max = max(abs(rhs)) / 1e-3 ## see ?glmnet::glmnet
    lam = logScaleSeq(lam.min,lam.max,nlam,FALSE)
    if(dx == p && d[dx] > tol.lam0) lam = c(0, lam)
  }else{
    lam = sort(lam)
  }
  return(lam)
}

ahat_linComb <- function(Xs, qxy, qxqy, Ys, nlam, tol){
  Yhat = rbind(Xs %*% qxy$bhatsc, Xs %*% qxqy$bhatsc)
  dim(Yhat) = c(nrow(Xs),2,nlam)
  ahat = apply(Yhat, 3, function(Yh){ # a 2 x nlam matrix
    fit = lm.fit(Yh, Ys)
    if(svd(qr.R(fit$qr),0,0)$d[2]<tol){
      return(c(.5,.5))
    }else{
      return(fit$coefficients)
    }
  }
  )
  aqxy = matrix(ahat[1,], nrow=ncol(Xs), ncol=nlam, byrow=TRUE)
  aqxqy = matrix(ahat[2,], nrow=ncol(Xs), ncol=nlam, byrow=TRUE)
  return(list(aqxy = aqxy, aqxqy = aqxqy))
}

ahat_convexComb <- function(Xs, qxy, qxqy, Ys, nlam, tol){
  yhat1 = Xs %*% qxy$bhatsc
  yhat2 = Xs %*% qxqy$bhatsc
  A = yhat1 - yhat2
  b = drop(Ys) - yhat2
  AtA = colSums(A^2)
  Atb = colSums(A*b)
  ahat = Atb / AtA
  ahat[ahat < 0] = 0
  ahat[ahat > 1] = 1
  aqxy = matrix(ahat, nrow=ncol(Xs),ncol=nlam,byrow=TRUE)
  aqxqy = 1- aqxy
  return(list(aqxy = aqxy, aqxqy = aqxqy))
}

compress_Comb <- function(S, scaled, comp, lam, lam.max, lam.min, nlam, ahatfun, tol.lam0, tol.lc){
  n = length(scaled$ys)
  qxy = ridge_svd(S, scaled, 'qxy', comp, lam, lam.max, lam.min, nlam, tol.lam0)
  if(is.null(lam)) lam=qxy$lam
  qxqy = ridge_svd(S, scaled, 'qxqy', comp, lam)
  ahat = ahatfun(scaled$Xs, qxy, qxqy, scaled$ys, length(lam), tol.lc)
  bhatsc = qxy$bhatsc * ahat$aqxy + qxqy$bhatsc * ahat$aqxqy
  bhat = bhatsc/scaled$Xscale
  bhat0 = scaled$ym - colSums(bhat*scaled$Xm)
  df = pmin(drop(qxy$df * ahat$aqxy[1,] + qxqy$df * ahat$aqxqy[1,]),n)
  fitted = scaled$Xs %*% bhatsc
  residuals = c(scaled$ys) - scaled$Xs %*% bhatsc
  train = colMeans( (residuals)^2 )
  GCV = train / (1-df/n)^2
  out = list(intercept=bhat0, bhat=bhat, bhatsc = bhatsc,
             fitted=fitted, residuals=residuals, GCV=GCV, lam=lam, df=df,
             train=train)
  class(out) = 'cplr'
  return(out)
}