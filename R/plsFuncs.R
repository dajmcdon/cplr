
#' Generate a preconditioner
#'
#' @param q number of rows
#' @param n number of columns
#' @param rand quoted name of a random number generator, defaults to \code{'rsparseBern'}
#' @param normalize optional scaler to normalize the matrix, defaults
#' to \eqn{\sqrt{q/3}}, which makes \eqn{E[Q'Q]=I} for sparse Bernoullis
#' with \code{s=3}, see \code{\link{rsparseBern}}
#' @param ... optional additional arguments passed to the random number generator
#'
#' @return A qxn matrix of random numbers
#' #@export
#'
#' @examples
#' generateQ(5, 10)
generateQ <- function(q, n, rand='rsparseBern', normalize = sqrt(q/3),...){
  rand = get(rand)
  Q = matrix(rand(n*q,...), nrow=q, ncol=n) / normalize
  return(Q)
}


#' Generate sparse Bernoulli random variables
#'
#' The sparse Bernoulli distribution puts mass on the set \eqn{\{-1, 0, 1\}} with probabilities \eqn{\{1/2s, 1-1/s, 1/2s\}}.
#'
#' @param n the number of random variates
#' @param s calibrates the probability of non-zeros. Default value is 3. Should be greater than one or results in unexpected behavior.
#'
#' @return Vector of random sparse Bernoulli draws
#' #@export
#'
#' @examples
#' rsparseBern(10)
#' rsparseBern(10, 10)
rsparseBern <- function(n, s=3){
  x=sample.int(3, size=n, replace=TRUE, prob=c(1/(2*s), 1-1/s, 1/(2*s)))
  x=x-2
  return(x)
}

#' Computes regularized preconditioned least squares
#'
#' This is the main function in this package. It computes different forms of compressed least squares regression.
#'
#' @param X the design matrix
#' @param Y the response vector
#' @param compression either full compression ('qxqy') or partial
#' compression ('qxy')
#' @param rand random number generator to build the compression matrix
#' @param q columns in the compression matrix
#' @param regTerm determines how the regularization is imposed, default is ridge-like (\code{'Ident'}). See details below.
#' @param Q pass along an existing Q matrix?
#' @param lam the amount of regularization
#' @param ... optional arguments to the random number generator
#'
#' @details This function implements two different compression algorithms with
#' different regularization versions. The notation here uses \eqn{||.||} to be the standard Euclidean norm. Full compression \code{compression='qxqy'}
#' is the solution to
#' \deqn{min||Q(Y-Xb)||^2 = min -Y'Q'QXb + b'X'Q'QXb = (X'Q'QX)^{-1}X'Q'QY.}
#' Partial compression simply ignores the \eqn{Q'Q} term between \eqn{X'} and \eqn{Y}: \deqn{min -Y'Q'QXb + b'X'Q'QXb = (X'Q'QX)^{-1}X'Y.}
#'
#' Either version can be regularized by adding a ridge penalty
#' \eqn{\lambda ||b||^2} to the optimization problem. However, one can iterpret
#' standard ridge regression as data augmentation, extending Y with p zeros and
#' similarly extending X with \eqn{\lambda I}. With this interpretation, we would
#' compress the extension as well either fully or partially. Therefore, the
#' \code{regTerm} can be \code{'Ident'} (standard case), \code{'QtQ'} (data
#' augmentation version), or \code{'zero'} (no regularization).
#'
#' @return a list with components \code{bhat} (the
#'  coefficient vector), \code{hatmat} (the smoothing
#'  matrix), and \code{df} (the degrees of freedom, trace of \code{hatmat}).
#' #@export
#'
#' @examples
#' n = 100
#' p = 5
#' q = 50
#' X = generateX(n, diag(1,p), 'rnorm')
#' Y = generateY(X, p:1, 'rnorm')
#' bhat = plsBasic(X, Y, 'qxqy', q=q, regTerm='Ident', lam=1)
plsBasic <- function(X, Y, compression=c('qxqy','qxy'), rand='rsparseBern',
                     q=NULL, regTerm=c('Ident'), Q=NULL, lam=0,...){
  p = ncol(X)
  n = length(Y)
  ## Build compressed matrices
  if(is.null(Q)) Q = generateQ(q, n, rand,...)
  reg = switch(as.character(regTerm),
               QtQ = {
                 Qs = generateQ(q, p, rand,...)
                 lam * crossprod(Qs)
               },
               Ident = lam * diag(p),
               zero = 0,
               stop('Invalid ridge input'))
  Xn = Q %*% X
  XX = crossprod(Xn) / n
  inv = switch(compression,
         qxqy = solve(XX + reg) %*% crossprod(Xn, Q) / n,
         qxy = solve(XX + reg) %*% t(X) /n,
         stop('Invalid compression type')
         )
  hatmat = X %*% inv
  df = sum(diag(hatmat))
  bhat = inv %*% Y
  return(list(bhat=bhat, hatmat=hatmat, df=df))
}

#' Computes regularized preconditioned least squares
#'
#' This function generalized \code{\link{plsBasic}} by combining full and partial
#' compression if desired.
#'
#' @param X the design matrix
#' @param Y the response vector
#' @param compression either full compression ('qxqy') or partial
#' compression ('qxy') or 'linComb' or 'convexComb'
#' @param q columns in the compression matrix
#' @param rand random number generator to build the compression matrix
#' @param regTerm determines how the regularization is imposed
#' @param SameQ use the same Q matirx in both parts of the linear combination
#' @param lam the amount of regularization
#' @param divergencedf If true, the Stein approximator to the degrees of
#'  freedom is calculated and returned.
#' @param ... optional arguments to the random number generator
#'
#' @details This function implements four different compression algorithms with
#' different regularization versions. In addition to full and partial compression
#' as in \code{\link{plsBasic}}, it also combines these two in two ways.
#'
#' The first method, 'linComb', estimates both full and partial compression, calculates their predictions, and then generates a weighted combination of the two. The weights are given as the solution to
#' \deqn{min ||[Yhat_{qxqy}, Yhat_{qxy}]a - Y||.}
#' This method can use the same Q or different Q matrices and any type of regularization.
#'
#' The second method, 'convexComb', is similar, though it first finds the fully and partially compressed versions using the same Q with the identity regularization. It then finds a convex combination of the fitted values. That is, it also solves \deqn{min ||[Yhat_{qxqy}, Yhat_{qxy}]a - Y||} subject to the constraint that \eqn{a=[a1, 1-a1]} with \eqn{a1 \in [0,1]}.
#'
#' @return A list with components:
#'  \describe{
#'    \item{\code{bhat}}{the vector of estimated coefficients}
#'    \item{\code{df}}{The degrees of freedom of the procedure. If full or partial
#'     compression, this is the trace of the smoothing matrix. For the other cases,
#'     the procedure is not a linear smoother, but rather a weighted sum of
#'     two linear smoothers. In that case, this value is simply the weighted sum
#'     of the two linear procedures. Note that this likely underestimates the
#'     true degrees of freedom.}
#'    \item{\code{hatmat}}{If we use full or partial compression, the procedure
#'     is linear in the response and the smoothing matrix is returned here.}
#'    \item{\code{divdf}}{If we use the same Q matrix for both procedures and use
#'     lambda times the identity as a regularizer, then Stein's Lemma gives an
#'     alternative degrees-of-freedom approximation. If requested,
#'     (\code{divergencedf==TRUE}), this estimate is also returned.}
#' }
#' #@export
#'
#' @seealso plsBasic
#'
#' @examples
#' n = 100
#' p = 5
#' q = 50
#' X = generateX(n, diag(1,p), 'rnorm')
#' Y = generateY(X, p:1, 'rnorm')
#' bhat = pls(X, Y, 'linComb', q=q, regTerm='Ident', lam=1)
pls <- function(X, Y, compression, q, rand, regTerm, SameQ, lam,
                divergencedf=FALSE,...){
  p = ncol(X)
  n = length(Y)

  Q = NULL
  ## if(grepl('Same',compression)) Q = generateQ(q, n, rand,...)
  if(SameQ) Q = generateQ(q, n, rand,...)
  switch(compression,
         qxqy = return(plsBasic(X, Y, rand=rand, compression='qxqy',q=q,
                           Q=Q, regTerm=regTerm, lam=lam)),
         qxy = return(plsBasic(X, Y, rand=rand, compression='qxy',q=q,
                           Q=Q, regTerm=regTerm, lam=lam)),
         linComb = {
           qxqy = plsBasic(X, Y, rand=rand, compression='qxqy',q=q,
                           regTerm=regTerm, Q = Q, lam=lam)
           qxy  = plsBasic(X, Y, rand=rand, compression='qxy',q=q,
                           regTerm=regTerm, Q = Q, lam=lam)
           XbhatQxy  = X%*%qxy$bhat
           XbhatQxqy = X%*%qxqy$bhat
           #Form linear combination
           alphaHat = qr.solve(cbind(XbhatQxy,XbhatQxqy),Y)
           bhat = qxy$bhat*alphaHat[1] + qxqy$bhat*alphaHat[2]
           df = qxy$df*alphaHat[1] + qxqy$df*alphaHat[2]
         },
         convexComb = {
           qxqy = plsBasic(X, Y, rand=rand, compression='qxqy',q=q,
                           regTerm=regTerm,Q=Q, lam=lam)
           qxy  = plsBasic(X, Y, rand=rand, compression='qxy',q=q,
                           regTerm=regTerm,Q=Q, lam=lam)
           XbhatQxy  = X%*%qxy$bhat
           XbhatQxqy = X%*%qxqy$bhat
           #Form linear combination
           A = cbind(XbhatQxy, XbhatQxqy)
           initialVals = .5
           ahat.optim = optim(initialVals, fn=objF, gr=objGradF,
                              A=A, Y=Y, method = 'L-BFGS-B',
                        lower=c(0), upper=c(1), hessian=FALSE)
           alphaHat = c(ahat.optim$par, 1-ahat.optim$par)
           bhat = qxy$bhat*alphaHat[1] + qxqy$bhat*alphaHat[2]
           df = qxy$df*alphaHat[1] + qxqy$df*alphaHat[2]
         },
         stop('Invalid compression type')
  )
  out = list(bhat=bhat, hatmat=NULL, df=df)
  if(divergencedf && SameQ && regTerm=='Ident' && compression=='linComb'){
    W = X %*% solve(crossprod(Q %*% X) + diag(lam, p)) %*% t(X)
    out$divdf = divergenceF(X, Y, Q, W, qxy$bhat, qxqy$bhat, lam)$divergence
  }
  return(out)
}


#' Fast version of pls
#'
#' This version uses a C function to perform the compression without ever
#' generating the Q matrix. There are fewer options that with \code{\link{pls}}
#' for this reason. The divergence-based risk estimator is unavailable.
#' An estimate of the degrees-of-freedom is returned.
#'
#' @param X the design matrix (n x p)
#' @param Y the response vector (n x 1)
#' @param compression the compression type (\code{'qxqy'}, \code{'qxy'}, \code{'linComb'}, \code{'convexComb'})
#' @param q the compression parameter
#' @param lam the ridge penalty
#' @param s optional sparsity parameter for the compression matrix, defaults
#'  to 3
#'
#' @return A list with components \code{bhat} (estimated coefficients) and
#'  \code{df} (the trace of the hat matrix, possibly weighted)
#' #@export
#'
#' @seealso pls
plsFast <- function(X, Y, compression, q, lam, s=3){
  p = ncol(X)
  n = length(Y)
  switch(compression,
         qxqy = {
           comp = compressCpp(X, q, Y, s)
           Xn = comp$QX
           Yn = comp$QY
           XX = crossprod(Xn) / n
           inv = solve(XX + lam*diag(p))
           df = sum(colSums(XX * inv)) # already divided by n
           bhat = inv %*% crossprod(Xn, Yn) / n
         },
         qxy = {
           comp = compressCpp(X, q, Y, s)
           Xn = comp$QX
           XX = crossprod(Xn) / n
           inv = solve(XX + lam*diag(p))
           df = sum(colSums(crossprod(X) * inv)) / n
           bhat = inv %*% crossprod(X, Y) / n
         },
         linComb = {
           comp = compressCpp(X, q, Y, s)
           Xn = comp$QX
           Yn = comp$QY
           XX = crossprod(Xn) / n
           inv = solve(XX + lam*diag(p))
           qxqy = inv %*% crossprod(Xn, Yn) / n
           qxy  = inv %*% crossprod(X, Y) / n
           XbhatQxy  = X%*%qxy
           XbhatQxqy = X%*%qxqy
           #Form linear combination
           alphaHat = solve_robust(cbind(XbhatQxy,XbhatQxqy),Y)
           bhat = qxy*alphaHat[1] + qxqy*alphaHat[2]
           dfqxqy = sum(colSums(XX * inv))
           dfqxy = sum(colSums(crossprod(X) * inv)) / n
           df = dfqxy*alphaHat[1] + dfqxqy*alphaHat[2]
         },
         convexComb = {
           comp = compressCpp(X, q, Y, s)
           Xn = comp$QX
           Yn = comp$QY
           XX = crossprod(Xn) / n
           inv = solve(XX + lam*diag(p))
           qxqy = inv %*% crossprod(Xn, Yn) / n
           qxy  = inv %*% crossprod(X, Y) / n
           XbhatQxy  = X%*%qxy
           XbhatQxqy = X%*%qxqy
           #Form linear combination
           A = XbhatQxy - XbhatQxqy
           b = Y - XbhatQxqy
           # initialVals = .5
           # ahat.optim = optim(initialVals, fn=objF, gr=objGradF,
           #                    A=A, b=b, method = 'L-BFGS-B',
           #                  lower=c(0), upper=c(1), hessian=FALSE)
           #alphaHat = c(ahat.optim$par, 1-ahat.optim$par)
           alphaHat = max(min(crossprod(A,b)/crossprod(A),1),0)
           alphaHat = c(alphaHat,1-alphaHat)
           bhat = qxy*alphaHat[1] + qxqy*alphaHat[2]
           dfqxqy = sum(colSums(XX * inv))
           dfqxy = sum(colSums(crossprod(X) * inv)) / n
           df = dfqxy*alphaHat[1] + dfqxqy*alphaHat[2]
         },
         stop('Invalid compression type')
  )
  out = list(bhat=bhat, df=df)
  return(out)
}

solve_robust <- function(A, b, tol=1e-10){
  ## assumes A has 2 columns
  AA = crossprod(A)
  detAA = abs(AA[1]*AA[4] - AA[2]*AA[3])
  if(detAA > tol) return(solve(AA, crossprod(A,b)))
  else return(c(.5,.5))
}

# objF <- function(alpha, A, b) sum(A*alpha - b)^2

# objGradF <- function(alpha, A, b) 2 * crossprod(A*alpha - b, A)



# invX.SVD <- function(X){
#   S = svd(X, nu=0)
#   VDinv2Vt = scale(S$v, center=FALSE, scale=S$d^2) %*% t(S$v)
#   return(VDinv2Vt)
# }


