hF = function(Y,X,QtQ,W){
  h_1 = t(Y) %*% W %*% Y
  h_2 = t(Y) %*% QtQ %*% W %*% Y
  partial_h_1.partial_yi = 2 * W %*% Y
  partial_h_2.partial_yi = (QtQ %*% W + W %*% QtQ) %*% Y
  return(list('h'         = c(h_1,h_2),
              'partial_h' = cbind(partial_h_1.partial_yi,partial_h_2.partial_yi)))
  # partial h_j/partial y_i = partial_h[i,j],  j = 1,2
}

###g function stuff
###
### AtA = | a b |  So, (AtA)^(-1) = | d -b |
###       | c d |                   |-c  a |
innerProdForGfuncF = function(b_1,XtX,b_2=NULL){
  if(is.null(b_2)) b_2 = b_1
  return(drop(t(b_1) %*% XtX %*% b_2))
}
determinant_yF = function(a,b,c,d){
  return(a*d - b*c)
}

gOnlyF = function(b_PC,b_FC,X){
  XtX = t(X) %*% X
  a = innerProdForGfuncF(b_PC,XtX)
  b = innerProdForGfuncF(b_PC,XtX,b_FC)
  c = b
  d = innerProdForGfuncF(b_FC,XtX)
  det_y = determinant_yF(a,b,c,d)
  if(abs(det_y) < 1e-16) cat('zero determinant \n')
  g = matrix(c(d,-b,-c,a)/det_y,nrow=2,ncol=2,byrow=TRUE)
  return(list('g'=g,'a'=a,'b'=b,'c'=c,'d'=d,'det_y'=det_y))
}

gDerivativeF = function(Y,X,W,QtQ,a,b,c,d,det_y){
  n = length(Y)
  ##### Derivatives:
  Wsq    = W %*% W
  WsqQtQ = Wsq %*% QtQ
  QtQWsq = QtQ %*% Wsq
  Grad_a = 2* Wsq %*% Y
  Grad_b = (QtQWsq + t(QtQWsq)) %*% Y
  Grad_c = Grad_b
  Grad_d = 2* QtQWsq %*% QtQ %*% Y
  Grad_det_y = Grad_a * d + a * Grad_d - Grad_b * c - b * Grad_c
  #To get each component function:
  #Example for g_11:
  # partial g_11 partial y_i = partial [d(y)/det(y)] partial y_i
  # = partial (d(y) det(y) - d(y) partial det(y))/det(y)**2
  partial_g = array(0,dim=c(2,2,n))
  partial_g[1,1,] = Grad_d * det_y - d * Grad_det_y
  partial_g[2,2,] = Grad_a * det_y - a * Grad_det_y
  partial_g[1,2,] = - Grad_b * det_y + b * Grad_det_y
  partial_g[2,1,] = - Grad_c * det_y + c * Grad_det_y
  partial_g = partial_g/det_y**2
  return(partial_g)
}



#' Divergence calculation for degrees of freedom
#'
#' For approximating the degrees of freedom of the convex combination of full
#' and partial compression.
#'
#' @param X the \eqn{n \times p} design matrix
#' @param Y the response vector
#' @param Q the compression matris
#' @param W for faster computation, \eqn{W = X (X'Q'QX + lam I)^{-1} X'}
#' @param b_PC the coefficient vector estimated via partial compression
#' @param b_FC the coefficient vector estimated via full compression
#' @param lam the tuning parameter
#'
#' @return A list with two components: \code{divergence} gives the estimated
#' degrees of freedom, \code{partial_F.partial.y} gives the gradient with respect
#' to the response \code{Y}.
#'
#' @details
#'
#' \deqn{W = X( t(X)t(Q)QX + lam I)^(-1)t(X)}
#' \deqn{\hat{Y} = F(y) = W f(y) g(y) h(y)}
#' \deqn{\Rightarrow F_i(y) = W[i,] f(y) g(y) h(y)}
#' \deqn{\partial F_i(y) / \partial y_i = W[i,] (df gh + f(dg h + g dh))}
#'
#' @export
divergenceF = function(X,Y,Q,W,b_PC,b_FC,lam){
  n      = nrow(X);p = ncol(X)
  QtQ    = t(Q) %*% Q
  ## Compute partial Omega/partial y_i= partial Omega_1/partial y_i + partial Omega_2/partial y_i
  hStuff               = hF(Y,X,QtQ,W)
  h                    = hStuff$h
  partial_h.partial.yi = hStuff$partial_h
  gOnlyStuff           = gOnlyF(b_PC,b_FC,X)
  g = gOnlyStuff$g;  a = gOnlyStuff$a;b = gOnlyStuff$b;d = gOnlyStuff$d;det_y = gOnlyStuff$det_y
  c = b#symmetric matrix
  #OmegaBar_1 = matrix(0,nrow=n,ncol=n)
  #OmegaBar_2 = matrix(0,nrow=n,ncol=n)
  #f function stuff
  f_1 = Y
  f_2 = QtQ %*% Y
  f   = cbind(f_1,f_2)
  #Delta stuff
  Delta_1 = f[,1]*g[1,1] + f[,2]*g[2,1]
  Delta_2 = f[,1]*g[1,2] + f[,2]*g[2,2]
  #g function stuff
  partial_g = gDerivativeF(Y,X,W,QtQ,a,b,c,d,det_y)
  OmegaBar_1 = outer(partial_h.partial.yi[,1],Delta_1) +
    h[1]*(diag(g[1,1],n) +
            outer(f[,1],partial_g[1,1,]) +
            outer(f[,2],partial_g[2,1,]) +
            QtQ*g[2,1])
  OmegaBar_2 = outer(partial_h.partial.yi[,2],Delta_2) +
    h[2]*(diag(g[1,2],n) +
            outer(f[,1],partial_g[1,2,]) +
            outer(f[,2],partial_g[2,2,]) +
            QtQ*g[2,2])
  OmegaBar   = OmegaBar_1 + OmegaBar_2
  partial_F.partial.y = rowSums(W * t(OmegaBar))
  return(list(divergence = sum(partial_F.partial.y),
              partial_F.partial.y = partial_F.partial.y))
}
