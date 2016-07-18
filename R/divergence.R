###h function stuff
HmatF = function(X,QtQ,lam){
  p = ncol(X)
  return(X %*% solve( t(X)%*%t(Q)%*%Q%*%X + lam*diag(p)) %*% t(X))
}
hF = function(Y,X,QtQ,H){
  h_1 = t(Y) %*% H %*% Y
  h_2 = t(Y) %*% QtQ %*% H %*% Y
  partial_h_1.partial_yi = 2 * H %*% Y
  partial_h_2.partial_yi = 2 * QtQ %*% H %*% Y
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
  return(t(b_1) %*% XtX %*% b_2)
}
determinant_yF = function(a,b,c,d){
  return(a*d - b*c)
}

gOnlyF = function(b_PP,b_FP,X){
  XtX = t(X) %*% X
  a = innerProdForGfuncF(b_PP,XtX)
  b = innerProdForGfuncF(b_PP,XtX,b_FP)
  c = b
  d = innerProdForGfuncF(b_FP,XtX)
  det_y = determinant_yF(a,b,c,d)
  if(abs(det_y) < 1e-16) cat('zero determinant \n')
  g = matrix(c(d,-b,-c,a)/det_y,nrow=2,ncol=2,byrow=TRUE)
  return(list('g'=g,'a'=a,'b'=b,'c'=c,'d'=d,'det_y'=det_y))
}

gDerivativeF = function(b_PP,b_FP,X,H,QtQ,i,g,a,b,c,d,det_y){
  ##### Derivatives:
  H_i    = H[,i]
  HQtQ   = H %*% QtQ
  HQtQ_i = HQtQ[,i]
  partial_a.partial_yi = 2* t(H_i) %*% H %*% Y
  partial_b.partial_yi = 2 * t(H_i) %*% HQtQ %*% Y
  partial_c.partial_yi = 2 * t(HQtQ_i) %*% H %*% Y
  partial_d.partial_yi = 2 * t(HQtQ_i) %*% HQtQ %*% Y
  partial_det_y.partial_yi = partial_a.partial_yi * d +
    a * partial_d.partial_yi -
    partial_b.partial_yi * c -
    b * partial_c.partial_yi
  #To get each component function:
  #Example for g_11:
  # partial g_11 partial y_i = partial [d(y)/det(y)] partial y_i
  # = partial (d(y) det(y) - d(y) partial det(y))/det(y)**2
  partial_g_11.partial_yi = (partial_d.partial_yi * det_y - d * partial_det_y.partial_yi)/det_y**2
  partial_g_12.partial_yi = -(partial_b.partial_yi * det_y - b * partial_det_y.partial_yi)/det_y**2
  partial_g_21.partial_yi = -(partial_c.partial_yi * det_y - c * partial_det_y.partial_yi)/det_y**2
  partial_g_22.partial_yi = -(partial_a.partial_yi * det_y - a * partial_det_y.partial_yi)/det_y**2
  partial_g.partial_yi = matrix(c(partial_g_11.partial_yi,partial_g_12.partial_yi,
                                  partial_g_21.partial_yi,partial_g_22.partial_yi),
                                nrow=2,byrow=TRUE)
  return(list('partial_g' = partial_g.partial_yi))
}
###f function stuff
###
### f = [Y,QtQ Y] \Rightarrow \partial f \partial y_i = [e_i,t(Q) %*% Q[,i]]
###   where e_i is i^{th} canonical basis vector
fF = function(Y,QtQ,i){
  e_i = rep(0,length(Y));e_i[i] = 1
  f_1 = Y
  f_2 = QtQ %*% Y
  partial_f_1.partial_yi = e_i
  partial_f_2.partial_yi = QtQ[,i]
  return(list('f' = cbind(f_1,f_2),
              'partial_f' = cbind(partial_f_1.partial_yi,partial_f_2.partial_yi)))
  # partial f_j/partial y_i = partial_f[i,j],  j = 1,2
}


#' Derivative of F
#'
#' For approximating the degrees of freedom of the convex combination of full
#' and partial compression.
#'
#' @param X the \eqn{n \times p} design matrix
#' @param Y the response vector
#' @param QtQ the crossproduct of the compression matris
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
divergenceF = function(X,Y,QtQ,b_PP,b_FP,lam){
  n      = nrow(X);p = ncol(X)
  W      = X %*% solve( t(X) %*% QtQ %*% X + diag(lam,p) ) %*% t(X)
  QtQ    = t(Q) %*% Q
  H      = HmatF(X,QtQ,lam)
  ## Compute partial Omega/partial y_i= partial Omega_1/partial y_i + partial Omega_2/partial y_i
  hStuff               = hF(Y,X,QtQ,H)
  h                    = hStuff$h
  partial_h.partial.yi = hStuff$partial_h
  gOnlyStuff           = gOnlyF(b_PP,b_FP,X)
  g = gOnlyStuff$g;a = gOnlyStuff$a;b = gOnlyStuff$b;d = gOnlyStuff$d;det_y = gOnlyStuff$det_y
  c = b#symmetric matrix
  Omega_1 = matrix(0,nrow=n,ncol=1)
  Omega_2 = matrix(0,nrow=n,ncol=1)
  partial_F.partial.y = rep(0,n)
  for(i in 1:n){
    fStuff = fF(Y,QtQ,i)
    gDerivativeStuff     = gDerivativeF(b_PP,b_FP,X,H,QtQ,i,g,a,b,c,d,det_y)
    f                    = fStuff$f
    partial_f.partial.yi = fStuff$partial_f
    partial_g.partial.yi = gDerivativeStuff$partial_g
    for(k in 1:n){
      Delta_k1 = f[k,1]*g[1,1] + f[k,2]*g[2,1]
      Delta_k2 = f[k,1]*g[1,2] + f[k,2]*g[2,2]
      Omega_k1 = partial_h.partial.yi[i,1]*Delta_k1 +
        h[1]*(partial_f.partial.yi[k,1]*g[1,1] + f[k,1]*partial_g.partial.yi[1,1] +
                partial_f.partial.yi[k,2]*g[2,1] + f[k,2]*partial_g.partial.yi[2,1])
      Omega_k2 = partial_h.partial.yi[i,2]*Delta_k2 +
        h[2]*(partial_f.partial.yi[k,1]*g[1,2] + f[k,1]*partial_g.partial.yi[1,2] +
                partial_f.partial.yi[k,2]*g[2,2] + f[k,2]*partial_g.partial.yi[2,2])
      Omega_1[k] = Omega_k1
      Omega_2[k] = Omega_k2
    }
    Omega = Omega_1 + Omega_2
    partial_Fi.partial.yi  = W[i,] %*% Omega
    partial_F.partial.y[i] = partial_Fi.partial.yi
  }
  return(list(divergence = sum(partial_F.partial.y),
              partial_F.partial.y = partial_F.partial.y))
}
