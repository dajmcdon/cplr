divergence_LC <- function(ahat, S, fitted, df, qxy, bqxy, qxqy, bqxqy, tol){

  # Lambdas with ahat_qxy=qhat_qxqy=.5 are easy
  ahat_p = ahat$aqxy[1,]
  ahat_f = ahat$aqxqy[1,]
  easy = (abs(ahat_p - ahat_f) < tol) # for which lambdas
  diverg = df
  if(all(easy)) return(diverg)
  n = nrow(qxy$XR)
  hard = (1:length(diverg))[!easy]

  L2 = qxqy$L2
  L2divLSQY = L2 * qxqy$divLSQY
  L2divRXY = L2 * qxy$divRXY

  dim(qxy$div) = dim(qxqy$divLSQY)

  for(i in hard){
    # tr(PzHp)
    GWA = svd(drop(ahat$yh[,,i]))
    GtXR = crossprod(GWA$u,qxy$XR)
    PZHPC =  sum(colSums(GtXR^2) / qxy$div[,i])
    # tr(PzHf)
    AWinv = sweep(GWA$v, 2, GWA$d, '/')
    B = cbind(bqxy[,i], bqxqy[,i])
    PZHFC = sum(colSums(S$v * B %*% AWinv %*% GtXR) * qxy$L2div[,i])
    # the rests
    Zinv = tcrossprod(GWA$u, AWinv)
    XRdiv = sweep(qxy$XR ,2, qxy$div[,i], '/')
    HpYhat = XRdiv %*% crossprod(qxy$XR,fitted[,i])
    HfYf = XRdiv %*% L2divLSQY[,i]
    HfYp = XRdiv %*% L2divRXY[,i]
    HHYhat = cbind(HfYp, HfYf)
    diverg[i] = diverg[i] + 2 - ahat_p[i]*PZHPC - ahat_f[i]*PZHFC -
      sum(colSums(Zinv * HHYhat))
  }
  return(diverg)
}


divergence_CC <- function(ahat, fitted, df, qxy, qxqy, tol){
  # Lambdas with ahat_qxy=qhat_qxqy=.5 are easy
  ahat_p = ahat$aqxy[1,]
  ahat_f = ahat$aqxqy[1,]
  easy = (ahat_f < tol) | (ahat_p < tol)
      # for which lambdas
  if(all(easy)) return(df)
  yh1 = qxy$fitted
  yh2 = qxqy$fitted
  #dim(yh) = c(nrow(qxy$fitted), 2, ncol(qxy$fitted))
  hard = (1:length(df))[!easy]

  # Precomputations
  L2 = qxqy$L2
  div = qxy$div
  L2divLSQY = L2 * qxqy$divLSQY
  L2divRXY = L2 * qxy$divRXY
  divL2divLSQY = L2divLSQY / div
  divL2divRXY = L2divRXY / div
  RXXR = crossprod(qxy$XR) # p x p
  RXXRdivRXY = RXXR %*% qxy$divRXY
  RXXRdivLSQY = RXXR %*% qxqy$divLSQY
  # Term 1
  YHpYp = colSums(qxy$divRXY * RXXRdivRXY) # YXRdivRXXRdivRXY
  YHpYf = colSums(qxy$divRXY * RXXRdivLSQY)  # YXRdivRXXRdivLSQY
  YfHpYp = colSums(RXXRdivLSQY * RXXRdivRXY / div) # YQSLdivRXXRdivRXXRdivRXY
  YfHpYf = colSums(RXXRdivLSQY^2 / div) # YQSLdivRXXRdivRXXRdivLSQY
  Term1 = YHpYp - YHpYf - YfHpYp + YfHpYf
  # Term 2
  YpYp = colSums(yh1^2)
  YfYp = colSums(yh1*yh2)
  YpHfYp = colSums(RXXRdivRXY * divL2divRXY) # YXRdivRXXRdivL2divRXY
  YpHfYf = colSums(RXXRdivRXY * divL2divLSQY) # YXRdivRXXRdivL2divLSQY
  Term2 = YpYp - YfYp - YpHfYp + YpHfYf
  # Term 3
  YfHfYf = colSums(RXXRdivLSQY * L2divLSQY / div) # YQSLdivRXXRdivL2divLSQY
  YpHpYp = colSums(RXXRdivRXY^2 / div) # YXRdivRXXRdivRXXRdivRXY
  Term3 = YpHpYp - 2*YfHpYp - YpHfYp + 2*YpHfYf +YfHpYf - YfHfYf
  # Multipliers
  YpmYf = yh1 - yh2
  YmYf = qxy$Ys - yh2
  d1 = colSums(YpmYf^2)
  d2 = d1^2
  num = 2*colSums(YmYf * YpmYf)
  divTemp = (Term1[hard] + Term2[hard])/d1[hard] - num[hard]*Term3[hard]/d2[hard]
  diverg = df
  diverg[hard] = df[hard] + divTemp
  return(diverg)
}
