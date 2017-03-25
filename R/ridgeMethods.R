#' @export
plot.cplr <- function(x, ...){
  matplot(x$lam, t(x$bhat), type = "l", ...)
}

#' @export
coef.cplr <- function(object, ...)
{
  drop(rbind(object$intercept, object$bhat))
}

#' @export
predict.cplr <- function(object, newdata, ...){
  drop(cbind(1,newdata) %*% coef.cplr(object))
}

#' @export
fitted.cplr <- function(object, ...){
  drop(object$fitted)
}

#' @export
residuals.cplr <- function(object, ...){
  drop(object$residuals)
}
