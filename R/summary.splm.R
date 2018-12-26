#' Summarize Spatial Linear
#' 
#' \code{summary} method for class 'splm'
#' 
#' @param object an object of class \code{splm}
#'   
#' @return data frame with spatial regression summary of class 'summary.splm'
#'   
#' @export
summary.splm <- function(object, ...) {
  if(class(object) != 'splm') {
    stop('Not a splm object.')
  }
  
  n <- length(object$y)
  p <- ncol(object$X) # includes intercept as a parameter
  
  tval <- object$betaHat / object$seBetaHat
  
  coefs <- data.frame(est = object$betaHat,
    se = object$seBetaHat,
    t = tval,
    pval = 2 * pt(-abs(tval), df = n-p, lower.tail = TRUE)) 
  
  out <- list(coefs=coefs)
  
  class(out) <- 'summary.splm'
  out
}
