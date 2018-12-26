#' Minus 2 time the log-likelihood for zero-mean simple kriging
#'
#' log-likelihood for residual rfrk model
#'
#' @param theta Numeric vector of covariance parameters (use log transformed
#'   nugget, partial sill, and range, in that order).
#' @param r Numeric vector of residuals
#' @param boundRange whether to bound effective range at max distance.
#'
#' @return Minus 2 times the log-likelihood for RFRK model.
m2LL_sk <- function(theta, r, dist_mtx, boundRange=T) {
  if(boundRange==T) {
    range <- exp(theta[3])
    if(3*range > max(dist_mtx)) {return(10^9)}
  }

  n <- length(r)
  covMat <- make_covmat(theta, dist_mtx, is_log=TRUE,
    include_nugg = TRUE) # exponential is default
  Vi <- solve(covMat)
  logDetV <- as.numeric(determinant(covMat, logarithm = TRUE)$modulus)

  minus2LL <- n*log(2*pi) + logDetV + t(r)%*%Vi%*%r
  return(as.numeric(minus2LL))
}