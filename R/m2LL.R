#' Minus 2 time the log-likelihood for spatial linear model.
#'
#' @param theta Numeric vector of covariance parameters (use log transformed
#'   nugget, partial sill, and range, in that order).
#' @param y Numeric response vector.
#' @param X Design matrix.
#' @param dist_mtx Matrix of pairwise distances between data locations.
#' @param est_meth character: Estimation method: 'REML' for restricted maximum
#'   likelihood (default), or 'ML' for maximum likelihood.
#' @param cov_type character; type of covariance type of covariance model:
#'   \code{exponential}, \code{gaussian}, or \code{spherical}; default is
#'   \code{exponential}.
#' @param boundRange whether to bound effective range at max distance.
#'
#' @return Minus 2 times the log-likelihood for spatial linear model.
#'
#' @author Eric W. Fox
#' @export
m2LL <- function(theta, y, X, dist_mtx, est_meth, cov_type = 'exponential', boundRange=F) {
  if(!(est_meth %in% c('REML', 'ML'))) {
    stop('est_meth must either be REML or ML')
  }

  if(boundRange==T) {
    range <- exp(theta[3])
    if(3*range > max(dist_mtx)) {return(10^9)}
  }

  n <- length(y)
  p <- sum(svd(X)$d>1e-10)

  covMat <- make_covmat(theta, dist_mtx, is_log=TRUE,
    include_nugg = TRUE, type = cov_type)
  Vi <- solve(covMat)

  logDetV <- as.numeric(determinant(covMat, logarithm = TRUE)$modulus)
  XViX <- t(X)%*%Vi%*%X
  logDetXViX <- as.numeric(determinant(XViX, logarithm = TRUE)$modulus)
  covBetaHat <- solve(XViX)
  betaHat <- covBetaHat%*%t(X)%*%Vi%*%y
  r <- y - X%*%betaHat
  rVir <- t(r)%*%Vi%*%r
  minus2LL <- logDetV + rVir + n*log(2*pi)
  if(est_meth == "REML") {
    minus2LL <- minus2LL + logDetXViX - p*log(2*pi)
  }

  return(as.numeric(minus2LL))
}
