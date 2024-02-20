#' Minus 2 time the log-likelihood for reduced rank spatial linear model.
#'
#' @param theta Numeric vector of covariance parameters (use log transformed
#'   nugget, partial sill, and range, in that order).
#' @param y Numeric response vector.
#' @param X Design matrix.
#' @param dist_K distance matrix for knot locations.
#' @param dist_S distance matric between data locations and knots (\code{nrow}
#'   is number of data points and \code{ncol} is number of knots).
#' @param est_meth character: Estimation method: 'REML' for restricted maximum
#'   likelihood (default), or 'ML' for maximum likelihood.
#' @param cov_type character; type of covariance type of covariance model:
#'   \code{exponential}, \code{gaussian}, or \code{spherical}; default is
#'   \code{exponential}.
#'
#' @return Minus 2 times the log-likelihood for reduced rank spatial linear
#'   model using Sherman-Morrison-Woodbury formula.  See Cressie FRK paper
#'   equation 2.16 for decomposition.
#'
#' @author Eric W. Fox
#' @export
m2LL_smw <- function(theta, y, X, dist_K, dist_S, est_meth, cov_type='exponential') {
  if(!(est_meth %in% c('REML', 'ML'))) {
    stop('est_meth must either be REML or ML')
  }
  n <- length(y)
  p <- sum(svd(X)$d>1e-10)
  nugg <- exp(theta[1])

  K <- make_covmat(theta, dist_K, is_log = T, include_nugg = F, type=cov_type)
  S <- make_covmat(theta, dist_S, is_log = T, include_nugg = F, type=cov_type)
  Ki <- solve(K)
  A <- nugg * diag(1, nrow(S))
  Ai <- (1 / nugg) * diag(1, nrow(S))
  AiS <- (1 / nugg) * S # same as Ai%*%S
  covMat <-  S%*%Ki%*%t(S) + A
  Vi <- Ai - AiS%*%solve(K + t(S)%*%AiS)%*%t(AiS) # applying SMW

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
