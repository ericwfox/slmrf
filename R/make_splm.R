#' Make a spatial linear model
#'
#' Make spatial linear model object using estimated covariance parameters from
#' \code{optim()}.
#'
#' @param params Numeric vector of fitted covariance parameter values (nugget,
#'   partial sill, and range, in that order).
#' @param X Design matrix.  Must include column of 1's for intercept.
#' @param y Numeric response vector.
#' @param cov_type character.  Type of covariance model: \code{exponential},
#'   \code{gaussian}, or \code{spherical}.
#' @param dist_mtx Distance matrix used to compute covariance matrix.
#' @param m2loglike Minus 2 times the log-likelihood (from \code{optim()}).
#' @param est_meth character. Estimation method used for \code{m2LL}: 'REML' for
#'   restricted maximum likelihood, or 'ML' for maximum likelihood.
#' @param hessian Hessian matrix from \code{optim()}.
#' @param dist_K Distance matrix for knot locations.
#' @param dist_S Distance matric between data locations and knots (\code{nrow}
#'   is number of data points and \code{ncol} is number of knots).
#' @param knot_coord data frame or matrix. Coordinates of knots used for reduced
#'   rank method.  Optional to include (coordinates are just saved in 'splm'
#'   object for later reference).
#' @param optim_out list output from \code{optim()}.
#'
#' @return List of class \code{splm} containg features of fitted spatial linear
#'   model (e.g., covariance matrix, predictor coefficients, ...).
#'
#' @seealso \code{\link{m2LL}}, \code{\link[stats]{optim}}
#'
#' @author Eric W. Fox
#'
#' @export
make_splm <- function(params, X, y, cov_type, dist_mtx = NULL, m2loglike = NULL, est_meth=NULL,
  hessian = NULL, dist_K = NULL, dist_S = NULL, knot_coord = NULL, optim_out = NULL) {
  if(!(cov_type %in% c('exponential', 'gaussian', 'spherical'))) {
    stop('Must specify covariance function as one of the following types: exponential, gaussian, or spherical.')
  }
  if(!is.matrix(X)) {
    warning('X is not a matrix, attempting to coerce.')
    X <- as.matrix(X)
  }

  if(!is.null(dist_mtx)) {
    covMat <-  make_covmat(params, dist_mtx, is_log = FALSE, include_nugg = TRUE, type = cov_type)
    Vi <- solve(covMat)
  } else if(!is.null(dist_K) & !is.null(dist_S)) {
    K <- make_covmat(params, dist_K, is_log = FALSE, include_nugg = FALSE, type = cov_type)
    S <- make_covmat(params, dist_S, is_log = FALSE, include_nugg = FALSE, type = cov_type)
    Ki <- solve(K)
    nugg <- params[1]
    A <- nugg * diag(1, nrow(S))
    Ai <- (1 / nugg) * diag(1, nrow(S))
    AiS <- (1 / nugg) * S # same as Ai%*%S
    covMat <-  S%*%Ki%*%t(S) + A
    Vi <- Ai - AiS%*%solve(K + t(S)%*%AiS)%*%t(AiS) # applying SMW
  } else {
    stop('Need to provide distance matrix.')
  }

  XViX <- t(X)%*%Vi%*%X
  covBetaHat <- solve(XViX)
  betaHat <- covBetaHat%*%t(X)%*%Vi%*%y
  seBetaHat <- as.numeric(sqrt(diag(covBetaHat)))
  sill <- params[1] + params[2]

  out <- list(
    y = y,
    X = X,
    cov_type = cov_type,
    params = params,
    nugg = params[1],
    parsil = params[2],
    range = params[3],
    sill = sill,
    covMat = covMat,
    Vi = Vi,
    XViX = XViX,
    betaHat = betaHat,
    covBetaHat = covBetaHat,
    seBetaHat = seBetaHat,
    m2loglike = m2loglike,
    est_meth = est_meth,
    hessian = hessian,
    optim_out = optim_out
  )
  if(!is.null(dist_K)) {
    out$num_knots <- nrow(dist_K)
    out$knot_coord <- knot_coord
  }

  class(out) <- 'splm'
  invisible(out)
}
