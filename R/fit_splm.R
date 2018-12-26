#' Estimate a spatial linear model
#'
#' Estimate parameters for spatial linear model and form object of class
#' \code{splm}.
#'
#' The purpose of this function is to automate the process of fitting a spatial
#' model using \code{optim} and forming an object of class \code{splm}, which
#' contains the attributes of that fitted model, using \code{make_splm}.
#'
#' @param X Design matrix.  Must include column of 1's for intercept.
#' @param y Numeric response vector.
#' @param cov_type character.  Type of covariance model: \code{exponential},
#'   \code{gaussian}, or \code{spherical}.
#' @param theta_ini numeric vector. Initial values for covariance parameters:
#'   logarithm transformed nugget, partial sill, and range, in that order.  If
#'   not specified uses default initialization.
#' @param dist_mtx Distance matrix used to compute covariance matrix.
#' @param dist_K Distance matrix for knot locations.
#' @param dist_S Distance matric between data locations and knots (\code{nrow}
#'   is number of data points and \code{ncol} is number of knots).
#' @param est_meth character. Estimation method: 'REML' for restricted maximum
#'   likelihood, or 'ML' for maximum likelihood.
#' @param include_hessian logical.  Whether or not to compute hessian matrix
#'   when running \code{optim}.
#' @param knot_coord data frame or matrix. Coordinates of knots used for reduced
#'   rank method. Optional to include (coordinates are just saved in \code{splm}
#'   object for later reference).
#'
#' @return List of class \code{splm} containg features of fitted spatial linear
#'   model (e.g., covariance matrix, predictor coefficients, ...).
#'
#' @author Eric W. Fox
#'
#' @export
fit_splm <- function(X, y, cov_type, theta_ini = NULL, dist_mtx = NULL, dist_K = NULL,
  dist_S = NULL, est_meth, include_hessian = F, knot_coord = NULL) {

  if(!is.matrix(X)) {
    warning('X is not a matrix, attempting to coerce.')
    X <- as.matrix(X)
  }

  # default parameter initialization
  if(is.null(theta_ini)) {
    if(!is.null(dist_mtx)) {
      theta_ini <- log(c('nugg' = 0.5*var(y), 'parsil' = 0.5*var(y), 'range' = mean(dist_mtx)))
    } else {
      stop('No default value for theta_ini, user must supply.')
    }
  }

  # estimate model using full rank covariance matrix
  if(!is.null(dist_mtx)) {
    # run optimization
    start_time <- Sys.time()
    parmest <- optim(theta_ini, m2LL, y = y, X = X, dist_mtx = dist_mtx,
      est_meth = est_meth, cov_type = cov_type, hessian = include_hessian)
    end_time <- Sys.time()

    # form splm object
    splm <- make_splm(exp(parmest$par), X = X, y = y, cov_type = cov_type, dist_mtx = dist_mtx,
      m2loglike = parmest$value, est_meth = est_meth, hessian = parmest$hessian, optim_out = parmest)
    splm$est_time <- end_time - start_time
  } else if(!is.null(dist_K) & !is.null(dist_S)) {
    # run optimization
    start_time <- Sys.time()
    parmest <- optim(theta_ini, m2LL_smw, y = y, X = X, dist_K = dist_K, dist_S = dist_S,
      est_meth = est_meth, cov_type = cov_type, hessian = include_hessian)
    end_time <- Sys.time()

    # form splm object
    splm <- make_splm(exp(parmest$par), X = X, y = y, cov_type = cov_type,
      m2loglike = parmest$value, est_meth = est_meth, hessian = parmest$hessian,
      dist_K = dist_K, dist_S = dist_S, knot_coord = knot_coord, optim_out = parmest)
    splm$est_time <- end_time - start_time
  } else {
    stop('Need to provide distance matrix.')
  }

  return(splm)
}
