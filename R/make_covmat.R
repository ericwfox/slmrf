#' Make covariance matrix.
#'
#' Make covariance matrix with different types of covariance models.
#'
#' @param theta numeric vector with 3 covariance parameters: nugget, partial
#'   sill, and range (in that order).
#' @param dist_mtx  Matrix of pairwise distances between data locations.
#' @param is_log logical. If \code{TRUE} then the values of \code{theta} are the
#'   natural logarithms of the parameters.
#' @param include_nugg logical. Whether or not to include the nugget along the
#'   diagonal.
#' @param type character. Type of covariance model: \code{exponential},
#'   \code{gaussian}, or \code{spherical}; default is \code{exponential}.
#'
#' @return Covariance matrix.
#'
#' @author Eric W. Fox
#' @export
make_covmat <- function(theta, dist_mtx, is_log, include_nugg, type = 'exponential') {
  if(!(type %in% c('exponential', 'gaussian', 'spherical'))) {
    stop('Covariance function must be one of the following types: exponential, gaussian, or spherical.')
  }
  if(include_nugg) {
    if(ncol(dist_mtx) != nrow(dist_mtx)) {
      warning('Are you sure you want to include the nugget in the covariance matrix?')
    }
  }

  if(is_log) {
    nugg <- exp(theta[1])
    parsil <- exp(theta[2])
    range <- exp(theta[3])
  } else {
    nugg <- theta[1]
    parsil <- theta[2]
    range <- theta[3]
  }

  if(type == 'exponential') {
    mtx <- parsil * exp(-1 * dist_mtx / range)
    if(include_nugg) {
      diag(mtx) <- parsil + nugg
    }
  }

  if(type == 'gaussian') {
    mtx <- parsil * exp(-1 * (dist_mtx / range)^2)
    if(include_nugg) {
      diag(mtx) <- parsil + nugg
    }
  }

  if(type == 'spherical') {
    mtx <- parsil * (1 - 1.5*(dist_mtx / range) + 0.5*(dist_mtx / range)^3)
    mtx[dist_mtx > range] = 0
    if(include_nugg) {
      diag(mtx) <- parsil + nugg
    }
  }

  return(mtx)
}
