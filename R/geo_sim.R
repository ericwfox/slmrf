#' Simulate geostatistical data
#'
#' Function that simulates zero-mean, spatially autocorrelated data with nugget
#' effect.  Simulation procedure uses the Cholesky decomposition method
#' (Cressie, 1993, p. 201).
#'
#' @param coord data.frame with x,y coordinates for simulation.
#' @param cov_type character.  Type of covariance model: \code{exponential},
#'   \code{gaussian}, or \code{spherical}; default is \code{exponential}.
#' @param nugg nugget parameter
#' @param parsil partial sill parameter
#' @param range range parameter
#'
#' @return numeric vector with simulated geostatistical data.
#'
#' @author Eric W. Fox
#'
#' @export
geo_sim <- function(coord, cov_type='exponential', nugg, parsil, range) {
  if(!(cov_type %in% c('exponential', 'gaussian', 'spherical'))) {
    stop('Must specify covariance function as one of the following types: exponential, gaussian, or spherical.')
  }

  # make covariance matrix
  # will not include nugget effect since this will be added later
  n <- nrow(coord)
  dist_mtx <- compute_distance(as.matrix(coord))
  covMat <- make_covmat(theta=c(0, parsil, range), dist_mtx=dist_mtx,
    is_log=F, include_nugg=F, type=cov_type)

  L <- chol(covMat)
  e <- rnorm(n)
  eps <- rnorm(n, sd=sqrt(nugg)) # vector of independent random errors
  z <- t(L) %*% e + eps

  return(as.numeric(z))
}
