#' Fit RFRK
#'
#' Estimate a random forest residual kriging model
#'
#' @param rf random forest model
#' @param y Numeric response vector.
#' @param dist_mtx Distance matrix used to compute covariance matrix.
#' @param theta_ini numeric vector. Initial values for covariance parameters:
#'   logarithm transformed nugget, partial sill, and range, in that order.  If
#'   not specified uses default initialization.
#'
#' @return object of type `rfrk'.
#'
#' @note assumes exponential covariance function
#'
#' @export
fit_rfrk <- function(rf, y, dist_mtx, theta_ini=NULL) {
  resid <- y - predict(rf)

  # default parameter initialization
  if(is.null(theta_ini)) {
    theta_ini <- log(c('nugg' = 0.5*var(resid), 'parsil' = 0.5*var(resid), 'range' = mean(dist_mtx)))
  }

  # run optimization
  start_time <- Sys.time()
  parmest <- optim(theta_ini, m2LL_sk, r=resid, dist_mtx=dist_mtx)
  end_time <- Sys.time()
  est_time <- end_time - start_time

  params <- exp(parmest$par)
  covMat <- make_covmat(params, dist_mtx, is_log=FALSE,
    include_nugg = TRUE) # exponential is default
  Vi <- solve(covMat)

  out <- list(
    resid = resid,
    params = params,
    nugg = params[1],
    parsil = params[2],
    range = params[3],
    sill = params[1] + params[2],
    covMat = covMat,
    Vi = Vi,
    optim_out = parmest,
    est_time = est_time
  )

  class(out) <- 'rfrk'
  return(out)
}