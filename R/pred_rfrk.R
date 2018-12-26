#' Kriging predictions and SEs for RFRK model
#'
#' Evaluate kriging predictions and standard errors for random forest regression
#' kriging model
#'
#' @param rf object of class 'randomForest'
#' @param rfrk object of class 'rfrk'
#' @param obsv_coord Data frame  or matrix with x,y coordinates for
#'   observations.
#' @param pred_coord Data frame or matrix with x,y coordinates for predictions.
#' @param Xp matrix of covariates at prediction sites.
#' @param row_names character or numeric vector specify row names for the
#'   predictions
#' @param computeSE logical; whether to compute standard errors, default is
#'   \code{TRUE}
#' @param scale scaling factor for distance matrix.  The distance matrix will be
#'   divided by this number. Dafault is 1 for no scaling.  Useful when
#'   converting Alber's coordinates from meters to kilometers; in this case, set
#'   scale to 1000.
#'
#' @return Matrix with RFRK predictions and SEs.
#' @export
pred_rfrk <- function(rf, rfrk, obsv_coord, pred_coord, Xp,
  row_names = NULL, computeSE = T, scale = 1) {

  pred_out <- matrix(NA, nrow = nrow(Xp), ncol = 2) #initialize
  # make kriging predictions for residuals
  r <- rfrk$resid
  Vi <- rfrk$Vi
  sill <- rfrk$sill
  obsv_coord <- as.matrix(obsv_coord)
  pred_coord <- as.matrix(pred_coord)
  dist_mtx_pred <- compute_distance(obsv_coord, pred_coord) / scale
  Vpred <- make_covmat(rfrk$params, dist_mtx_pred, is_log=F, include_nugg=F)
  resid_pred <- t(Vpred) %*% (Vi %*% (r))
  # make random forest predictions
  rf_pred <- predict(rf, newdata = Xp)
  pred_out[,1] <- rf_pred + resid_pred
  if(computeSE) {
    pred_out[, 2] <- sqrt(sill - apply((Vi %*% Vpred) * Vpred, 2, sum))
  }

  colnames(pred_out) <- c('pred', 'predSE')
  if(!is.null(row_names)) {
    rownames(pred_out) <- row_names
  }
  return(pred_out)
}