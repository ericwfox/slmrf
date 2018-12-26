#' Kriging predictions and errors
#'
#' Function that makes krigging predictions using a fitted spatial linear model
#' and an external prediction data set.
#'
#' @param object object of class 'splm'
#' @param obsv_coord Data frame  or matrix with Alber's x,y coordinates for
#'   observations.
#' @param pred_coord Data frame or matrix with Alber's x,y coordinates for
#'   prediction.
#' @param Xp matrix of covariates at prediction sites for universal kriging.
#'   Make sure that the first column of this matrix contains 1's for the
#'   intercept. The default is \code{NULL}, which assumes a spatial model with
#'   no predictors (the function will generate a matrix with a column of 1's for
#'   the intercept).
#' @param row_names character or numeric vector specify row names for the
#'   predictions
#' @param computeSE logical; whether to compute standard errors, default is
#'   \code{TRUE}
#' @param scale scaling factor for distance matrix.  The distance matrix will be
#'   divided by this number. Dafault is 1 for no scaling.  Useful when
#'   converting Alber's coordinates from meters to kilometers; in this case, set
#'   scale to 1000.
#'
#'
#' @return Matrix with kriging predictions and standard errors.
#'
#' @author Jay Ver Hoef (edited by Eric Fox)
#' @export
pred_krige_batch = function(object, obsv_coord, pred_coord, Xp = NULL,
  row_names = NULL, computeSE = TRUE, scale = 1) {
  if(class(object) != 'splm') {
    stop('Not a splm object.')
  }

  y <- object$y
  X <- object$X
  Vi <- object$Vi
  covb <- object$covBetaHat
  betaHat <- object$betaHat
  sill <- object$sill
  obsv_coord <- as.matrix(obsv_coord)
  pred_coord <- as.matrix(pred_coord)
  dist_mtx_pred <- compute_distance(obsv_coord, pred_coord) / scale
  if(is.null(Xp)) {
    Xp <- as.matrix(rep(1, nrow(pred_coord)))
  }

  Vpred <- make_covmat(object$params, dist_mtx_pred, is_log = FALSE,
    include_nugg = FALSE, type = object$cov_type)
  pred_out <- matrix(NA, nrow = nrow(Xp), ncol = 2)
  pred_out[,1] <- apply(as.vector((Vi %*% y)) * Vpred, 2, sum) +
    Xp %*% betaHat - t(Vpred) %*% (Vi %*% (X %*% betaHat))
  if(computeSE) {
    pred_out[,2] <- sqrt(rep(sill, times = nrow(Xp)) -
        apply((Vi %*% Vpred) * Vpred, 2, sum) +
        apply((covb %*% t(Xp)) * t(Xp), 2, sum) -
        2*apply((covb %*% t(Xp)) * (t(X) %*% Vi %*% Vpred), 2, sum) +
        apply((covb %*% t(X) %*% Vi %*% Vpred) * (t(X) %*% Vi %*% Vpred), 2, sum))
  }
  colnames(pred_out) <- c('pred', 'predSE')
  if(!is.null(row_names)) {
    rownames(pred_out) <- row_names
  }
  return(pred_out)
}
