#' Compute root-mean-square error
#'
#' @param y Numeric data vector.
#' @param y_pred Numeric vector of predictions.
#'
#' @return RMSE
#'
#' @author Eric Fox
#' @export
compute_rmse <- function(y, y_pred) {
  n <- length(y)
  sqrt((1 / n) * sum((y - y_pred)^2))
}
