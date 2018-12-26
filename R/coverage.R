#' Check coverage of prediction intervals
#'
#' @param y Numeric data vector.
#' @param pred Numeric vector of predictions.
#' @param se Numeric vector of prediction SEs.
#' @param cl Confidence level (default is 0.9).
#'
#' @return Coverage for prediction intervals.
#' @author Eric W. Fox
#' @export
coverage <- function(y, pred, se, cl = 0.9) {
  count <-  0
  alpha <- 1-cl
  critval <- qnorm(1 - alpha / 2)
  n <- length(y)
  for(i in 1:n) {
    ci_lb = pred[i] - critval * se[i]
    ci_ub = pred[i] + critval * se[i]
    if((y[i] >= ci_lb) & (y[i] <= ci_ub)) {
      count <- count + 1
    }
  }
  count / n
}
