#' Check coverage of prediction intervals
#' 
#' In this version user must supply prediction intervals (e.g., from 
#' \code{lm.predict})
#' 
#' @param y numeric data vector
#' @param ci_lb numveric vector of lower bounds for prediction intervals
#' @param ci_ub numberic vector of upper bounds for prediction intervals
#'   
#' @return Coverage for prediction intervals.
#' 
#' @author Eric W. Fox
#' 
#' @export  
coverage2 <- function(y, ci_lb, ci_ub) {
  count <-  0
  n <- length(y)
  for(i in 1:n) {
    if((y[i] >= ci_lb[i]) & (y[i] <= ci_ub[i])) {
      count <- count + 1
    }
  }
  count / n
}
