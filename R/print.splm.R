#' Print method of spatial linear model
#' 
#' @param x an object of class \code{splm}
#' @param digits how many significant digits to be used for numeric values. 
#'   Note that at least 2 places after the decimal will always be printed.
#'   
#' @author Eric W. Fox
#' @export
print.splm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat('Estimated parameters for', x$cov_type, 'covariance model:\n',
    'nugget:', format(x$nugg, digits = digits, nsmall=2), '\n',
    'partial sill:', format(x$parsil, digits = digits, nsmall=2), '\n',
    'range:', format(x$range, digits = digits, nsmall=2), '\n\n')
  
  cat('Estimation method:', x$est_meth, '\n\n')
  if(!is.null(x$num_knots)) {
    cat('Using reduced rank method with ', x$num_knots, ' knots.\n')
  }
}
