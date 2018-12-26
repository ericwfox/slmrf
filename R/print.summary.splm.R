#' Print method of spatial linear model
#' 
#' @param x an object of class \code{summary.splm}
#'   
#' @author Eric W. Fox
#' @export
print.summary.splm <- function(x, ...) {
  printCoefmat(x$coefs, P.values=T, has.Pvalue = T)
}
