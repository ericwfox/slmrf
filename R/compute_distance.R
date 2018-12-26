#' Compute Euclidian distance matrix.
#'
#' @param x1 Matrix containing coordinates for the first set of locations.
#' @param x2 Matrix containing coordinates for the second set of locations.
#'
#' @return A matrix with Euclidean distances between two sets of coordinates
#'   \code{x1} and \code{x2}.  If the second set of coordinates \code{x2} is not
#'   specified then a matrix pairwise distances for the first set of coordinates
#'   \code{x1} are returned.
#'
#' @examples
#' x1 <- matrix(1:6, ncol =2, byrow = TRUE)
#' x2 <- matrix(7:12, ncol = 2, byrow = TRUE)
#' compute_distance(x1)
#' compute_distance(x1,x2)
#'
#' @author Eric W. Fox
#' @export
compute_distance <- function(x1, x2 = NULL) {
  if(!is.matrix(x1)) {
    warning('x1 is not a matrix, attempting to coerce.')
    x1 <- as.matrix(x1)
  }
  if(!is.null(x2) & !is.matrix(x2)) {
    warning('x2 is not a matrix, attempting to coerce.')
    x2 <- as.matrix(x2)
  }
  
  if(is.null(x2)) {
    a <- outer(x1[,1], x1[,1], '-')
    b <- outer(x1[,2], x1[,2], '-')
  } else {
    a <- outer(x1[,1], x2[,1], '-')
    b <- outer(x1[,2], x2[,2], '-')
  }
  d <- sqrt(a^2 + b^2)
  return(d)
}
