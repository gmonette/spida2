#' setdiff as a binary operator
#'
#' The set difference of x less y
#'
#' @param x,y vectors of the same mode containing sequences of items
#' @return a vector of the same mode containing the difference of \code{x} less \code{y}
#' @export
"%less%" <- function(x,y) setdiff(x,y)
