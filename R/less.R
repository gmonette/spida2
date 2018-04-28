#' setdiff as a binary operator
#'
#' The set difference of x less y
#'
#' @param x,y vectors of the same mode containing sequences of items
#' @return a vector of the same mode containing the difference of \code{x} less \code{y}
#' @seealso \code{\link[spida2]{\%and\%}} and \code{\link[spida2]{\%or\%}} 
#' @export
"%less%" <- function(x,y) setdiff(x,y)
#' intersect as a binary operator
#'
#' The intersection of x and y
#'
#' @param x,y vectors of the same mode containing sequences of items
#' @return a vector of the same mode containing the intersection of \code{x} less \code{y}
#' @seealso \code{\link[spida2]{\%less\%}} and \code{\link[spida2]{\%or\%}} 
#' @export
"%and%" <- function(x,y) intersect(x,y)
#' union as a binary operator
#'
#' The union of x and y
#'
#' @param x,y vectors of the same mode containing sequences of items
#' @return a vector of the same mode containing the union of \code{x} less \code{y}
#' @seealso \code{\link[spida2]{\%and\%}} and \code{\link[spida2]{\%less\%}} 
#' @export
"%or%" <- function(x,y) union(x,y)
