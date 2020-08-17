#' Version of lapply that turns the input list into a list with elements of length one
#'
#' Allows lapply to work on elements that are themselves
#' lists of length one so that the name of each element
#' is available within the function called by \code{\link{lapply}}.
#' 
#' @param x a vector (atomic or list each of whose elements
#'      will be made into a list of length one and inserted
#'      into a list to be passed to \code{\link{lapply}}.
#' @param FUN a function to be applied to each element
#'      of the list formed from \code{x}. If \code{FUN} is missing,
#'      it defaults to the identity function, thus returning
#'      x turned into a list of elements of length one.
#' @param \dots additional variables passed to \code{\link{lapply}}
#' if \code{FUN} returns a vector of length 1.
#'
#' @export
llapply <- function(x, FUN = function(x) x, ...) {
  listify <- function(x) {
    ret <- vector(length(x), mode = 'list')
    for(i in seq_along(x)) ret[[i]] <- x[i]
    names(ret) <- names(x)
    ret
  }
  x <- listify(x)
  lapply(x, FUN, ...)
} 
