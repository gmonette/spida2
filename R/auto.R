#' Creates an autonym
#'
#' Uses x to name x. Useful in when applying \code{\link{lapply}} to a character vector. Should probably be called 'autonym'.
#'
#' @param x object that can be used to name itself.
#' @param force (default FALSE) should x replace its name if it already has one.
#' @return x with names equal to x
#' @export
auto <- function(x, force = FALSE) {
  if(is.null(names(x)) || force) names(x) <- x
  x
}
#' @rdname auto
#' @export
autonym <- auto
