#' Creates an autonym
#'
#' Uses x to name x. Useful in when applying \code{\link{lapply}} to a character vector. Should probably be called 'autonym'.
#'
#' @param x object that can be used to name itself.
#' @param force (default FALSE) should x replace its name if it already has one.
#' @return x with names equal to x
#' @seealso \code{\link{anon}}
#' @export
auto <- function(x, force = FALSE) {
  if(is.null(names(x)) || force) names(x) <- x
  x
}
#' @rdname auto
#' @export
autonym <- auto
#'
#' Add or modify a name of an object
#'
#' @param x object
#' @param nam name
#' @return x with a modified name
#' @examples
#' x <- as.list(letters[1:3])
#' x
#' library(magrittr)
#' x  %>% name
#' gsub_ <- function(x,pat,to,...) gsub(pat,to,x,...)
#' x %>% name()
#' x  %>% name %>% name(toupper)
#' x  %>% name %>% name(gsub_, '(.*)', 'time_\\1')
#' @export
name <- function(x, nam = x,...) {
  if(is.function(nam)) nam <- nam(names(x),...)
  names(x) <- nam
  x
}
