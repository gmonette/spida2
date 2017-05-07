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
#' Modify the name of an object and return the renamed object
#' 
#' This function is well suited to modifying the name of an object in a pipeline.
#'
#' @param x object
#' @param nam name
#' @return x with a modified name
#' @examples
#' \dontrun{
#' x <- as.list(letters[1:3])
#' x
#' library(magrittr)
#' x  %>% name
#' gsub_ <- function(x,pat,to,...) gsub(pat,to,x,...)
#' x %>% name()
#' x  %>% name %>% name(toupper)
#' x  %>% name %>% name(sub_, '(.*)', 'time_\\1')
#' }
#' @export
name <- function(x, nam = x,...) {
  if(is.function(nam)) nam <- nam(names(x),...)
  names(x) <- nam
  x
}
#' Methods providing versions of sub and gsub for pipelines and factors
#' 
#' The first argument of sub_  and gsub_ is the object to 
#' be modified instead of the pattern to be matched. Thus \code{sub_}
#' and \code{gsub_} can be used as generic functions that dispatch on
#' the class of the first argument. This allows them to recognize factors
#' and work with the levels attribute instead of coercing the
#' factor to a character vector. Two consequences are that a
#' factor, instead of a character vector, is returned, and that the
#' original order of the levels is preserved. 
#' 
#' @param x object to change
#' @param pattern the 'from' regular expression
#' @param replacement the 'to' regular expression
#' @param ... other arguments for  \code{\link{sub}} or  \code{\link{gsub}}
#' @return Value of \code{sub(pattern, replacement, x, ...)} or \code{gsub(pattern, replacement, x, ...)} 
#' @seealso \code{\link{sub}}
#' @rdname sub_
#' @examples
#' \dontrun{
#' library(magrittr)
#' x <- as.list(1:4)
#' x  %>% name(letters[1:4])
#' x  %>% name(letters[1:4]) %>% name(toupper)
#' x  %>% name(letters[1:4]) %>% 
#'        name(sub_, '(.*)', 'Time_\\1')  %>% unlist
#' }
#' @export
sub_ <- function(x, pattern, replacement, ...) UseMethod("sub_")
#' @rdname sub_
sub_.default <- function(x, pattern, replacement,...) {
  replacement <- rep_len(replacement, length.out = length(pattern) )
  for(i in seq_along(pattern)){
    x <- sub(pattern[i],replacement[i],x, ...)
  }
  x
}
#' @rdname sub_
sub_.factor <- function(x, pattern, replacement, ...) {
  levels(x) <- sub_(levels(x), pattern, replacement, ...)
  x
}
#' @rdname sub_
#' @export
gsub_ <- function(x, pattern, replacement, ...) UseMethod("gsub_")
#' @rdname sub_
gsub_.default <- function(x, pattern, replacement,...) {
  replacement <- rep_len(replacement, length.out = length(pattern) )
  for(i in seq_along(pattern)){
    x <- gsub(pattern[i],replacement[i], x, ...)
  }
  x
}
#' @rdname sub_
gsub_.factor <- function(x, pattern, replacement, ...) {
  levels(x) <- gsub_(levels(x),pattern,replacement, ...)
  x
}
