#' Creates an autonym
#'
#' Uses x to name x. Useful in when applying \code{\link{lapply}} to a character vector. 
#' Should probably be called 'autonym'. Superseded by \code{\link{name}}.
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
#' It can be thought of as the verb corresponding to the noun 'names'.
#'
#' @param x object to named with a 'names' attribute
#' @param nam vector of names or a function
#' @param ... additional arguments to used by \code{nam} if it is a function
#' @return x with a modified name. If \code{nam} is a function then the names attribute
#'         is modified with \code{names(x) <- nam(names(x), ...)}
#' @examples
#' \dontrun{
#' x <- as.list(letters[1:3])
#' x
#' name(x)
#' library(magrittr)
#' x %>% name %>% name(toupper)
#' x %>% name %>% name(sub_, '(.*)', 'time_\\1')
#' }
#' @export
name <- function(x, nam = x,...) {
  if(is.function(nam)) nam <- nam(names(x),...)
  names(x) <- nam
  x
}
#' @rdname name
#' @export
col_name <- function(x, nam = x,...) {
  if(is.function(nam)) nam <- nam(colnames(x),...)
  colnames(x) <- nam
  x
}
#' @rdname name
#' @export
row_name <- function(x, nam = x,...) {
  if(is.function(nam)) nam <- nam(rownames(x),...)
  rownames(x) <- nam
  x
}
#' 
#' Methods providing versions of sub and gsub for pipelines and factors
#' 
#' The first argument of sub_  and gsub_ is the object to 
#' be modified instead of the pattern to be matched. Thus \code{sub_}
#' and \code{gsub_} can be used as generic functions that dispatch on
#' the class of the first argument. This allows them to recognize factors
#' and work with the levels attribute instead of coercing the
#' factor to a character vector. Two consequences are that a
#' factor, instead of a character vector, is returned, and that the
#' original order of the levels is preserved unless the number of 
#' levels is reduced. 
#' 
#' @param x object to change. If x is a factor, the substitution is performed on its levels attribute.
#' @param pattern the 'from' regular expression
#' @param replacement the 'to' regular expression
#' @param ... other arguments for  \code{\link{sub}} or  \code{\link{gsub}}
#' @return Value of \code{sub(pattern, replacement, x, ...)} or \code{gsub(pattern, replacement, x, ...)} 
#' @seealso \code{\link{sub}} \code{\link{sub_}}
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
#' @describeIn sub_ default method
#' @export
sub_.default <- function(x, pattern, replacement,...) {
  replacement <- rep_len(replacement, length.out = length(pattern) )
  for(i in seq_along(pattern)){
    x <- sub(pattern[i],replacement[i],x, ...)
  }
  x
}
#' @describeIn sub_ method for factor objects
#' @export
sub_.factor <- function(x, pattern, replacement, ...) {
  levels(x) <- sub_(levels(x), pattern, replacement, ...)
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
#' original order of the levels is preserved unless the number of 
#' levels is reduced. 
#' 
#' @param x object to change. If x is a factor, the substitution is performed on its levels attribute.
#' @param pattern the 'from' regular expression
#' @param replacement the 'to' regular expression
#' @param ... other arguments for  \code{\link{sub}} or  \code{\link{gsub}}
#' @return Value of \code{sub(pattern, replacement, x, ...)} or \code{gsub(pattern, replacement, x, ...)} 
#' @seealso \code{\link{sub}} \code{\link{gsub_}}
#' @examples
#' \dontrun{
#' library(magrittr)
#' x <- as.list(1:4)
#' x  %>% name(letters[1:4])
#' x  %>% name(letters[1:4]) %>% name(toupper)
#' x  %>% name(letters[1:4]) %>% 
#'        name(gsub_, '(.*)', 'Time_\\1')  %>% unlist
#' }
#' @export
gsub_ <- function(x, pattern, replacement, ...) UseMethod("gsub_")
#' @describeIn gsub_ default method
#' @export
gsub_.default <- function(x, pattern, replacement,...) {
  replacement <- rep_len(replacement, length.out = length(pattern) )
  for(i in seq_along(pattern)){
    x <- gsub(pattern[i],replacement[i], x, ...)
  }
  x
}
#' @describeIn gsub_ method for factors
#' @export
gsub_.factor <- function(x, pattern, replacement, ...) {
  levels(x) <- gsub_(levels(x),pattern,replacement, ...)
  x
}
#' Extract or remove a substring
#' 
#' Extract of remove a substring matching a regular expression.
#' 
#' \code{getex} and \code{getex_} can work as a pair to gradually extract information
#' from a string.
#' 
#' @param x object from which strings are extracted. Will be coerced to a
#'        character by \code{\link{regexpr}}.
#' @param pattern a regular expression whose matches in \code{x} are returned.
#' @param ... other arguments passed to \code{\link{regexpr}}.
#' @return \code{getex} returns the substrings matched by \code{pattern} in \code{x}. 
#'         \code{getex_} returns the portion of \code{x} that is not matched.
#' @examples
#' library(magrittr)
#' # The following could be done with 'strsplit' but 'getex' can handle messier cases
#' progs <- c('SC BSc Honours Biology', 'Arts BA Ordinary Economics')
#' faculty <- progs  %>%  getex('^[^ ]+') %>% toupper
#' progs <- progs  %>% getex_('^[^ ]+ ')
#' degree <- progs %>% getex('^[^ ]+')
getex <- function(x, pattern,...) {
  m <- regexpr(pattern, x, ...)
  substring(x,m, attr(m,'match.length'))
}
getex_ <- function(x, pattern,...) {
  sub(pattern, '', x)
}