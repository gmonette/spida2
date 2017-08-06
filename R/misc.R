#' Change NAs to FALSE
#' 
#' @param x vector, possibly with NAs
#' @export
na2f <- function(x) {
  x[is.na(x)] <- FALSE
  x
}
#' Change NAs to TRUE
#' 
#' @param x vector, possibly with NAs
#' @export
na2t <- function(x) {
  x[is.na(x)] <- TRUE
  x
}
#' Pipe from magrittr
#' 
#' @export
"%>%" <- magrittr::`%>%`

#' Transform NAs to 0
#'
#' @param x vector, possibly with NAs
#' @export
na20 <- function(x) {
  x[is.na(x)] <- 0
  x
}
