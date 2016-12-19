#' Change NAs to FALSE
#' 
#' @param x
#' @export
na2f <- function(x) {
  x[is.na(x)] <- FALSE
  x
}
#' Change NAs to TRUE
#' 
#' @param x
#' @export
na2t <- function(x) {
  x[is.na(x)] <- TRUE
  x
}