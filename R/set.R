###
### set: part
###

#' Intersection and symmetric differences of two sets
#'
#' @param A,B vectors to be coerced to character
#' @export
part <- function( A, B) {
  # partitions of two sets
  if ( is.factor(A)) A <- as.character(A)
  if (is.factor(B)) B <- as.character(B)
  ret <- list( "A and B" = intersect(A,B),
               "A - B" = setdiff( A, B),
               "B - A" = setdiff(B, A))
  attr(ret,"N") <- sapply(ret, length)
  ret
}
