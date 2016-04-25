#' Displaying p-values
#'
#' @param x a n x 2 matrix or data frame
#' @param digits (default: 5) number of digits for rounding coefficient
#' @param pdigits (default: 5) number of digits for rounding p-values
#' @examples
#' rpfmt( cbind( 'estimate' = rnorm(7), "p-values" = 10^c(-1,-2,-3,-4,-5,-6,-7)))
#' @export
pfmt <- function(x, digits = 5, scientific = FALSE) {
  x <- format(xx <- round(x, digits), scientific = scientific)
  x[as.double(xx) == 0] <- paste(c("<.", rep("0", digits -
                                               1), "1"), collapse = "")
  x
}
