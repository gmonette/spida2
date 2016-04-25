#' Displaying p-values with estimates
#'
#' Takes a 2 column matrix or data frame in which the first column has estimated
#' values and the second has p-values and formats them into a single column
#' with p-values in parentheses
#'
#' @param x a n x 2 matrix or data frame
#' @param digits (default: 5) number of digits for rounding coefficient
#' @param pdigits (default: 5) number of digits for rounding p-values
#' @examples
#' rpfmt( cbind( 'estimate' = rnorm(7), "p-values" = 10^c(-1,-2,-3,-4,-5,-6,-7)))
#' @export
rpfmt <- function(x, ...) UseMethod('rpfmt')
#' @rdname rpfmt
#' @export
rpfmt.default <- function(x, digits = 3, pdigits = 5) {
  # x is a n x 2*k matrix or data frame with
  # values in first column and p-values in second
  k <- ncol(x)/2
  if( (k %% 1) != 0) stop('x must have an even number of columns')
  digits <- rep(digits, length.out = k)
  pdigits <- rep(pdigits, length.out = k)
  ret <- matrix("", nrow(x), k)
  for( i in 1:k) {
    vals <- rnd(x[,2*(i-1)+1], digits = digits[i])
    ps <- pfmt(x[,2*(i-1)+2], digits = pdigits[i])
    co <- cbind(paste0(vals," (",ps,")"))
    ret[,i] <- co
  }
  colnames(ret) <- colnames(x)[seq(1,ncol(x),2)]
  rownames(ret) <- rownames(x)
  ret
}
#' @rdname rpfmt
#' @details works on first element of a \code{\link{wald}} list
#' @export
rpfmt.wald <- function(w, ...) {
  w <- w[[1]][[2]]
  # disp(w)
  x <- w[,c("Estimate",'p-value')]
  rpfmt(x, ...)
}



