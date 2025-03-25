#' Round and format
#'
#' Round to a given number of significant digits to show variability in a vector of numbers.
#'
#' Instead of rounding to a given number of significant digits, or rounding
#' to a given significance, 'rnd' rounds so that the variability in a set of
#' number is exhibited to a given number of digits. For example, 100.1233, 100.2345, 100.3456,
#' rounded to 3 significant digits would be 100., 100., 100., but rounded to 3 significant
#' digits for variability is 100.123, 100.234, 100.345.
#'
#' @param x an object to be rounded and formatted
#' @param digits significant digits for variation in x
#' @param ... passed on to \code{\link{format}} function
#' @seealso \code{\link[spida2]{pfmt}}, \code{\link[spida2]{rpfmt}}
#' @examples
#' rnd(c(0.00001111111,0.000022222,0.000015), scientific = F)
#' rnd(c(0.00001111111,0.000022222,0.000015)+100, scientific = F)
#' round(c(0.00001111111,0.000022222,0.000015))
#' round(c(0.00001111111,0.000022222,0.000015)+100)
#' rnd(123)
#' rnd(123.45678901)
#' @export
rnd <- function(x, digits, ...) UseMethod('rnd')
#' @describeIn rnd default method
#' @export
rnd.default <- function(x, digits = 3, ..., verbose = 0) {
  if (is.numeric(x)){
    ran <- diff(range(x, na.rm = T))
    mea <- mean(x, na.rm = T)
    # show variability
    dig <- max(round(-log(ran, base = 10) + digits),0)
    if(verbose) disp(dig)
    if(is.finite(dig)) {
      x <- round(x, digits = dig)
      ret <- format(x, nsmall = dig, ...)
    } else {
      ret <- format(signif(x, digits = digits + 3),...)
    }
  } else {
    ret <- format(x, ...)
  }
  ret
}
#' @export
rnd.default <- function(x,...) x
#' @export
rnd.data.frame <- function(x, ...) {
  x[] <- lapply(x, rnd, ...)
  x
}
#' @export
rnd.numeric <- function(x,...) round(x, ...)
