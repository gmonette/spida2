#' Remove names
#'
#' Anonymize a vector or list
#'
#' Useful is a list has names but will be passed as the argument list to \code{\link{do.call}} 
#' if the names are not intended to be parameter names.
#' 
#' @param x object that may have names
#' @return x without names
#' @seealso \code{\link{auto}}
#' @export
#' @rdname anon
anon <- function(x) {
  names(x) <- NULL
  x
}

