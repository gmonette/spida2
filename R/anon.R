#' Remove names
#'
#' Anonymize a vector or list
#'
#' Useful is a list has names but will be passed as the argument list to \code{\link{do.call}} 
#' if the names are not intended to be parameter names.
#' 
#' @param x object that may have names
#' @return x without names
#' @seealso \code{\link[spida2]{name}}
#' @seealso \code{\link[spida2]{auto}}
#' @export
#' @rdname anon
anon <- function(x) {
  names(x) <- NULL
  x
}
#' test1 to see how to combine different function
#' 
#' test1 second paragraph
#' 
#' test1 third paragraph
#' @section test1:
#' special info for test1
#' 
#' @param x the x parameter of test1
#' @param y the y parameter of test1
#'
test1 <- function(x, y) {
  x + y
}
#' test2 to see how to combine different function -- won't see this
#' 
#' test2 second paragraph -- goes in description
#' 
#' test2 third paragraph -- goes in details
#' 
#' @param x the x parameter of test2
#' @param y the y parameter of test2
#' @param z the z parameter of test2
#'
#' @section test2:
#' special info for test2
#' @rdname test1
test2 <- function(x, y, z) {
  x + y
}


