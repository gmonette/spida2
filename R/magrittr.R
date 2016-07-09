###
### Extra pipes for magrittr
###

#' Pipe for magrittr that passes LHS invisibly
#'
#' @param lhs
#' @param rhs
#' @examples
#' data(hs)
#' hs %T% dim  %T% names  %>% head 
#' @export

`%T%` <- function(lhs,rhs) {
  print(lhs %>% rhs)
  cat("----------\n")
  invisible(lhs)
}