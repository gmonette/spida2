#' Grep with value = T
#'
#' Equivalent of grep(..., value = TRUE)
#'
#' Now obsolete since the addition of the same command to base R.
#' 
#' @param \dots arguments to \code{\link{grep}}
#' @seealso \code{\link{grepl}}
#' @export
grepv2 <- function(...) grep( ..., value = TRUE)
