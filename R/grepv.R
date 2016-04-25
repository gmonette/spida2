#' Grep with value = T
#'
#' Equivalent of grep(..., value = TRUE)
#'
#' @param \dots
#' @seealso \code{\link{grepl}}
#' @export
grepv <- function(...) grep( ..., value = TRUE)
