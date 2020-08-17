#' Grep with value = T
#'
#' Equivalent of grep(..., value = TRUE)
#'
#' @param \dots arguments to \code{\link{grep}}
#' @seealso \code{\link{grepl}}
#' @export
grepv <- function(...) grep( ..., value = TRUE)
