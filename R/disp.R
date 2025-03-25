#' Print the name and value of a variable in a function
#'
#' Useful for debugging as a substitute for printing the value of a variable in a function. It also prints the name of the variable.
#'
#' Prints object if \code{options(verbose=TRUE)}
#' 
#' @param x value to print
#' @param head (default: deparse(substitute(x))) heading preceding printed value
#' @export
disp <- function(x, head = deparse(substitute(x)))
{
  if(isTRUE(options('verbose')[['verbose']])) {
    cat("::: ", head, " :::\n")
    print(x)
    cat("======================\n")
  }
  invisible(x)
}
