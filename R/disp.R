#' Print the name and value of a variable in a function
#'
#' Useful for debugging as a substitute for printing the value of a variable in a function. It also prints the name of the variable.
#'
#' @param x value to print
#' @param head (default: deparse(substitute(x))) heading preceding printed value
#' @export
disp<-function(x, head = deparse(substitute(x)))
{
  cat("::: ", head, " :::\n")
  print(x)
  cat("======================\n")
  invisible(x)
}