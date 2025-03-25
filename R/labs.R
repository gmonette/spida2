#' Assign, extract and print labels
#'
#' Assigns, extracts and prints labels for various objects
#'
#' @param x an object
#' @param \dots additional arguments (currently not used)
#' @export
labs <- function(x,...) UseMethod("labs")
#' @rdname labs
#' @export
labs.default <- function(x,...) names(dimnames(x))
#' @rdname labs
#' @export
labs.data.frame.lab <- function( x ,...) attr(x,"labs")
#' @rdname labs
#' @export
"labs<-" <- function(x, value, ...) UseMethod("labs<-")
#' @param digits to print
#' @param quote (default FALSE) put quotes around strings
#' @param right (default TRUE) justfication
#' @rdname labs
#' @export
print.data.frame.lab <-
  function (x, ..., digits = NULL, quote = FALSE, right = TRUE)
  {
    labs <- attributes(x)$labs
    if (length(x) == 0) {
      cat("NULL data frame with", length(row.names(x)), "rows\n")
    }
    else if (length(row.names(x)) == 0) {
      print.default(names(x), quote = FALSE)
      cat("<0 rows> (or 0-length row.names)\n")
    }
    else {
      mat <- as.matrix(format.data.frame(x, digits = digits,
                                         na.encode = FALSE))
      labs <- c(labs,"","")
      labs <- labs[1:2]
      names(dimnames(mat)) <- labs
      print(mat , ..., quote = quote, right = right)
    }
    invisible(x)
  }
#' @rdname labs
#' @export
"[.data.frame.lab" <- function(x, ...){
  lab <- labs(x)
  ret <- get("[.data.frame")(x,...)
  if( inherits(ret, "data.frame")) labs(ret) <- lab
  ret
}
#' @rdname labs
#' @export
"labs<-.data.frame" <- function( x, value, ... ) {
  value <- c( value, "", "") [ 1:2 ]
  attr(x,"labs") <- value
  if( !inherits(x,"data.frame.lab")) class(x) <- c( "data.frame.lab", class(x))
  x
}
#' @rdname labs
#' @export
"labs<-.default" <- function(x, value, ...) {
  nd <- length(dim(x))
  value <- c( value, rep("",nd))[1:nd]
  names(dimnames(x)) <- value
  x
}
