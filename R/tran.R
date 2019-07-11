###
### tran, tr and map
###

#' Translate elements of a vector
#'
#' This should return something like 'to' whenever all of 'x'
#' is in 'from'
#' Otherwise, if some of x is not modified then it needs to
#' either to change those to NAs or to leave them as is.
#' For a strict function, we can simply index to which will
#' leave it unchanged
#' which conflicts with my versions of tr, ergo this is called tran
#'
#' @param from numerical, logical or character vector of values or factor levels to change 
#' @param to vector of corresponding values to change into
#' @param x numeric, logical, character vector or factor to modify
#' @param tofactor if TRUE, value will be a factor, default: \code{is.factor(x)}
#' @examples
#' tran(FALSE, TRUE, c(FALSE, TRUE, NA))
#' tran(c(FALSE, NA), c(TRUE, FALSE), c(FALSE, TRUE, NA))
#' tran(c(FALSE, NA), c(TRUE, FALSE), c(FALSE, TRUE, NA)) %>% mode
#' # numeric
#' tran(c(1,2), c(11,12), 1:5)
#' tran(NA, 99, c(1:5, NA))
#' tran(NA, 99, c(1:5, NA)) %>% mode
#' # mixed
#' tran(c(1,2,10), c('a','b', 'z'), 1:5)
#' # factor
#' tran(letters, 1:26, factor(c('a','b','z')) )
#' tran('a', 1, factor(c('a','b','z')) )
#' @export
tran <- function(from, to, x, tofactor = is.factor(x)) {
	if(is.factor(from)) from <- as.character(from)
	if(is.factor(to)) to <- as.character(to)
	to <- rep(to, length=length(from))
	ret <- x
	if( is.factor(x) ) {
		ret <- as.character(ret)
		levs <- levels(x)
	}
	ret <- c(to,unique(ret)) [ match( ret, c(from, unique(ret)))]
	if (tofactor) {
		if(is.factor(x)) tolevs <- tran(from,to,levs)
		else tolevs <- to
		tolevs <- c(tolevs,unique(ret))
		tolevs <- intersect( tolevs, unique( ret ))
		ret <- factor(ret,levels = tolevs)
	}
  ret
}
#' Translate elements of a vector
#'
#' \code{tr} changes the order of arguments in \code{tran}.
#' 
#' @rdname tran
#' @export
tr <- function( x, from , to ) tran( from, to, x)

#' Apply a mapping to a vector
#'
#' A mnemonic for \code{to[match(x, from)]}. 
#' Compare with the more elaborate \code{\link{tr}} or \code{\link{tran}}.
#'
#' @param from vector of values to be matched against
#' @param to vector from which selected values are returned 
#' @param x vector of values to be matched
#' @seealso tran
#' @seealso tr
#' @examples
#' x.from <- c('A',"FromNA",'Z','B',"N")
#' x.to   <- factor(c('a',NA,'z','b',NA),levels=c('z','a','b',NA))
#' x <- c("Z","B",NA,"N","M")
#' map(x.from, x.to, x)
#' @export
map <- function( from, to, x, ...) {
  to[ match( x, from, ...) ]
}

