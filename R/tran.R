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
#' @param from
#' @param to
#' @param x
#' @param tofactor
#' @export
tran <- function( from, to, x, tofactor = is.factor(x)) {
	if(is.factor(from)) from <- as.character(from)
	if(is.factor(to)) to <- as.character(to)
	to <- rep(to, length=length(from))
	ret <- x
	if( is.factor(x) ) {
		ret <- as.character(ret)
		levs <- levels(x)
	}
	ret <- c(to,unique(ret)) [ match( ret, c(from, unique(ret)))]
	# DEBUG: print(ret)
	if (tofactor) {
		if(is.factor(x)) tolevs <- tran(from,to,levs)
		else tolevs <- to
		tolevs <- c(tolevs,unique(ret))
		tolevs <- intersect( tolevs, unique( ret ))
		ret <- factor(ret,levels = tolevs)
	}
  ret
}
#' Translate elements of a character vector
#'
#' Changes the order of arguments in \code{\link{tran}}.
#'
#' @param x
#' @param from
#' @param to
#' @export
tr <- function( x, from , to ) tran( from, to, x)

#' Apply a mapping to a vector
#'
#' A mnemonic for \code{to[match(x, from)]}. Compare with the more elaborate \code{\link[yscs]{tran}}.
#'
#' @param from
#' @param to
#' @param x
#' @examples
#' x.from <- c('A',"FromNA",'Z','B',"N")
#' x.to   <- factor(c('a',NA,'z','b',NA),levels=c('z','a','b',NA))
#' x <- c("Z","B",NA,"N","M")
#' map(x.from, x.to, x)
#' @export
map <- function( from, to, x ) {
  to[ match( x, from) ]
}

