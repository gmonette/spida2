#' Contextual lag with intrapolation
#'
#' Lags within clusters based on a time vector with possibly unequal intervals
#' 
#' @param x vector to be lagged
#' @param id clustering variable
#' @param idx time variable
#' @param lag (default 1) 
#' @param at alternative to 'lag': evaluate x at a particular value of 'idx'
#' @param check (default TRUE) check for unique values of idx within each cluster
#' @seealso \code{\link[spida2]{LagI}}, \code{\link[spida2]{DiffI}}, \code{\link[spida2]{Diff}}, 
#' @return return the value of 'x' in the position 'idx - lag' in the same cluster
#' @export
Lag <- function(x, id = rep(1, length(x)), idx = 1:length(x), lag=1, at = NULL, check = T) {
  ## Takes approx. 30% of time of version Lag.0 on a small problem
  if (check) {
    comb <- paste(id,idx)
    if (any(duplicated(comb))) stop("idx not unique in each level of id")
  }
  ret <- x
  names(ret) <- paste(id,idx,sep='::')
  retnn <- if(is.null(at)) paste(id,idx - lag,sep='::')  else paste(id,at,sep="::")
  ret [ retnn ]
}
#' Contextual diff
#'
#' Diff within clusters based on a time vector with possibly unequal intervals
#' 
#' @param x vector to be lagged
#' @param id clustering variable
#' @param idx time variable
#' @param lag (default 1) 
#' @param at alternative to 'lag': evaluate x at a particular value of 'idx'
#' @param check (default TRUE) check for unique values of idx within each cluster
#' @seealso \code{\link[spida2]{LagI}}, \code{\link[spida2]{DiffI}}, \code{\link[spida2]{Diff}}, 
#' @return return the value of x - x in the position 'idx - lag' in the same cluster
#' @export
Diff <- function(xx,...) xx - Lag(xx,...)
#' Contextual lag with intrapolation
#'
#' Lags within clusters based on a time vector with possibly unequal intervals
#' 
#' @param x vector to be lagged
#' @param id clustering variable
#' @param idx time variable
#' @param lag (default 1) 
#' @param delta (default .01) 
#' @param check (default TRUE) check for unique values of idx within each cluster
#' @seealso \code{\link[spida2]{LagI}}, \code{\link[spida2]{DiffI}}, \code{\link[spida2]{Diff}}, 
#' @return return the value of 'x' in the position 'idx - lag' in the same cluster
#'         calculated by intrapolation.
#' @export
LagI <- function(x,id,time,lag=1,delta=.01,check=T) {
  # lags by intrapolating 
  # with complete data at each value of index, this does the same thing
  # as Lag
  # If values of Lag are skipped then we get linear intrapolations.
  # Note that 'delta' should be small enough so that values of x are
  # at least delta apart. However, too small a value for delta introduces
  # numerical error
  #  	
  if (check) {
    comb <- paste(id,time)
    if (any(duplicated(comb))) stop("id not unique in each level of idx")
  }
  ret <- x
  id <- as.character(id)
  names(x) <- id
  names(time) <- id
  for (nn in unique(id)){
    pos <- id == nn
    xx <- x[pos]
    tt <- time[pos]
    topred <- tt-delta
    drop <- is.na(xx)|is.na(tt)
    xxc <- xx[!drop]
    ttc <- tt[!drop]
    nl <- length(xxc)
    if ( nl > 0) {
      if ( nl > 1 ) xx.p <- approx(ttc,xxc,topred)$y
      else xx.p <- NA 
      xx.lag <- xx - lag*(xx - xx.p)/delta
      ret[pos] <- xx.lag
    }	
  }
  ret
}
#' Contextual diff with intrapolation
#'
#' Diffs within clusters based on a time vector with possibly unequal intervals
#' 
#' @param x vector to be lagged
#' @param id clustering variable
#' @param idx time variable
#' @param lag (default 1) 
#' @param at alternative to 'lag': evaluate x at a particular value of 'idx'
#' @param check (default TRUE) check for unique values of idx within each cluster
#' @seealso \code{\link[spida2]{LagI}}, \code{\link[spida2]{DiffI}}, \code{\link[spida2]{Diff}}, 
#' @return return the value of  x minus x in the position 'idx - lag' in the same cluster 
#'         using intrapolation
#' @export
DiffI <- function(xx,...) {
  xx - LagI(xx,...)
}
