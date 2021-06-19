#' Substitute last occurrence of a pattern in a string
#' 
#' @param pattern character string containing a regular expression whose
#'        whose last match in 'x' will be replaced
#' @param replacement character string that is a replacement for the
#'        matched 'pattern' in 'x'
#' @param x a character vector where matches are sought.
#' @param ... arguments passed to \code{\link{gsub}}.
#' 
#' @examples
#' string <- 'a_b_c_D'
#' sublast('_','__',string)
#' dd <- data.frame(id = 1:3, X_a_1 = 1:3, X_a_2 = 1:3, X_b_1 = 1:3, X_b_2 = 1:3)
#' dd
#' names(dd) <- sublast('_','__',names(dd))
#' tolong(dd, sep = '__')
#' tolong(dd, sep = '__') %>% tolong(sep = '_', idvar = 'id2', timevar = 'time2')
#' 
#' @export
sublast <- function(pattern, replacement, x, ...) {
  pat <- paste0('(', pattern, ')(?!.*\\1)')
  # disp(pat)
  sub(pat, replacement, x, ..., perl = TRUE)
}
#' Change NAs to FALSE
#' 
#' @param x vector, possibly with NAs
#' @export
na2f <- function(x) {
  x[is.na(x)] <- FALSE
  x
}
#' Change NAs to TRUE
#' 
#' @param x vector, possibly with NAs
#' @export
na2t <- function(x) {
  x[is.na(x)] <- TRUE
  x
}
# Pipe from magrittr
#
# Removed because of conflict with pipe defined in tidyverse
# @importFrom magrittr %>%
# @export
# magrittr::`%>%`

#' Transform NAs to 0
#'
#' @param x vector, possibly with NAs
#' @export
na20 <- function(x) {
  x[is.na(x)] <- 0
  x
}
#' Directory or filename of active R script
#' 
#' @param dir (default: TRUE) return the directory of the current R script, otherwise the full file name
#' 
#' @export
here <- function(dir = FALSE) {
  path <- rstudioapi::getActiveDocumentContext()$path
  if(dir) path <- dirname(path)
  path
}
#' Set working directory to directory of active R script
#' 
#' Set working directory to directory of active R script
#' if script is not being knitted
#' 
#' @export
setwd_here <- function() {
  path <- here(TRUE)
  if(is.null(knitr::opts_knit$get('output.dir'))) setwd(path)
  invisible(here)
}
#' Vectorized ifelse with multiple conditions
#' 
#' Avoids nested ifelse statements when the action depends on
#' the value of a variable
#' 
#' @param .select. a variable whose values determine the 
#'        argument to be used
#' @param \dots named arguments and one possibly unnamed argument. 
#' 
#' 
#' @details 
#' Each argument in \dots evaluates
#' to a vector whose value is returned where the name of the 
#' argument matches a value of \code{.select.}. 
#' 
#' The vectors in \dots are combined into a matrix with
#' \code{\link{cbind}} and the names of the arguments
#' are used as values of \code{.select.} to select which
#' vector value is returned.  See the examples. 
#' 
#' If there is an unnamed argument, its value is used
#' as a value in \code{.select.} is not matched by
#' an argument name.
#' 
#' See an alternative: \code{\link[dplyr]{case_when}}
#'
#' @examples
#' x <- c(letters[1:4],NA)
#' case(x, a = 'was an a', b = 'was a b', z = 'was a z')
#' case(x, a = 'was an a', x) # x is returned as default
#' # numerical 'select' is coerced to character
#' case(1:4, '1' = 'was a 1', '2' = 'was a 2')
#' 
#' location <- c('England','England','France','France',
#'      'Germany','Spain')
#' xvar <- c('yes','no','non','oui','nein','si')
#' case(location,
#'    'England' = tr(xvar, c('yes','no'), c(1,0)),
#'    'France'  = tr(xvar, c('oui','non'), c(1,0)),
#'    'Germany' = tr(xvar, c('nein','ja'), c(0,1)))
#' case(location,
#'    'England' = tr(xvar, c('yes','no'), c(1,0)),
#'    'France'  = tr(xvar, c('oui','non'), c(1,0)),
#'    'Germany' = tr(xvar, c('nein','ja'), c(0,1)),
#'    xvar) 
#' case(location,
#'    'England' = tr(xvar, c('yes','no'), c(1,0)),
#'    'France'  = tr(xvar, c('oui','non'), c(1,0)),
#'    'Germany' = tr(xvar, c('nein','ja'), c(0,1)),
#'    'no match')
#' @export
case <- function(.select., ...) {
  nas <- is.na(.select.)
  replace <- list(...)
  levels <- names(replace)
  # if "" is in levels, i.e. if there is an unnamed argument
  # then this is the default for non-matches
  # otherwise non-matches return NA
  which <- match(as.character(.select.), levels)
  if(length(default <- grep("^$", levels))) which[is.na(which)] <- default
  # But NAs in select nevertheless return NAs
  which[nas] <- NA
  what <- do.call(cbind, replace)
  what[cbind(1:nrow(what), which)]
}
#' Left Cholesky factor
#' 
#' Decomposes positive-definite G = L'L where L is lower-triangular.
#' 
#' In R, \code{\link{chol}} returns a upper-triangular matrix \code{R}
#' such that G = R'R. \code{lchol} return a lower-triangular matrix.
#' 
#' @param x a positive-definite matrix
#' @examples
#' mm <- cbind( c(8,2,1), c(2,10,2), c(1,2,5))
#' mm
#' chol(mm)
#' lchol(mm)
#' crossprod(chol(mm))
#' t(chol(mm)) %*% chol(mm)
#' crossprod(lchol(mm))
#' t(lchol(mm)) %*% lchol(mm)
#' @export
lchol <- function(x) {
  rind <- rev(1:nrow(x))
  xret <- x[rind,][,rind]
  ret <- chol(xret)
  ret[rind,][,rind]
}
#' Vovk-Sellke Maximum p-Ratio
#' 
#' Calculates the Vovk-Sellke Maximum p-Ratio
#' 
#' @param p p-values
#' @aliases vovk sellke
#' @export 
vs <- function(p) {
  -1/(exp(1) * p * log(p))
}
#' Decomposes positive-definite G = L'L where L is lower-triangular.
#' 
#' In R, \code{\link{chol}} returns a upper-triangular matrix \code{R}
#' such that G = R'R. \code{lchol} return a lower-triangular matrix.
#' 
#' @param x a positive-definite matrix
#' @examples
#' mm <- cbind( c(8,2,1), c(2,10,2), c(1,2,5))
#' mm
#' chol(mm)
#' lchol(mm)
#' crossprod(chol(mm))
#' t(chol(mm)) %*% chol(mm)
#' crossprod(lchol(mm))
#' t(lchol(mm)) %*% lchol(mm)
#' @export
lchol <- function(x) {
  rind <- rev(1:nrow(x))
  xret <- x[rind,][,rind]
  ret <- chol(xret)
  ret[rind,][,rind]
}
