#' Easy pch characters
#'
#' Fifteen pch symbols and six line types can be generated
#' mnemonic names.
#'
#' @param shape 'square', 'circle', 'triangle', 'diamond', 'ltriangle' (or 1, 2, 3, 4, 5 respectively). Default: 1
#' @param fill 'empty', 'solid', 'fill' (with different color) (or 1, 2, 3 respectively). Default: 1
#' @param line 'solid','dashed','dotted','dash dot','long dash','long dash dot' (or 1, 2, 3, 4, 5, 6 respectively). Default: 1
#' @return value(s) of pch or lty to produce corresponding character or line type
#' @export
pch <- function(shape = 1, fill = 1){
  shapenames <- c('square', 'circle', 'triangle', 'diamond', 'ltriangle')
  fillnames <- c('empty', 'solid', 'fill')
  if(!is.numeric(shape)) shape <- sapply(shape,
                                         function(s) grep(s,shapenames)[1])
  if(!is.numeric(fill)) fill <- sapply(fill,
                                       function(s) grep(s,fillnames)[1])
  mat <- t(rbind( c(0,1,2,5,6), c(15,16,17,18,25),c(22,21,24,23,25)))
  mat[cbind(shape,fill)]
}
#' @describeIn pch Easy line types
#' @export
lty <- function(line = 1){
  linetypes <- c('solid', 'dashed', 'dotted', 'dash dot', 'long dash', 'long dash dot')
  if(!is.numeric(line)) line <- sapply(line, function(s) grep(s,linetypes)[1])
  line
}
