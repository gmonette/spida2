##
## functions for lattice graphs
##
#' Plot '100 x log scale' data with real scale on right
#'
#' Plots axis on right side of lattice graph showing pretty selection of 
#' real values corresponding to 100 x log values in graph. The function has no 
#' parameters and is intended to be supplied as an argument to the 
#' `yscale.components` parameter in `xyplot``
#' @concepts lattice
#' @examples
#' \dontrun{
#' library(car)
#' library(lattice)
#' log100 <- function(x) 100*log(x)
#' (ob1 <- xyplot(log100(income) ~ education | type, Prestige, groups = type))
#' (ob2 <- update(ob1, ylab = expression(plain(log)[e](income) %*% 100)))
#' (ob3 <- update(ob2, 
#'               scales = list(y = list(alternating=3)), 
#'               ylab.right = 'income',
#'               yscale.components = 
#'                 yscale.components.log100real
#'    )
#' )
#' }
#' @export
yscale.components.log100real <-
  function(...)  {
    # y on 100 x log scale on left axix and real values on right axis.
    ans <- yscale.components.default(...)
    ans$right <- ans$left
    #    ans$left$labels$labels <-
    #      parse(text = sprintf("%s ~ degree * F", ans$left$labels$at))
    stretch <- function(ran, fac = 0) ran + fac * c(-1,1) * diff(ran)
    prettyRight <- pretty(n=7, exp(stretch(ans$num.limit/100)))
    ans$right$ticks$at <- 100* log(prettyRight)
    ans$right$labels$at <- 100* log(prettyRight)
    ans$right$labels$labels <- prettyRight
    #      parse(text = sprintf("%s ~ degree * C", prettyC))
    #    ans$scales <- list(y = list(xlab.right = 'AORK'))
    ans
}
