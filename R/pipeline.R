## Pipeline Convenience functions ####

## Sort a data frame ####
# GM 2014 08 09



#' Order the rows of a data frame
#'
#' Order the rows of a data frame using variables identified in a formula
#' evaluated in the data frame. Convenient in a pipeline.
#'
#' @param data a data frame
#' @param form a formula identifying variables to be used in sorting or an
#' object that can be used directly in \code{data[order(form,,drop=FALSE)]}.
#' @return The formula is evaluated using \code{model.frame} and the
#' result is used as the argument of \code{order} which, in turn is used to
#' order the data frame. The ordered data frame is returned.
#' @examples
#' \dontrun{
#' require(car)
#' require(magrittr)
#' Prestige  %>% sortdf(~type+income) -> Prestige.ordered
#' }
#' @export
sortdf <- function(data, form = formula(data)) {
  if(is.fomula(form)) xx <- as.list(model.frame(form, data, na.action=NULL))
  ord <- do.call(order, xx)
  data[ ord,,drop = FALSE]
}

## Assign in a pipeline

#' Assign in a pipeline
#'
#' Same as \code{assign} with the first two arguments reversed and the default
#' position referrring to the global environment so it can be used
#' interactively in a pipeline
#'
#' @param value to be assigned
#' @param x name of the object to be assigned
#' @param pos where the object is to be assigned. Default is \code{pos=0} which
#' saves in the global environment
#' @param ...  other arguments are passed to \code{assign}
#' @return The formula is evaluated using \code{model.frame} and the
#' result is used as the argument of \code{order} which, in turn is used to
#' order the data frame. The ordered data frame is returned.
#' @examples
#' \dontrun{
#' # Using the Prestige data set in package:car and the pipeline
#' # function in package:magrittr, we can sort the data frame and save the sorted
#' # version as "Prestige.ordered" with:
#' require(magrittr)
#' require(car)
#' Prestige  %>% sortdf(~type+income)  %>% assn("Prestige.ordered")  %>%  dim
#' }
#' @export
assn <- function(value, x, pos = 0, ...) assign(x = x, value = value, pos = pos, ...)
