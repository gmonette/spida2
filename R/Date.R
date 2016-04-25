##
##  Functions to complete Dates
##
#' Get the year of a Date object
#'
#' Works like 'months', 'weekdays', etc. for years.
#'
#' @param x a Date object
#' @return the year as an integer
#' @export
years <- function(x,...) UseMethod("years")
#' @describeIn years method for Date objects
#' @export
years.Date <- function(x,...) as.numeric(format(x,"%Y"))
