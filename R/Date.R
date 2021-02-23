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
#' Extract month as a factor
#' 
#' @param date an object that can be converted to a 'Date'
#' @export
Months <- function(date,...) {
  date <- as.Date(date)
  factor(months(date), levels = c('January','February','March','April','May','June','July','August','September','October',
                                  'November','December'))
}
if(FALSE) {
  Months('2021_02_23')
}
#' Extract day of the week as a factor
#' 
#' @param date an object that can be converted to a 'Date'
#' @param first day of the week, default: "Sunday"
#' @export
Weekdays <- function(date,
                     first = 'Sunday',...){
  data <- as.Date(date)
  levs <-  c('Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday')
  if(is.character(first)) first <- match(first, levs)
  
  factor(weekdays(date), levels = c('Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'))
}
#' 