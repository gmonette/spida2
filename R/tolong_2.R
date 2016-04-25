#' Long data set from wide data set
#'
#' Creates a long data set from a wide data set identifying 'time'-varying
#' variables by a separator
#'
#'
#' @param data a data frame in 'wide' form.
#' @param sep (default: '_') the separator that identifies time-varying
#' variables separating the 'root' containing the name of the variable in the
#' long data set from the 'occasion' name. e.g. the separator '_' in the
#' variable names 'y_1', 'y_2', separate the 'root', 'y', from the occasions,
#' '1' and '2'. Default: '_'.
#' @param expand (default TRUE) If different root names have inconsistent
#' occasion names, should variable names be created (with NAs as values) so all
#' root names have consistent occasion names.
#' @param timevar (default: 'time') the variable containing the occasion names
#' in the long file. Typically, this variable name should not exist in 'data'.
#' If these is a conflict, an integer is appended to the name to create a new
#' name.
#' @param idvar (default: 'id') the variable identifying the row in the wide
#' file from which each row in the long file is obtained. The variable should
#' have distinct values in the wide file. If the variable does not exist in the
#' wide file (or if it exists but does not have distinct values) a new variable
#' is created with values equal to the row number in the wide frame.
#' @param safe_sep is used to create unique variable name combinations when
#' expanding time-varying variables to have consistent occasion names.
#' @param ... Other arguments are passed to \link{list("stats::reshape")}.
#' @return a data frame in long form in which each 'root' of a time-varying
#' variables is a variable and the occasions are the values of the variable
#' with name 'timevar'.
#' @examples
#'
#' \dontrun{
#' dd <- data.frame( y.a = 1:3, y.b = 1:3, x.a= 1:3, time = 1:3,
#'                   x.b = 11:13, x.c = 21:23, id = c('a','a','b'))
#' tolong_2(dd, sep = '.')
#' tolong_2(dd, sep = '.', timevar = "type", idvar = 'patient')
#' }
#'
#' @export
tolong_2 <-
function (data, sep = "_", expand = TRUE,
          timevar = 'time',
          idvar = 'id',
          safe_sep = "#%@!",
          ...)
{
  if (any(timevar == names(data))) {
    i <- 1
    new_timevar <- paste(timevar, i, sep = '')
    while( any(new_timevar == names(data))) {
        i <- i + 1
        new_timevar <- paste(timevar, i, sep = '')
    }
    warning(paste("Variable '",new_timevar, "'used to index occasions", sep = ""))
    timevar <- new_timevar
  }
  if (any(idvar == names(data))) {
    # check if unique
    if ( length(unique(data[[idvar]])) < nrow(data)) {
      i <- 1
      new_idvar <- paste(idvar, i, sep = '')
      while( any( new_idvar == names(data))) {
        i <- i + 1
        new_idvar <- paste(idvar, i, sep = '')
      }
      warning(paste("Variable '", new_idvar, "' used as id for original rows", sep = ""))
      idvar <- new_idvar
    }
  } else {
    data[[idvar]] <- 1:nrow(data)
  }
  if (expand) {
    namessafe <- sub(sep, safe_sep, names(data), fixed = TRUE)
    varnames <- grep(safe_sep, namessafe, value = TRUE, fixed = TRUE)
    names <- unique(sub(paste(safe_sep, ".*$", sep = ""),
                        "", varnames))
    times <- unique(sub(paste("^.*", safe_sep, sep = ""),
                        "", varnames))
    allnames <- paste(rep(names, each = length(times)), sep,
                      rep(times, length(names)), sep = "")
    z <- data
    for (nn in allnames) {
      z[[nn]] <- if (is.null(data[[nn]]))
        NA
      else data[[nn]]
    }
    data <- z
  }
  namessafe <- sub(sep, safe_sep, names(data), fixed = TRUE)
  namestimes <- sub(paste("^.*", safe_sep, sep = ""), "ZZZZZ",
                    namessafe)
  ord <- order(namestimes)
  data <- data[, ord]
  stats::reshape(data, direction = "long",
          sep = sep, idvar = idvar, timevar = timevar,
          varying = grep(sep, names(data), fixed = TRUE), ...)
}
