#' Wrapper for read_excel in 'readxl' package
#' 
#' read_excel, in the 'readxl' package, guesses the type of each column based on
#' some initial number of values.  If non-numeric 
#' values first occur quite far into the spread sheet, read_excel
#' will misguess the type and the non-numeric values will
#' be returned as missing when coerced to numeric.
#' Using the argument "col_types = 'text'" in read_excel does
#' not work, since read_excel expects the length of 'col_types'
#' to equal the number of columns. 'Read_excel' first
#' counts the number number of columns and recycles the
#' 'col_types' argument so it has the right length.
#' 
#' @param path path to xls/xlsx file
#' @param sheet string or integer to identify sheet to read. Defaults to first sheet.
#' @param col_names Either TRUE to use first row as
#'        as column names, FALSE to number sequentially from
#'        X1 to Xn, or a character vector giving a name
#'        for each column.
#' @param col_types Either NULL to guess from the spreadsheet
#'        or a character vector containing 'blank', 'numeric',
#'        'date' or 'text' which is recycled to match the 
#'        number of columns (this behavior is different from
#'        that of 'read_excel' that does not recycle)
#' @param na Missing value. By default readxl converts blank cells to missing data. Set this value if you have used a sentinel value for missing values.
#' @param skip Number of rows to skip before reading any data.  
#' @param stringsAsFactors defaults to 'default.stringsAsFactors()'
#' @return a data frame. The value returned by readxl::read_excel 
#' @export
Read_excel <- function(path, sheet = 1, col_types = NULL, 
                       stringsAsFactors = default.stringsAsFactors(),
                       ...)
{
  library("readxl")
  if(is.null(col_types)) ret <- read_excel(path, sheet = sheet, ...)
  else {
    num.columns <- length(readxl:::xlsx_col_types(path, sheet = sheet, n = 1))
    ret <- readxl::read_excel(path, sheet = sheet,
                     col_types = rep(col_types, num.columns),...)
  }
  as.data.frame(as.list(ret),stringsAsFactors = default.stringsAsFactors())
}
