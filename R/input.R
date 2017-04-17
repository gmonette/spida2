#' Wrapper for read_excel in 'readxl' package
#' 
#' Facilitates reading excel file so each variable is read as a
#' character without automatic coercion when it fails with 'read_excel'
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
#' @param col_types Either NULL to use guess from the spreadsheet made
#'        by 'read_excel'
#'        or a character vector containing 'blank', 'numeric',
#'        'date' or 'text' which is recycled to match the 
#'        number of columns in the spreadsheet. 
#'        This behavior is different from
#'        that of 'read_excel' that does not recycle. 
#'        If \code{col_types} is equal to 'text', all variables in
#'        the returned data frame have class character and the user
#'        needs to handle conversions manually.
#'         
#' @param na Missing value. By default readxl converts blank cells to missing data. Set this value if you have used a sentinel value for missing values.
#' @param skip Number of rows to skip before reading any data.  
#' @param stringsAsFactors defaults to 'default.stringsAsFactors()'
#' @return a data frame. The value returned by readxl::read_excel 
#' @seealso \code{\link{read_excel}}
#' @examples
#' \dontrun{
#' library(spida2)
#' library(readxl)  # install from CRAN if necessary
#' url <- 'http://nross626.math.yorku.ca/MATH4939/2017/files/Read_excel_example.xlsx'
#' fname <- 'temp_file.xlsx'
#' if(!file.exists(fname)) download.file(url, fname, mode = 'wb') 
#' dd1 <- Read_excel(fname)
#' class(dd1)
#' sapply(dd1, class)
#' # Note on the POSIXct date class:
#' # - adding 1 to a "POSIXct" object adds one second
#' dd1$date_var
#' dd1$date_var + 1
#' # POSIXct objects keep track of day and time.
#' # If you only need dates,
#' # you can coerce a POSIXct object to a Date object:
#' z <- as.Date(dd1$date_var)
#' # adding 1 to a 'Date' object adds one day
#' z
#' z + 1
#' #
#' #  Reading all columns as text:
#' #
#' dd2 <- Read_excel(fname, col_types = 'text')
#' sapply(dd2, class)
#' dd2 # all characters
#' # 
#' # After fixing problematic variables:
#' #
#' # Date conversions:
#' # Internally dates are stored as a numeric variable denoting 
#' # the number of days from a 'origin'.  Different systems use
#' # different origins. Excel, unfortunately, uses 2, either
#' # '1899-12-30' or '1904-12-31'. (see https://support.microsoft.com/en-ca/help/214330/)
#' # Thus, it is safest to do a manual conversion followed by
#' # sanity checks to make sure the right origin was used.
#' # Note that '1899-12-30' is '1899-12-31' within Excel because
#' # they didn't take into account that 1900 was not a leap year
#' # in the Gregorian calendar. But in converting Excel's integer
#' # to a system that knows about the Gregorian calendar (such
#' # the POSIX standard used in Unix and R) you need to use
#' # '1899-12-30'.
#' # So, assuming the xlsx file was created in Windows:
#' (date_var <- as.Date( as.numeric(dd2$date_var), origin = '1899-12-30'))
#' # do a sanity check. If it was created on a Mac you might have to use:
#' (date_var <- as.Date( as.numeric(dd2$date_var), origin = '1903-12-31'))
#' # Numeric conversions:
#' (num_var <- as.numeric(dd2$num_var))
#' # Create a new data frame with converted variables:
#' dd_new <- with(dd2,
#'        data.frame(date_var = as.Date(as.numeric(data_var), 
#'                              origin = '1899-12-30'),
#'                   num_var = as.numeric(num_var),
#'                   char_var = as.factor(char_var))
#' dd_new                   
#' }
#' @export
Read_excel <- function(path, sheet = 1, col_types = NULL, skip = 0,
                       stringsAsFactors = default.stringsAsFactors(),
                       ...)
{
  if(!require("readxl")) stop('Install the "readxl" package to use "Read_excel"')
  if(is.null(col_types)) ret <- read_excel(path, sheet = sheet, skip = skip, ...)
  else {
    col_types[col_types == 'character'] <- 'text'
    opts <- options(warn=-1)
    num.columns <- ncol(read_excel(path, sheet = sheet, skip = skip))
    options(opts)
    ret <- read_excel(path, sheet = sheet,
                     col_types = rep_len(col_types, num.columns),
                     skip = skip, ...)
  }
  if(all(col_types == 'text')) as.data.frame(as.list(ret),
                                             stringsAsFactors = FALSE)
  else as.data.frame(as.list(ret),
                stringsAsFactors = default.stringsAsFactors())
}
