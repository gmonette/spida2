#' Print array using kableExtra
#' 
#' Print a higher-dimensional array with subheadings for the levels
#' defined by the third and higher dimensions of the array.
#' 
#' @param a array
#' @param caption an optional caption for the printed table
#' @param ... functions from kable extra and their arguments
#' 
#' @return hmtl code to include in a 'asis' chunk
#' @examples
#' kable_array(Titanic)
#' kable_array(Titanic, row_spec = list(0, angle = -30), 
#'   add_header_above =list(c(' '= 1,'Gender'=2,Freq = 1)),
#'   column_spec(1, bold = T, border_right = T)) 
#' \dontrun{  
#' kable_array(Titanic, row_spec = list(0, angle = -30), 
#'   add_header_above =list(c(' '= 1,'Gender'=2,Freq = 1)),
#'   column_spec(1, bold = T, border_right = T)) %>%
#'   save_kable(file = 'test.html', self_contained = T)
#' }
#' @export
kable_array <- function(a,caption = '', ...) {
  require(kableExtra)
  # ...: row_spec
  dots <- list(...)
  print(length(args))
  dd <- dim(a)
  ddn <- dimnames(a)
  ret <-
    if(length(dd) > 2) {
      ap <- aperm(a, c(1,3:length(dd),2))
      ncols <- dd[2]
      nelts <- prod(dd)
      nrows <- nelts / ncols
      dim(ap) <- c(nrows,ncols)
      dnap <- list(rep(ddn[[1]], nrows/dd[1]), ddn[[2]])
      names(dnap) <- names(ddn)[1:2]
      dimnames(ap) <- dnap
      ap
    } else a
  args <- c(ddn[-(1:2)], stringsAsFactors = F)
  gnams <- do.call(expand.grid, args)
  gnams[] <- lapply(seq_along(gnams), function(ii) {
    paste0(names(gnams)[ii],': ', gnams[[ii]])
  })
  gnams <- c(rev(gnams), sep = ' / ')
  gnams <- do.call(paste, gnams)
  gnams
  ret <- kable(ret)  
  ret <- kable_styling(ret, full_width=F)
  #  aarg <- list(ret, 0 , angle = +45)
  #  ret <- row_spec(ret, 0, angle = -45)
  #  ret <- do.call('row_spec', aarg)
  for(nn in names(dots)) {
    #   print(nn)
    aargs <- list(ret)
    aargs <- c(aargs, dots[[nn]])
    ret <- do.call(nn, aargs)
  }
  for( i in seq_along(gnams)) {
    ret <- pack_rows(ret, gnams[i], (i - 1) * dd[1] + 1, i * dd[1])
  }
  ret <- kable_styling(ret, fixed_thead = T)
  ret
}