#' Levels of a factor or unique values of a vector
#'
#' Useful as an argument to lapply
#'
#' @param x vector
#' @param na.rm (default TRUE) remove NAs if x is not a factor, otherwise just return levels(x)
#' @examples
#' \dontrun{
#' head(hs)
#' fits <- lapply(Levels(hs$school), function(sid)
#'            lm(mathach ~ ses, subset(hs, school == sid)))
#' coefs <- lapply(fits, coefs)
#' df <- do.call(rbind, coefs)
#' df$school <- names(df)
#' }
#' @export
Levels <- function(x, na.rm = TRUE) {
  if(is.factor(x)) ret <- levels(x)
  else {
    if(na.rm) x <- na.omit(x)
    ret <- unique(x)
  }
  names(ret) <- ret
  ret
}