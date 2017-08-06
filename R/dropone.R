#' Drop each row or cluster from a fitted object
#' 
#' Drop each row or each cluster in a fitted object and 
#' collect the resulting linear estimates in a matrix.
#' 
#' @param fit a fitted object with an \code{\link{update}} method
#'        and a \code{\link{getData}} method.
#' @param form a formula evaluated in the data frame for 'fit'
#'        that defines clusters to be dropped one at a time. If NULL, the 
#'        default, each row is dropped successively.
#' @param FUN function used to extract coefficients. Default:
#'        \code{\link{fixef}} for class 'lme' and 
#'        \code{\link{coef}} otherwise. 
#' @param data  data frame in which to evaluate the refitted object. Default:
#'        'getData(fit)'.
#' @export
dropone <- function(fit, form = NULL, 
                    FUN = if(inherits(fit, 'lme')) fixef else coef, 
                    data = getData(fit),...) {
  if(is.null(form)) {
    by <- factor(1:nrow(data))
    data$by <- by
    dframe <- data
  } else {
    by <- model.frame(form, data)
    by <- do.call(paste, c(by, sep='/'))
    data$by <- factor(by)  # bad if by already exists in data
    dframe <- up(data, form)
  }
  values <- dframe$by
  names(values) <- values
  ret <- lapply(values, function(v) {
    ret <- try(update(fit, data = data[by != v, ,drop=FALSE],...))
    if(class(ret) == 'try-error') NA else FUN(ret)
  })
  ret <- do.call(rbind,ret)
  colnames(ret) <- paste0('b_',colnames(ret))
  cbind(dframe, ret)
}
