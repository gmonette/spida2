#' Cartesian product of model predictors
#' 
#' Cartesian product of levels of predictor values in a model with
#' the option to specify desired values for some terms of the
#' model
#' 
#' @param fit a model for which \{code\{getD}} returns the 
#'        model frame, assuming that the response is 
#'        in the first column
#' @param ...  predictors in the model for which a different
#'        set of values will be used.
#' @examples
#' fit <- lm(mathach ~ ses * Sector, hs)
#' dd <- makegrid(fit, ses = seq(-2,2,.1))
#' library(car)
#' fit <- lm(income ~ type * education, Prestige)
#' dd <- makegrid(fit, type = 'prof', education = 6:18)
#' # note that dd$type will have the correct form, 
#' # i.e. a factor with the same three levels as those in 
#' # the original data set.
#' levels(dd$type)
#'        
#' @export     
makegrid <- function(fit, ...) {
  dots <- list(...)
  dd <- getD(fit)[,-1] # remove response
  dd <- as.list(dd)
  dd <- lapply(dd, unique)
  if(length(dots)) {
    for(nn in names(dots)){
      dd[nn] <- dots[nn] 
    }
  }
  dd <- c(dd, stringsAsFactors = FALSE)
  dd <- do.call(expand.grid, dd)
  dd <- getD(fit, add = dd )
  dd <- subset(dd, .source == 'add')
  dd$.source <- NULL
  dd
}