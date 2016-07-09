###
### agg
###

#' Create a data frame at a higher level of aggregation with a possible incidence matrix for categorical factors
#'
#' Produce a higher level data set with one row per cluster. The data set 
#' contains variables that are invariant in each cluster and, optionally, summaries
#' of other variables (means for numeric variables and blocks of variables
#' corresponding to incidence matrices for factors.
#' Adapted
#' from \code{gsummary} in the \code{nlme} package and from \code{up} in
#' the \code{spida2} package.
#'
#' @param object a data frame to be aggregated.
#' @param form a one-sided formula or a list or data frame identifying the variable(s) in \code{object}
#' that identifies clusters. e.g. ~ school/Sex to get a summary within each Sex
#' of each school.
#' @param agg  a one-sided formula or a list or data frame identifying the variable(s) in \code{object}
#' to be aggregated.
#' @param sep (default _) separator to separate variable name from variable value when
#' aggregating a factor variable.
#' @param sep.clust (default /) separator to form cluster names combining more than one
#' clustering variables.  If the separator leads to the same name for distinct
#' clusters (e.g. if var1 has levels 'a' and 'a/b' and var2 has levels 'b/c'
#' and 'c') the function produces an error and a different separator should be
#' used.
#' @param FUN (default cvar) function to be used for summaries.
#' @return a data frame with one row per value of the variable in \code{form} and 
#' aggregate variabless for each variable name in \code{agg}
#' @examples
#' # a labor force survey with individual level data
#' surv <- read.table(header = TRUE, text = "
#' year sex region status
#' 2010   M      A employed
#' 2011   F      B unemployed
#' 2012   M      C employed
#' 2010   F      A employed
#' 2011   F      B employed
#' 2012   M      C out_of_labor_force
#' 2010   F      A employed
#' 2011   F      B employed
#' 2012   M      C out_of_labor_force
#' 2010   M      A employed
#' 2011   M      B unemployed
#' 2012   M      C employed
#' 2010   M      A employed
#' 2011   M      B unemployed
#' 2012   M      C employed
#' 2010   M      A employed
#' 2011   F      B unemployed
#' 2012   F      C unemployed
#' ")
#' surv
#'
#' @author largely adapted from gsummary in Bates & Pinheiro
#' @export
agg <- function ( object, form = formula(object),
                 agg = NULL,
                 sep = '_', 
                 sep.clus = "/",
                 na.rm = TRUE,
                 ...)
{
  sel.mf <- model.frame(form , object , na.action = na.include )
  if(!is.null(agg)) agg.mf <- model.frame(agg, object, na.action = na.include)
  else return(up(object, form, sep = sep.clus, na.rm = na.rm))
  ret <- object
  for (i in seq_along(agg.mf)) {
    x <- agg.mf[[i]]
    if(is.factor(x)) {
      mat <- cvar(x, sel.mf, all = T, na.rm = na.rm, ...)
      colnames(mat) <- paste0(names(agg.mf[i]),sep,colnames(mat))
    }
    else {
      mat <- cvar(x, sel.mf, na.rm = na.rm, ...)
      mat <- data.frame(x=mat)
      names(mat) <- names(agg.mf[,i])
    }
    ret <- cbind(ret, mat)
  }
  up(ret, sel.mf, sep = sep.clus, na.rm = na.rm)
}








    