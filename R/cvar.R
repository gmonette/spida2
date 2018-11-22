##
##
##   cvar: V0.1 August 15, 2006
##   Creating contextual variables for categorical variables
##
##   cvar is designed to create contextual variables
##   for factors, as well as for numerical variables.
##   If a factor has g levels, convar will create a
##   matrix with g-1 columns each of which is the within group
##   mean value of the correponding column of the "contrast"
##   matrix for the factor.
##
##   Modified 2018-11-21: added freq argument to up to summarize factor using frequencies
##
#' Create a contextual variable for regression
#'
#' cvar and  See also \code{\link{dvar}} are 
#' designed to be used in regression formulas to
#' create a contextual mean of a cluster-varying variable and
#' a 'centered-within-groups' version.
#'
#' @param x variable to be centered or residualized within groups. If x is a
#' factor, cvar and dvar return matrices whose columns are named consistently
#' with the names of coding variables for factors.
#' @param id identifies clusters
#' @param all (default FALSE) if TRUE cvar.factor returns the columns means of
#' an incidence matrix including the first level. Otherwise, the first level is
#' dropped for use in a linear model.
#' @param FUN (default mean) function to be applied to create contextual values
#' @param na.rm (default TRUE) whether to drop missing values
#' @examples
#' \dontrun{
#' dd <- data.frame(x= 1:100, id = rep( LETTERS[1:10], each = 10))
#' dd$a <- factor(sample( c('a','b','c'), 100, replace = T))
#' dd$y <- dd$x + rep(rnorm(10), each = 10) + rnorm(100) + as.numeric(dd$a)
#' library(nlme)
#' fit <- lme( y ~ x + cvar(x,id), dd, random = ~ 1 + dvar(x,id) | id)
#' anova( fit , type = 'm')
#' # The output of 'anova' can be used to test whether a contextual variable
#' # should be included in the model
#' }
#' @export
cvar <- function( x, id , all, na.rm , FUN = mean, ... ) UseMethod("cvar")
#' @describeIn cvar method for class 'factor'
#' @export
cvar.factor <- function(x, id, all = FALSE, na.rm = TRUE, FUN = mean, ... ) {
  if(all) mat <- contrasts(x, contrasts = FALSE) [ x,]
  else mat <- contrasts(x) [x, ]
  ret <- cvar(mat, id, na.rm = na.rm, FUN = FUN, ...)
  colnames(ret) <- colnames(mat)
  ret
}
#' @describeIn cvar default method
#' @export
cvar.default <- function( x, id, all , na.rm = TRUE, FUN = mean, ... ) {
  if ( is.matrix (x) ) {
    if ( dim(x)[2] == 1) return( cvar( x[,], id, na.rm = na.rm, FUN = FUN, ...))
    else {
      ret <-  cbind( cvar(x[,1], id, na.rm = na.rm, FUN = FUN, ...), cvar(x[,-1],id, na.rm = na.rm, FUN = FUN, ...))
      colnames(ret) <- colnames(x)
      return( ret )
    }
  } else {
    capply( x, id, FUN = FUN , na.rm = na.rm, ...)
  }
}
#' Create a centered-within-groups variable for regression
#'
#' cvar and dvar are designed to be used in regression formulas to
#' create a contextual mean of a cluster-varying variable and
#' a 'centered-within-groups' version. See also \code{\link{cvar}}.
#'
#' @param x variable to be centered or residualized within groups. If x is a
#' factor, cvar and dvar return matrices whose columns are named consistently
#' with the names of coding variables for factors.
#' @param id identifies clusters
#' @param all (default FALSE) if TRUE cvar.factor returns the columns means of
#' an incidence matrix including the first level. Otherwise, the first level is
#' dropped for use in a linear model.
#' @param na.rm (default TRUE) whether to drop missing values
#' @examples
#' \dontrun{
#' dd <- data.frame(x= 1:100, id = rep( LETTERS[1:10], each = 10))
#' dd$a <- factor(sample( c('a','b','c'), 100, replace = T))
#' dd$y <- dd$x + rep(rnorm(10), each = 10) + rnorm(100) + as.numeric(dd$a)
#' library(nlme)
#' fit <- lme( y ~ x + cvar(x,id), dd, random = ~ 1 + dvar(x,id) | id)
#' anova( fit , type = 'm')
#' # The output of 'anova' can be used to test whether a contextual variable
#' # should be included in the model
#' }
#' @export
dvar <- function( x, id , all , na.rm , ... ) {
  help = "
  dvar: produces group mean centering: x - cvar(x, id)
  See 'cvar'
  "
  UseMethod("dvar")
}
#' @describeIn dvar method for class 'factor'
#' @export
dvar.factor <- function( x, id, all = FALSE, na.rm = TRUE, ... ) {
  if(all) mat <- contrasts( x, contrasts= FALSE) [ x,]
  else mat  <- contrasts( x ) [ x,]
  ret <- mat - cvar(mat, id, all = all, na.rm = na.rm, ...)
  colnames(ret) <- colnames(mat)
  ret
}
#' @describeIn dvar default method
#' @export
dvar.default <- function( x, id, all, na.rm = TRUE, ... ) {
  if ( is.matrix (x) ) {
    if ( dim(x)[2] == 1) return( dvar( x[,], id, na.rm = na.rm,...))
    else {
      ret <-  cbind( dvar(x[,1], id, na.rm = na.rm, ...), dvar(x[,-1], id, na.rm = na.rm, ...))
      colnames(ret) <- colnames(x)
      return( ret )
    }
  } else {
    x - capply( x, id, mean, na.rm = na.rm)
  }
}

# ##
# ##  sum
# ##
#' @export
cvars <- function(  x, by, ...) {
   if ( length(x) == 1 && x == 1) {
     n <- nrow(as.data.frame(by))
     capply( rep(1,n), by, sum)
   } else {
     capply( x, by, sum, ...)
   }
}

