###
### capply
###

#' Apply a function within each cluster of multilevel data
#'
#' Apply a function to each cell of a ragged array, that is to each (non-empty)
#' group of values given by a unique combination of the levels of certain
#' variables and, in contrast with \code{tapply}, return within each cell a
#' vector of the same length as the cell, which are then ordered to match the
#' corresponding positions of the cells in the input.
#'
#' \code{capply} is very similar to \code{\link[stats]{ave}} in \code{package:stats}. They
#' differ in the way they treat missing values in the clustering variables.
#' \code{ave} treats missing values as if they were legitimate clustering
#' levels while \code{capply} returns a value of NA within any cluster formed
#' by a combination of clustering variable values that includes a value of NA.
#'
#' \code{capply} extends the function of \code{tapply(x, by, FUN)[ tapply(x,
#' by) ]}. The function \code{FUN} is applied to each cell of \code{x} defined
#' by each value of \code{by}. The result in each cell is recycled to a vector
#' of the same length as the cell. These vectors are then arranged to match the
#' input \code{x}.  Thus, if the value returned within each cell is a scalar,
#' the effect of \code{capply(x, by, FUN)} is the same as \code{tapply(x, by,
#' FUN)[ tapply(x, by) ]}.  \code{capply} extends this use of \code{tapply} by
#' allowing the value returned within each cell to be a vector of the same
#' length as the cell.
#'
#' The \code{capply.formula} method allows the use of two-sided formula of the
#' form \code{x ~ a + b} or \code{cbind(x, y) ~ a + b} where the variables on
#' the left-hand side are used to create a data frame that is given as a first
#' argument to \code{FUN}. If there is a single variable on the left-hand side
#' then that variable can be treated as a vector by \code{FUN}.
#'
#' @aliases capply capply.default capply.formula cvar cvars dvar dvar.factor
#' dvar.default cvar.factor cvar.default
#' @param x a vector or data frame that provides the first argument of
#' \code{FUN}
#' @param by If \code{x} is a vector: a 'factor' of the same lenth as \code{x}
#' whose levels identify clusters.  If \code{x} is a data frame, a one-sided
#' formula that identifies the variable(s) within \code{x} to be used to
#' clusters.
#' @param FUN a function to be applied to \code{x} within each cluster.
#' \code{FUN} can return a single value, or a vector whose length is equal to
#' the number of elements in each cluster.
#' @param fmla in \code{capply.formula}, fmla is a two-sided formula as in
#' \code{\link{aggregate.formula}}. The left-hand side identifies the
#' variable(s) in \code{data} to be include in a data.frame that is clusterd using
#' the variables in the right-hand side of the formula.
#' @param \dots additional variables to be supplied to \code{FUN}
#' @return When the result in each cell is a scalar, \code{capply} can be used
#' to for multilevel analysis to produce 'contextual variables' computed within
#' subgroups of the data and expanded to a constant over elements of each
#' subgroup.
#'
#' \code{capply( x , by, FUN , ...)} where \code{x} is a vector
#'
#' is equivalent to
#'
#' \code{unsplit ( lapply ( split ( x , by ), FUN, ...), by )}
#'
#' which has the same effect as
#'
#' \code{tapply( x, by, FUN, ...) [ tapply( x, by) ]}
#'
#' if \code{FUN} returns a vector of length 1.
#'
#' If \code{FUN} returns a vector, it is recycled to the length of the input
#' value.
#'
#' When the first argument is a data frame:
#'
#' \code{capply ( dd, by, FUN, ...)}
#'
#' uses unsplit - lapply - split to apply \code{FUN} to each sub data frame. In
#' this case, \code{by} can be a formula that is evaluated in 'dd'.
#'
#' This syntax makes it easy to compute formulas involving more than one
#' variable in 'dd'. An example:
#'
#' \code{capply( dd, ~gg, function(x) with( x, mean(Var1) / mean(Var2) ) )}
#'
#' where 'Var1' and 'Var2' are numeric variables and 'gg' a grouping factor in
#' data frame 'dd'.  Or, using the \code{with} function:
#'
#' \code{capply( dd, ~gg, with , mean(Var1) / mean(Var2) )}
#'
#' \code{cvar} and \code{cvars} are intended to create contextual variables in
#' model formulas. If 'x' is numerical, \code{cvar} is equivalent to
#' \code{capply(x,id,mean)} and \code{cvars} is equivalent to
#' \code{capply(x,id,sum)}.
#'
#' If \code{x} is a factor, \code{cvar} generates the equivalent of a model
#' matrix for the factor with indicators replaced by the proportion within each
#' cluster.
#'
#' \code{dvar} is equivalent to \code{x - cvar(x,by)} and creates what is
#' commonly known as a version of 'x' that is 'centered within groups' (CWG).
#' It creates the correct matrix for a factor so that the between group
#' interpretation of the effect of \code{cvar(x,by)} is that of the 'between
#' group' or 'compositional' effect of the factor.
#' @note \code{capply} tends to be slow when there are many cells and \code{by}
#' is a factor. This may be due to the need to process all factor levels for
#' each cell. Turning \code{by} into a numeric or character vector improves
#' speed: e.g. \code{capply( x, as.numeric(by), FUN)}.
#' @examples
#' \dontrun{
#'      data( hs )
#'      head( hs )
#'
#'      # FUN returns a single value
#'      hs$ses.mean <- capply( hs$ses, hs$school, mean, na.rm = T)
#'      hs$ses.hetero <- capply ( hs$ses, hs$school, sd , na.rm = T)
#'      hs.summ <- up( hs, ~school )
#'      head( hs.summ )   # variables invariant within school
#'
#'      # FUN returns a vector
#'      # with 'x' a data frame
#'      # Note how the 'with' function provides an easy way to write use a
#'      #   formula as the '...' variable.
#'
#'      hs$minority.prop <- capply( hs, ~ school, with, mean( Minority == "Yes"))
#'
#'      # equivalently:
#'
#'      hs$minority.prop <- capply( hs$Minority, hs$school, mean)
#'
#'      # on very large data frames with many columns that are not used, the 'data frame'
#'      # version of 'capply' can be very slow in comparison with 'vector' version.
#'
#'      # In contrast with 'tapply' 'FUN' can return a vector, e.g. ranks within groups
#'
#'      hs$mathach.rank <- capply( hs, ~ school, with , rank(mathach))
#'
#'      # cvar and dvar in multilevel models
#'
#'      library( nlme )
#'      data ( hs )
#'      fit <- lme( mathach ~ Minority * Sector, hs, random = ~ 1 | school)
#'      summary ( fit )
#'
#'      fit.contextual <- lme( mathach ~ (Minority + cvar(Minority, school)) * Sector,
#'                        hs, random = ~ 1| school)
#'      summary(fit.contextual) # contextual effect of cvar(Minority)
#'
#'      fit.compositional <- lme( mathach ~ (dvar(Minority,school) + cvar(Minority, school)) * Sector,
#'                        hs, random = ~ 1| school)
#'      summary(fit.compositional) # compositional effect of cvar(Minority)
#' }
#' @export
capply <- function ( x ,... ) UseMethod("capply")
#' @describeIn capply
#' @export
capply.formula <- function(formula, data, FUN, ...) {
  # the first portion of this code is from stats:::aggregate.formula
  # to avoid evaluating formula
  FUN <- match.fun(FUN)
  if (missing(formula) || !inherits(formula, "formula"))
    stop("'formula' missing or incorrect")
  if (length(formula) != 3L)
    stop("'formula' must have both left and right hand sides")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m$... <- m$FUN <- NULL
  m[[1L]] <- as.name("model.frame")
  m$na.action <- quote(na.include)
  if (as.character(formula[[2L]] == ".")) {
    rhs <- unlist(strsplit(deparse(formula[[3L]]), " *[:+] *"))
    lhs <- sprintf("cbind(%s)", paste(setdiff(names(data),
                                              rhs), collapse = ","))
    m[[2L]][[2L]] <- parse(text = lhs)[[1L]]
  }
  mf <- eval(m, parent.frame())
  # disp( str(m))
  #  disp( dim(mf))
  # disp( str(mf))
  if (is.matrix(mf[[1L]])) {
    lhs <- as.data.frame(mf[[1L]])
    names(lhs) <- as.character(m[[2L]][[2L]])[-1L]
    ret <- capply.default(lhs, mf[-1L], FUN = FUN, ...)
  }
  else ret <- capply.default(mf[1L], mf[-1L], FUN = FUN, ...)
  ret
}
#' @export
capply.default <- function ( x, by, FUN , ...) {
  FUN <- match.fun(FUN)
  if (inherits(by,'formula')) by <- model.frame( by , x , na.action = na.include)
  if (is.character(by)) by <- factor(by)
  if (is.factor(by)) by <- as.numeric(by)
  ret <- unsplit ( lapply ( split ( x , by ), FUN, ...), by )
  if ( !is.null( dim(ret)) && length(dim(ret)) ==1) ret <- c(ret)
  ret
}

# test on large data frame
#
# zh <- data.frame( a <-factor( sample(1:1000, 100000, rep = T) ), x = rnorm(100000))
# system.time(
#       ret <- capply( zh, ~a, with, x )  #     1.11    0.01    1.35
# )
# system.time(
#       ret <- capply( zh, ~as.character(a), with, x )  #  1.32    0.02    1.54
# )
# system.time(
#       ret <- capply( zh, ~as.vector(a), with, x )  #  1.38    0.00    1.49
# )
# system.time(
#       ret <- capply(x~a, zh, with , x )  # 1.09    0.05    1.26
# )
# system.time(
#       ret <- capply(cbind(x,a)~a, zh, with , x )  #  1.22    0.04    1.45
# )
#
# zh <- data.frame( a = factor(1:10000), x = 1:10000)
# system.time(
#       ret <- capply( zh, ~a, with, x )  #     9.37    0.14    9.64
# )
# system.time(
#       ret <- capply( zh, ~as.character(a), with, x )  #   10.44    0.06   10.60
# )
# system.time(
#       ret <- capply( zh, ~as.vector(a), with, x )  #   11.77    0.13   12.11
# )
# system.time(
#       ret <- capply(x~a, zh, with , x )  # 4.36    0.07    4.63
# )
# system.time(
#       ret <- capply(cbind(x,a)~a, zh, with , x )  #  5.04    0.06    5.36
# )
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
#' Create a contextual variable for regression
#'
#' cvar and dvar are designed to be used in regression formulas to
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
#' @param na.rm (default TRUE) whether to drop missing values
#' @export
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
cvar <- function( x, id , all, na.rm , ... ) UseMethod("cvar")
#' @describeIn cvar
#' @export
cvar.factor <- function(x, id, all = FALSE, na.rm = TRUE, ... ) {
  if(all) mat <- contrasts(x, contrasts = FALSE) [ x,]
  else mat <- contrasts(x) [x, ]
  ret <- cvar(mat, id, na.rm = na.rm, ...)
  colnames(ret) <- colnames(mat)
  ret
}
#' @describeIn cvar
#' @export
cvar.default <- function( x, id, all , na.rm = TRUE, ... ) {
  if ( is.matrix (x) ) {
    if ( dim(x)[2] == 1) return( cvar( x[,], id, na.rm = na.rm, ...))
    else {
      ret <-  cbind( cvar(x[,1], id, na.rm = na.rm, ...), cvar(x[,-1],id, na.rm = na.rm, ...))
      colnames(ret) <- colnames(x)
      return( ret )
    }
  } else {
    capply( x, id, mean, na.rm = na.rm)
  }
}
#' @describeIn cvar
#' @export
dvar <- function( x, id , all , na.rm , ... ) {
  help = "
  dvar: produces group mean centering: x - cvar(x, id)
  See 'cvar'
  "
  UseMethod("dvar")
}
#' @describeIn cvar
#' @export
dvar.factor <- function( x, id, all = FALSE, na.rm = TRUE, ... ) {
  if(all) mat <- contrasts( x, contrasts= FALSE) [ x,]
  else mat  <- contrasts( x ) [ x,]
  ret <- mat - cvar(mat, id, all = all, na.rm = na.rm, ...)
  colnames(ret) <- colnames(mat)
  ret
}
#' @describeIn cvar
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


#' Transform NAs to 0
#'
#' @param x
#' @export
na20 <- function(x) {
  x[is.na(x)] <- 0
  x
}



