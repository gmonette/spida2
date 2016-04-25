





#' Centre of an object
#'
#'
#'
#'
#'
#' @param obj
#' @param \dots

#' @author Georges Monette


#' @keywords ~kwd1 ~kwd2
#' @examples
#'
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#'
#' ## The function is currently defined as
#' function (obj, ...)
#' UseMethod("center")
#'
#' @export
center <- function( obj, ... ) UseMethod("center")





#' Center of an ellipse
#'
#'
#'
#'
#'
#' @param obj
#' @param \dots

#' @author Georges Monette


#' @keywords ~kwd1 ~kwd2
#' @examples
#'
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#'
#' ## The function is currently defined as
#' function (obj, ...)
#' attr(obj, "parms")$center
#'
#' @export
center.ell <- function( obj, ...) attr(obj, 'parms') $ center


# line between two ellipsoids ( would like to extend to cylinders)




#'
#' ellipses
#'
#'
#' Generates points on the curve of osculation between the centers of two
#' families of ellipses
#'
#'
#'
#' @aliases osculant.default osculant
#' @param center1
#' ellipses.
#' @param shape1
#' first family of ellipses.
#' @param center2 of second family.
#' @param shape2
#' @param n n + 1 is the number of points to generate along the locus. \code{n
#' = 1} generates the two centers, \code{n=0} generates to point on the first
#' ellipse that lines on the locus of osculation provided the centre of the
#' second family lines outside the first ellipse.
#' @param range of values of \code{u} to use to generate points. (See the
#' algorithm in the code)

#' @author Georges Monette (georges@@yorku.ca)

#' \code{\link{cell}}, \code{\link{dell}},
#' \code{\link{ellplus}},\code{\link{dellplus}},

#' @keywords ellipse ellipse geometry
#' @examples
#'
#' v1 <- 36*diag(2)
#' v2 <- .5 * (diag(2) - .4)
#' v2[2,2] <- 2
#' plot( 0:10,0:10,type = 'n')
#' lines( ell( c(2,2), v1))
#' lines( ell( c(4,4), v2), col = 'red')
#' osculant(  c(2,2), v1, c(4,4), v2, n = 3)
#' osculant(  c(2,2), v1, c(4,4), v2, n = 1)
#' lines( osculant( c(2,2), v1, c(4,4), v2), col = 'red')
#'
#' lines( ell( c(8,8), v2), col = 'blue')
#' lines( osculant( c(2,2), v1, c(8,8), v2), col = 'blue')
#' points( osculant( c(2,2), v1, c(8,8), v2, n=1),pch = 16, col = 'blue')
#' points( osculant( c(2,2), v1, c(8,8), v2, n=0),pch = 16, col = 'blue')
#' points( osculant( c(8,8), v2, c(2,2), v1,  n=0),pch = 16, col = 'blue')
#'
#'
#' @export
osculant <- function(x, ...) UseMethod("osculant")

osculant.default <- function( center1, shape1, center2, shape2, n = 100, range =c(0,1), maxu = 100) {
  # Use solution from Lagrangean:
  # p = ( shape1^-1 + lam * shape1^-1)^-1 %*% lam2 shape2^-1 delta
  pt <- function(u)  u* solve( u*diag(p) + (1-u)* shape, delta)


  #' p - quick paste with sep = ''
  #'
  #' Works like \code{paste}, using an empty separator string.
  #'
  #'
  #' @param \dots one or more R objects, to be converted to character vectors.
  #' @return A character vector of the concatenated values.
  #' @author Georges Monette
  #' @seealso \code{\link[base]{paste}}
  #' @keywords manip
  #' @examples
  #'
  #' p(letters[1:5], 1:5)
  #'
  p <- nrow(shape1)
  delta <- center2 - center1
  shape <- t(solve(shape1,shape2))
  shape <- shape/mean(diag(shape))   # attempt to equalize intervals over range
  if( n == 0) {
    norm1 <- function( u ) {
      pp <- pt(u)
      sqrt( sum( pp  * solve( shape1, pp))) -1
    }
    if ( norm1(1) < 0 ) {
      warning( "Center of second ellipse inside first ellipse")
      return( NULL)
    }
    u <- uniroot( norm1, c(0,1))$root
    rbind( pt(u) + center1)

  } else {
    vec <- sapply( seq(range[1],range[2],diff(range)/n), pt)
    t( vec + center1 )
  }
}

#
#
#
# v1 <- 36*diag(2)
# v2 <- .5 * (diag(2) - .4)
# v2[2,2] <- 2
# plot( 0:10,0:10,type = 'n')
# lines( ell( c(2,2), v1))
# lines( ell( c(4,4), v2), col = 'red')
# osculant(  c(2,2), v1, c(4,4), v2, n = 3)
# osculant(  c(2,2), v1, c(4,4), v2, n = 1)
# lines( osculant( c(2,2), v1, c(4,4), v2), col = 'red')
#
# lines( ell( c(8,8), v2), col = 'blue')
# lines( osculant( c(2,2), v1, c(8,8), v2), col = 'blue')
# points( osculant( c(2,2), v1, c(8,8), v2, n=1),pch = 16, col = 'blue')
# points( osculant( c(2,2), v1, c(8,8), v2, n=0),pch = 16, col = 'blue')
# points( osculant( c(8,8), v2, c(2,2), v1,  n=0),pch = 16, col = 'blue')
#











# Replaced with version below from p3d
# ell <- function(center = rep(0,2) , shape = diag(2) , radius = 1, n = 100,
#         angles = (0:n)*2*pi/n) {
#        circle <- radius * cbind( cos(angles), sin(angles))
#        t( c(center) + t( circle %*% fac(shape)))
# }




#' Old confidence ellipse
#'
#'
#'
#'
#'
#' @param model
#' @param which.coef
#' @param levels
#' @param Scheffe
#' @param dfn
#' @param center.pch
#' @param center.cex
#' @param segments
#' @param xlab
#' @param ylab
#' @param las
#' @param col
#' @param lwd
#' @param lty
#' @param add
#' @param \dots

#' @author Georges Monette


#' @keywords ~kwd1 ~kwd2
#' @examples
#'
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#'
#' ## The function is currently defined as
#' function (model, which.coef, levels = 0.95, Scheffe = FALSE,
#'     dfn = 2, center.pch = 19, center.cex = 1.5, segments = 51,
#'     xlab, ylab, las = par("las"), col = palette()[2], lwd = 2,
#'     lty = 1, add = FALSE, ...)
#' {
#'     help <- "\nSee help for car::confidence.ellipse.lm\nexcept that 'cell' returns the points to form the ellipse\nwhich must be plotted with plot(...,type='l') or lines(...)\n-- Use dfn to determine Sheffe dimension, i.e. dfn = 1 to generate ordinary CIs, dfn = 2 for 2-dim CE, etc.\n"
#'     require(car)
#'     which.coef <- if (length(coefficients(model)) == 2)
#'         c(1, 2)
#'     else {
#'         if (missing(which.coef)) {
#'             if (has.intercept(model))
#'                 c(2, 3)
#'             else c(1, 2)
#'         }
#'         else which.coef
#'     }
#'     coef <- coefficients(model)[which.coef]
#'     xlab <- if (missing(xlab))
#'         paste(names(coef)[1], "coefficient")
#'     ylab <- if (missing(ylab))
#'         paste(names(coef)[2], "coefficient")
#'     if (missing(dfn)) {
#'         if (Scheffe)
#'             dfn <- sum(df.terms(model))
#'         else 2
#'     }
#'     dfd <- df.residual(model)
#'     shape <- vcov(model)[which.coef, which.coef]
#'     ret <- numeric(0)
#'     for (level in rev(sort(levels))) {
#'         radius <- sqrt(dfn * qf(level, dfn, dfd))
#'         ret <- rbind(ret, c(NA, NA), ell(coef, shape, radius))
#'     }
#'     colnames(ret) <- c(xlab, ylab)
#'     ret
#'   }
#'
old.cell <-
  function (model, which.coef, levels = 0.95, Scheffe = FALSE, dfn = 2,
            center.pch = 19, center.cex = 1.5, segments = 51, xlab, ylab,
            las = par("las"), col = palette()[2], lwd = 2, lty = 1,
            add = FALSE, ...)
  {
    help <- "
    See help for car::confidence.ellipse.lm
    except that 'cell' returns the points to form the ellipse
    which must be plotted with plot(...,type='l') or lines(...)
    -- Use dfn to determine Sheffe dimension, i.e. dfn = 1 to generate ordinary CIs, dfn = 2 for 2-dim CE, etc.
    "
    require(car)
    which.coef <- if (length(coefficients(model)) == 2)
      c(1, 2)
    else {
      if (missing(which.coef)) {
        if (has.intercept(model))
          c(2, 3)
        else c(1, 2)
      }
      else which.coef
    }
    coef <- coefficients(model)[which.coef]
    xlab <- if (missing(xlab))
      paste(names(coef)[1], "coefficient")
    ylab <- if (missing(ylab))
      paste(names(coef)[2], "coefficient")
    if(missing(dfn)) {
      if (Scheffe) dfn <- sum(df.terms(model))
      else 2
    }
    dfd <- df.residual(model)
    shape <- vcov(model)[which.coef, which.coef]
    ret <- numeric(0)
    for (level in rev(sort(levels))) {
      radius <- sqrt(dfn * qf(level, dfn, dfd))
      ret <- rbind(ret, c(NA,NA), ell( coef, shape, radius) )
    }
    colnames(ret) <- c(xlab, ylab)
    ret
  }

# from Plot3d.R





#' Calculate coordinates of a data ellipse
#'
#'
#' \code{dell} to calculates the coordinates of a 2D data ellipse
#' (concentration ellipse) from (X, Y) variables.
#'
#' \code{dellplus} can produce, in addition to the points of an ellipse, the
#' conjugate axes corresponding to a \code{chol} or other decomposition and the
#' surrounding parallelogram defined by these axes.
#'
#' These functions simply calculate the mean vector and covariance matrix and
#' call \code{ell} or \code{ellplus}.
#'
#' @aliases dell dellplus
#' @param x,y Either a two-column matrix or numeric vectors of the same length
#' @param radius Radius of the ellipse-generating unit circle.  The default,
#' \code{radius=1} corresponds to a "standard" ellipse.
#' @param \dots Other arguments passed down to \code{ell} or \code{ellplus}.
#' @return Returns a 2-column matrix of (X,Y) coordinates suitable for drawing
#' with \code{lines()}.
#'
#' For \code{dellplus}, when more than one of the options \code{ellipse},
#' \code{diameters}, and \code{box} is \code{TRUE}, the different parts are
#' separated by a row of \code{NA}.
#' @author Georges Monette
#' @seealso \code{\link{cell}}, \code{\link{ell}}, \code{\link{ellplus}},
#' @references Monette, G. (1990). Geometry of Multiple Regression and
#' Interactive 3-D Graphics. In Fox, J. & Long, S. (ed.)  \emph{Modern Methods
#' of Data Analysis}, Sage Publications, 209-256.
#' @keywords dplot aplot
#' @examples
#'
#' data(Prestige)   # from car
#' attach(Prestige)
#' fit.simple <- lm( prestige ~ education, Prestige)
#'
#' plot(prestige ~ education, type='p')
#' lines(dell(education, prestige), col="blue", lwd=3)
#' lines(bbox <- dellplus(education, prestige, box=TRUE))
#' lines(dellplus(education, prestige, diameter=TRUE, radius=2), col="gray")
#' detach(Prestige)
#'
#'
#' @export
dell <- function( x, y, radius = 1, ...) {
  if ( (is.matrix(x) && (ncol(x) > 1))|| is.data.frame(x)) mat <- as.matrix(x[,1:2])
  else if (is.list(x)) mat <- cbind(x$x, x$y)
  else mat <- cbind( x,y)
  ell( apply(mat,2,mean), var(mat), radius = radius, ...)

}














"</pre>

== Diagnostics ==
<pre>"




#' Generic diagnostics
#'
#' Generic diagnostics
#'
#'
#'
#' @param x
#' @param \dots

#' @author Georges Monette


#' @keywords ~kwd1 ~kwd2
#' @examples
#'
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#'
#' ## The function is currently defined as
#' function (x, ...)
#' UseMethod("diags")
#'
#' @export
diags <- function(x, ...) UseMethod("diags")




#' Standard diagnostics for lm objects
#'
#' Standard diagnostics for lm objects
#'
#'
#'
#' @param x
#' @param y
#' @param \dots
#' @param ask
#' @param labels
#' @param showlabs

#' @author Georges Monette


#' @keywords ~kwd1 ~kwd2
#' @examples
#'
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#'
#' ## The function is currently defined as
#' function (x, y, ..., ask, labels = names(residuals(x)), showlabs = text)
#' {
#'     if (!missing(ask)) {
#'         op <- par(ask = ask)
#'         on.exit(par(op))
#'     }
#'     form <- formula(x)
#'     f <- predict(x)
#'     r <- residuals(x)
#'     nams <- names(r)
#'     if (!missing(labels)) {
#'         nams <- names(residuals(x))
#'         if (length(nams) != length(labels))
#'             labels <- labels[nams]
#'     }
#'     ret <- NULL
#'     if (missing(y)) {
#'         y <- f + r
#'         yname <- deparse(form[[2]])
#'     }
#'     else yname <- deparse(substitute(y))
#'     fname <- paste("Fitted:", deparse(form[[3]]), collapse = " ")
#'     plot(f, y, xlab = fname, ylab = yname, main = "Dependent var. vs. Predicted",
#'         ...)
#'     abline(0, 1, lty = 1)
#'     lines(supsmu(f, y))
#'     showlabs(f, y, labels, ...)
#'     lmi <- lm.influence(x)
#'     hat <- lmi$hat
#'     sigma <- lmi$sigma
#'     mm <- scale(model.matrix(x), scale = F)
#'     mp <- predict(x, type = "terms")
#'     comp.res <- mp + r
#'     plot(f, abs(r), xlab = fname, ylab = deparse(substitute(abs(resid(x)))),
#'         main = "Absolute Residual vs. Predicted", ...)
#'     showlabs(f, abs(r), labels, ...)
#'     zq <- qqnorm(r, main = "Normal Quantile Plot", ylab = "Residual",
#'         sub = fname)
#'     qqline(r)
#'     showlabs(zq, labels, ...)
#'     n <- length(r)
#'     r.o <- sort(r)
#'     half <- (n + 1)/2
#'     if (n%%2 == 1) {
#'         med <- r.o[half]
#'         below <- med - r.o[half:1]
#'         above <- r.o[half:n] - med
#'     }
#'     else {
#'         med <- sum(r.o[c(half, half + 1)])/2
#'         below <- med - r.o[(n/2):1]
#'         above <- r.o[(n/2 + 1):n] - med
#'     }
#'     opt <- par(pty = "s")
#'     ran <- range(c(below, above))
#'     plot(below, above, main = "Symmetry plot of residuals", xlab = "Distance below median",
#'         ylab = "Distance above median", xlim = ran, ylim = ran)
#'     abline(0, 1, lty = 2)
#'     par(opt)
#'     std.r <- r/(sigma * sqrt(1 - hat))
#'     plot(hat, std.r, xlab = "Leverage (hat)", ylab = yname, sub = fname,
#'         main = "Studentized residual vs. Leverage", ...)
#'     showlabs(hat, std.r, labels, ...)
#'     nams <- dimnames(lmi$coefficients)[[1]]
#'     pairs(lmi$coefficients)
#'     pairs(lmi$coefficients, panel = function(x, y, nams) {
#'         points(x, y)
#'         text(x, y, nams)
#'     }, nams = nams)
#'     invisible(0)
#'   }
#'
#' @export
diags.lm <- function(x, y, ..., ask, labels = names(residuals(x)), showlabs = text)
{
  # diags.lm
  # graphical diagnostics for lm, locally first-order for glm
  # enlarged version of plot.lm with emphasis on diagnostics
  # G. Monette, Dec. 94
  # modified Nov. 97, May 98
  # Slight modification to pairs adding labels, Jan 03
  if(!missing(ask)) {
    op <- par(ask = ask)
    on.exit(par(op))
  }
  form <- formula(x)
  f <- predict(x)
  r <- residuals(x)
  nams <- names(r)
  if(!missing(labels)) {
    nams <- names(residuals(x))	#
    # if labels not same length as residuals assume it's a vector
    # of len == original data and select elements included in residuals
    if(length(nams) != length(labels))
      labels <- labels[nams]
  }
  ret <- NULL
  if(missing(y)) {
    y <- f + r
    yname <- deparse(form[[2]])
  }
  else yname <- deparse(substitute(y))
  fname <- paste("Fitted:", deparse(form[[3]]), collapse = " ")
  plot(f, y, xlab = fname, ylab = yname, main = "Dependent var. vs. Predicted",
       ...)
  abline(0, 1, lty = 1)
  lines(supsmu(f,y))
  showlabs(f, y, labels,...)
  #
  # get influence diags and model matrix while looking at first plot
  #
  lmi <- lm.influence(x)
  hat <- lmi$hat
  sigma <- lmi$sigma	# drop 1 sigma
  mm <- scale(model.matrix(x), scale = F)	# centres each column
  mp <- predict(x, type = "terms")
  comp.res <- mp + r	# effect + residual
  #
  # Absolute residual vs. predicted
  #
  plot(f, abs(r), xlab = fname, ylab = deparse(substitute(abs(resid(x)))),
       main = "Absolute Residual vs. Predicted", ...)
  showlabs(f, abs(r), labels, ...)	#
  #
  # Normal quantile plot
  #
  zq <- qqnorm(r, main = "Normal Quantile Plot", ylab = "Residual", sub
               = fname)
  qqline(r)
  showlabs(zq, labels,...)	#
  #
  # Symmetry plot of residuals (Lawrence C. Hamilton, Regression with
  #       Graphics, Duxbury, 1992)
  n <- length(r)
  r.o <- sort(r)
  half <- (n + 1)/2
  if (n%%2 == 1) {    # n is odd
    med <- r.o[half]
    below <- med - r.o[half:1]
    above <- r.o[half:n] - med
  }
  else {
    # n is even
    med <- sum(r.o[c(half, half + 1)])/2
    below <- med - r.o[(n/2):1]
    above <- r.o[(n/2 + 1):n] - med
  }
  opt <- par(pty = "s")
  ran <- range(c(below, above))
  plot(below, above, main = "Symmetry plot of residuals", xlab =
         "Distance below median", ylab = "Distance above median", xlim
       = ran, ylim = ran)
  abline(0, 1, lty = 2)
  par(opt)	#

  #
  # Studentized residual vs. leverage
  #

  std.r <- r/(sigma * sqrt(1 - hat))
  plot(hat, std.r, xlab = "Leverage (hat)", ylab = yname, sub = fname,
       main = "Studentized residual vs. Leverage", ...)
  showlabs(hat, std.r, labels,...)	#	plot(lmi$sig, std.r)	#

  #
  # effect of dropping one observation DFBETA
  #

  nams <- dimnames(lmi$coefficients)[[1]]
  pairs(lmi$coefficients)
  pairs(lmi$coefficients, panel = function(x,y,nams){
    points(x,y)
    text(x,y,nams)
  }, nams = nams)

  # main = "Effect of dropping one case", sub = fname)
  invisible(0)
}

#' @export
model.matrix.lme <- function( fit , data = fit$data,
                              na.action = fit$na.action,
                              ...){
  mCall <- fit$call
  fixed <- eval(eval(mCall$fixed)[-2])
  data <- model.frame(fixed, data = data)
  naresid( na.action, model.matrix(fixed, data = data))
}
# model.matrix(fit, data=pred) %>% dim
# model.frame(fit, data=pred) %>% dim
# model.matrix(fit, na.action = na.pass) %>% dim
# model.frame(fit, na.action = na.pass) %>% dim


# BUG: not working as it should for na.action=na.exclude
# model.frame.lme <- function( fit , data = fit$data,
#                              na.action = fit$call$na.action,...)
#   model.frame(formula(fit), data = data, na.action = na.action)

#' @export
model.frame.lme <- function (object, data =object$data, na.action = object$na.action,
                             ...)
{
  # adapted from portions of predict.lme
  mCall <- object$call
  fixed <- eval(eval(mCall$fixed)[-2])
  Terms <- object$terms
  data <- as.data.frame(data)
  mfArgs <- list(formula = fixed,
                 data = data, na.action = na.action, drop.unused.levels = TRUE)
  dataMix <- do.call("model.frame", mfArgs)
  dataMix
}





#' Standard diagnostics for lme objects
#'
#'
#'
#'
#'
#' @param \dots

#' @author Georges Monette


#' @keywords ~kwd1 ~kwd2
#' @examples
#'
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#'
#' ## The function is currently defined as
#' function (...)
#' cat("Being implemented")
#'
#' @export
diags.lme <- function( ... ) cat("Being implemented")



"
</pre>
== vplot -- a plot function for matrix algebra ==
<pre>
"

##
##  vplot  plots columns of a 2 x n matrix
##  Transferred to coursfun: Nov. 15, 2005



#' Collection of functions to help teach matrix geometry in 2 dimensions
#'
#'
#'
#' vplot - plots the columns of a 2 x n matrix or a vector of length 2 - vplot
#' adds to the current plot resizing it to include all plotted objects in a
#' 'euclidean' frame - to start a new plot, use 'new = T' - to remove the last
#' element added use 'vplot(pop=1)' Associated functions: - vell( mean, var)
#' generates an ellipse, default = unit circle - vbox() generates a box -
#' vobj() generates a circle in a box - orthog(theta) generates an orthog
#' matrix rotating through angle theta - orthog.proj generates the matrix of an
#' orthog. projection into span (x) - vmat( .... ) generates a 2 by n matrix
#' Examples: vplot( new = T ) vplot( vell(), 'l' ) vplot( cbind(c(3,1),c(1,4))
#' \%*\% vell()) vplot( pop = 1) vplot( cbind(c(3,1),c(1,4)) \%*\% vell(), type
#' = 'l', col = 'red')
#' above ~~
#'
#' @param mat
#' @param type
#' @param new
#' @param pch
#' @param pop
#' @param \dots

#' @author Georges Monette


#' @keywords ~kwd1 ~kwd2
#' @examples
#'
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#'
#' ## The function is currently defined as
#' function (mat, type = "p", new = F, pch = 16, pop = 0, ...)
#' {
#'
#'
#'     if (new || !exists(".vplot"))
#'         assign(".vplot", list(list(x = 0, y = 0, type = "n")),
#'             pos = 1)
#'     a <- .vplot
#'     if (!missing(mat)) {
#'         mat <- cbind(mat)
#'         if (type == "v") {
#'             zz <- rbind(0 * mat, mat, mat/0)
#'             mat <- matrix(c(zz), nrow = 2)
#'             type = "b"
#'         }
#'         d <- dim(mat)
#'         if (d[1] != 2 && d[2] == 2) {
#'             mat <- t(mat)
#'             warning("mat is n x 2 and has been transposed")
#'         }
#'         a <- c(a, list(list(x = mat[1, ], y = mat[2, ], type = type,
#'             pch = pch, ...)))
#'     }
#'     dat <- NULL
#'     for (i in seq(along = a)) {
#'         dat <- c(dat, a[[i]]$x, a[[i]]$y)
#'     }
#'     par(pty = "s")
#'     plot(range(na.omit(dat)), range(na.omit(dat)), type = "n",
#'         xlab = "", ylab = "")
#'     if (pop > 0) {
#'         keep <- 1:max(1, (length(a) - (pop + 1)))
#'         a <- a[keep]
#'     }
#'     abline(h = 0, v = 0)
#'     for (i in seq(along = a)) do.call("points", a[[i]])
#'     assign(".vplot", a, pos = 1)
#'     invisible(a)
#'   }
#'
#' @export
vplot <- function( mat , type = 'p', new = F,  pch = 16, pop = 0, ...) {
  help <- "
  vplot    - plots the columns of a 2 x n matrix or a vector of length 2
  - vplot adds to the current plot resizing it to include all plotted
  objects in a 'euclidean' frame
  - to start a new plot, use 'new = TRUE'
  - to remove the last element added use 'vplot(pop=1)'
  Associated functions:
  - vell( mean, var) generates an ellipse, default = unit circle
  - vbox() generates a box
  - vobj() generates a circle in a box
  - orthog(theta) generates an orthog matrix rotating through angle theta
  - orthog.proj generates the matrix of an orthog. projection into span (x)
  - vmat( .... ) generates a 2 by n matrix
  Examples:
  vplot( new = TRUE )
  vplot( vell(), 'l' )
  vplot( cbind(c(3,1),c(1,4)) %*% vell())
  vplot( pop = 1)
  vplot( cbind(c(3,1),c(1,4)) %*% vell(), type = 'l', col = 'red')
  "
  if (  new || !exists(".vplot")) assign(".vplot", list(list(x=0,y=0,type='n')),pos=1)
  a <- .vplot
  if ( ! missing(mat) ) {
    mat <- cbind(mat)
    if ( type == 'v' ) {
      zz <- rbind( 0*mat, mat, mat/0)
      mat <- matrix( c(zz), nrow = 2)
      type = 'b'
    }
    d <- dim(mat)
    if ( d[1] != 2 && d[2] == 2){
      mat <- t(mat)
      warning("mat is n x 2 and has been transposed")
    }
    a <- c(a,list( list(x=mat[1,],y = mat[2,],type=type, pch = pch, ...)))
  }
  dat <- NULL
  for ( i in seq( along = a )) {
    dat <- c( dat, a[[i]]$x, a[[i]]$y)
  }
  # print(a)
  par ( pty = 's')
  plot( range(na.omit(dat)), range(na.omit(dat)), type = 'n', xlab = '', ylab ='')
  if ( pop > 0 ) {
    keep <- 1:max(1,(length(a)-(pop+1)))
    a <- a[keep]
  }
  abline( h = 0, v = 0)
  for ( i in seq( along = a)) do.call('points', a[[i]])
  assign(".vplot", a, pos = 1)
  invisible(a)
}



#' Vector around an ellipse
#'
#'
#'
#'
#'
#' @param \dots

#' @author Georges Monette


#' @keywords ~kwd1 ~kwd2
#' @examples
#'
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#'
#' ## The function is currently defined as
#' function (...)
#' t(ell(...))
#'
#' @export
vell <- function(...) t( ell(...))


#' Unit box
#'
#'
#'
#'
#'
#' @param \dots

#' @author Georges Monette


#' @keywords ~kwd1 ~kwd2
#' @examples
#'
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#'
#' ## The function is currently defined as
#' function (...)
#' cbind(c(-1, -1), c(-1, 1), c(1, 1), c(1, -1), c(-1, -1))
#'
#' @export
vbox <- function(...) cbind( c(-1,-1), c(-1,1), c(1,1), c(1,-1), c(-1,-1))


#' Combine ellipse with subtending box
#'
#'
#'
#'
#'
#' @param \dots

#' @author Georges Monette


#' @keywords ~kwd1 ~kwd2
#' @examples
#'
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#'
#' ## The function is currently defined as
#' function (...)
#' {
#'     cbind(vell(), NA, vbox(), NA, c(0, -1), c(0, 1), NA, c(-1,
#'         0), c(1, 0))
#'   }
#'
#' @export
vobj <- function(...) {
  cbind( vell(), NA, vbox(), NA, c(0,-1),c(0,1), NA, c(-1,0), c(1,0))
}


#' Square in 2 dimensions
#'
#'
#'
#'
#'
#' @param \dots

#' @author Georges Monette


#' @keywords ~kwd1 ~kwd2
#' @examples
#'
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#'
#' ## The function is currently defined as
#' function (...)
#' vmat(0, 0, 0, 1, 1, 1, 1, 0, 0, 0)
#'
#' @export
vsquare <- function(...) vmat( 0,0,0,1,1,1,1,0,0,0)



#' Create a matrix entering vectors column by column
#'
#'
#'
#'
#'
#' @param \dots

#' @author Georges Monette


#' @keywords ~kwd1 ~kwd2
#' @examples
#'
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#'
#' ## The function is currently defined as
#' function (...)
#' {
#'     help <- "\nvmat creates a matrix entering data column by column\n"
#'     aa <- list(...)
#'     aa <- do.call("c", aa)
#'     matrix(aa, nrow = 2)
#'   }
#'
#' @export
vmat <- function(...) {
  help <- "
  vmat creates a matrix entering data column by column
  "
  aa <- list(...)
  aa <- do.call('c', aa)
  matrix(aa, nrow = 2)
}

