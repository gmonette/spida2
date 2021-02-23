## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, comment='', message = FALSE)
library(spida2)
library(carEx)
library(magrittr)
basis  <- function(X, tol = 1e-9) {
  # returns linear independent columns
  # with possible pivoting
  q <- qr(X, tol = tol)
  sel <- q$pivot[seq_len(q$rank)]
  ret <- X[, sel, drop = FALSE]
  colnames(ret) <- colnames(X)[sel]
  ret
}

Xmat <-	function(x,
                 degree,
                 D = 0,
                 signif = 3) {
  # Returns rows of design matrix if D = 0
  # or linear hypothesis matrix for D-th derivative
  if (length(D) < length(x))
    D = rep(D, length.out = length(x))
  if (length(x) < length(D))
    x = rep(x, length.out = length(D))
  xmat = matrix(x, nrow = length(x), ncol = degree + 1)
  expvec <- 0:degree
  coeffvec <- rep(1, degree + 1)
  expmat <- NULL
  coeffmat <- NULL
  
  for (i in 0:max(D)) {
    expmat <- rbind(expmat, expvec)
    coeffmat <- rbind(coeffmat, coeffvec)
    coeffvec <- coeffvec * expvec
    expvec <- ifelse(expvec > 0, expvec - 1, 0)
  }
  X = coeffmat[D + 1,,drop = FALSE] * xmat ^ expmat[D + 1,, drop = FALSE]
  
  xlab = signif(x, signif)
  rownames(X) = ifelse(D == 0,
                       paste('f(', xlab , ')', sep = ''),
                       paste("D", D, "(", xlab, ")", sep = ""))
  colnames(X) = paste("X", 0:(ncol(X) - 1), sep = "")
  class(X) <- c('gspline_matrix',class(X)) 
  X
}
Xf <-
  function(x,
           knots,
           degree = 3,
           D = 0,
           right = TRUE,
           periodic = FALSE,
           signif = 3) {
    # Returns block diagonal matrix (if x ordered)
    # with design matrix or linear hypothesis matrix
    # for the 'full' model with contrained parameters.
    #
    # With the default, right == TRUE, 
    # if x is at a knot then it is included in
    # in the lower interval. 
    # With right = FALSE, it is included in the higher
    # interval. This is needed when building 
    # derivative constraints at the knot
    if(periodic) {
      period <- max(knots)
      xx <- x %% period
      if(right) xx[xx==0] <- period
      x <- xx
      knots <- knots[-length(knots)]
    }
    xmat = Xmat (x, degree, D , signif)
    k = sort(knots)
    g = cut(x, c(-Inf, k, Inf), right = right)
    ret <- do.call('cbind',
                   lapply(seq_along(levels(g)), function(i)
                     (g == levels(g)[i]) *  xmat))
    if(periodic) rownames(ret) <- 
      sub(')',paste0(' mod ',period,')'), rownames(ret))
  class(ret) <- c('gspline_matrix',class(ret)) 

        ret
  }

# 
# Value and derivatives at 0
#
# And all discontinuities at knots
#
Dmat <- function(knots, degree, periodic = FALSE, signif = 3) {
  dm <- max(degree)
  cmat <- Xf(0, knots, dm, D=0:dm, periodic = periodic)
  n_knots <- length(knots)
  for (i in seq_len(n_knots - 1) ) {
    k <- knots[i]
    dmat <- Xf(k, knots, dm, D = seq(0,dm), F, periodic = periodic) -   
      Xf(k, knots, dm, D = seq(0,dm), T, periodic = periodic)
    rownames( dmat ) <- paste( "C(",signif(k, signif),").",
                               seq(0,dm), sep = '')
    cmat <- rbind( cmat,  dmat)
  }
  k <- knots[length(knots)]
  if(periodic) {
    dmat <- Xf(0, knots, dm, D = seq(0,dm) , F ,periodic = periodic) -   
      Xf(k, knots, dm, D = seq(0,dm) , T ,periodic = periodic)
    rownames( dmat ) <- paste( "C(0 mod ",signif(k, signif),").",
                               seq(0,dm), sep = '')
    cmat <- rbind(cmat, dmat)
  } else {
    dmat <- Xf(k, knots, dm, D = seq(0,dm) , F ,periodic = periodic) -
      Xf(k, knots, dm, D = seq(0,dm) , T ,periodic = periodic)
    rownames( dmat ) <- paste( "C(",signif(k, signif),").",
                               seq(0,dm), sep = '')
    cmat <- rbind(cmat, dmat)
  }
    class(cmat) <- c('gspline_matrix',class(cmat)) 

  cmat
}

# 
# Parameters to constrain
#
# TODO: Change so degree can be a list
#
Pcon <- function(knots, degree, periodic) {
  degree <- rep( degree, length.out = length(knots) + 1)
  if(periodic) {
    if(degree[length(degree)] != degree[1]) warning("For periodic splines, the degree of the last and first intervals should match")
    knots <- knots[-length(knots)]
    degree <- degree[-length(degree)]
  }
  dm <- max(degree)
  cmat <- NULL
  for ( i in seq_along(degree)) {
    di <- degree[i]
    if ( dm > di ) {
      dmat <- diag( (length(knots) + 1) * (dm +1)) [  
        (i - 1)*(dm + 1) + 1 + seq( di+1,dm), , drop = F]
      rownames( dmat ) = paste( "I.", i,".",seq(di+1,dm),sep = '')
      cmat = rbind( cmat, dmat)
    }
  }
    class(cmat) <- c('gspline_matrix',class(cmat)) 

  cmat
}

## ------------------------------------------------------------------------
Xf(0:9, knots = c(3,7), degree = 3)

## ------------------------------------------------------------------------
Cmat(knots = c(3, 7), degree = c(2, 3, 2), smooth = c(1, 2))

## ------------------------------------------------------------------------
Emat(knots = c(3, 7), degree = c(2, 3, 2), smooth = c(1, 2))

## ------------------------------------------------------------------------
sp <- gspline(knots = c(3, 7), degree = c(2, 3, 2), smoothness = c(1, 2))
sp(0:9)

## ------------------------------------------------------------------------
df <- data.frame(x = 0:10)
set.seed(123)
df <- within(df, y <- -2* (x-5) + .1 * (x-5)^3 + rnorm(x))
df <- rbind(df, data.frame(x = seq(0,10,.1), y = NA))
df <- sortdf(df, ~ x)
plot(y~x, df, pch = 16)
fit <- lm(y ~ sp(x), data = df)
summary(fit)
lines(df$x , predict(fit, df))

## ------------------------------------------------------------------------
sp <- gspline(knots = c(3,7), degree = c(2,3,2), smoothness = c(1,2))
sp(0:9)

## ------------------------------------------------------------------------
sp(c(2, 3, 7), D = 1)

## ------------------------------------------------------------------------
sp(c(3, 3, 3), D = 2, limit = c(-1,0,1))

## ------------------------------------------------------------------------
# xpred <- seq(0,10, .05)
# A.1 <- cbind(0, sp(xpred, D = 1))
# A.2 <- cbind(0, sp(xpred, D = 2))
# ww.1 <- as.data.frame(wald(fit, A.1))
# ww.2 <- as.data.frame(wald(fit, A.2))
# 
# plot(xpred, ww.1$coef, type = 'l', lwd = 2)
# plot(xpred, ww.2$coef, type = 'l', lwd = 2)
# library(latticeExtra)
# ww.1$x <- xpred
# xyplot(coef ~ x, ww.1, type = 'l',
# 	   lower = ww.1$L2, upper = ww.1$U2,
# 	   subscripts = TRUE) +
# 	layer(gpanel.fit(...))
# head(ww.1)

## ------------------------------------------------------------------------
Cmat(knots = 5, degree = c(3,3), smooth = 3)

## ------------------------------------------------------------------------
C <- Cmat(5, c(3,3), 3)[c(1,3,4),]
C

## ------------------------------------------------------------------------
sp2 <- function(x) gsp(x, knots = 5, smooth = -1, degree = c(3,3), lin = C)
sp2(seq(0,2,.5))

## ------------------------------------------------------------------------
fit <- lm(y ~ sp2(x), df)
summary(fit)
df$fit <- predict(fit, df)
df$fit[df$x == 5] <- NA
plot(fit ~ x, df, type = 'l')

