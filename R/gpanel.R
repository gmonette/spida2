##
## panel functions:
##
## gpanel.fit incorporating gpanel.band
## gpanel.text
##
## designed to work easily with 'layer' and 'glayer' in latticeExtra
## Added: 2014 09 18
## Author: GM
##
#' Panel functions for predicted values and SE bands
#'
#' Panel functions for predicted values and SE bands using 'layer' and 'glayer' in the package latticeExtra
#'
#' With 'layer' and 'glayer' in 'latticeExtra', these functions can be used to easily generate fitted values and
#' confidence or prediction bands that have a reasonable appearance whether a plot uses 'groups' or not.
#'
#
#' @param data data frame to be used to add additional values of numeric variable
#' @param form formula evaluated in data. The first term defines the variable with values to be filled in and the remaining terms define the variables to be used for grouping determining the minima and maxima within which values are added.
#' @param n the number of values to be added between the global mininum and maximum.
#'        Values falling outside conditional minima and maxima are culled. Default 200.
#' @param xpd expansion factor to add points beyond minima and maxima. Default 1.0.
#' @param y
#' @param fit fitted values of a model, generally passed through 'layer' from a call to 'xyplot': e.g. \code{xyplot( y ~ x, data, groups = g, fit = data$yhat, lower = with(data, yhat - 2*se), upper = with(data, yhat + 2*se), subscripts = T)}
#' @param lower,upper
#' @param subscripts
#' @param dots
#' @param col
#' @param group.number
#' @param alpha
#' @param col.symbol is used to control color when using 'groups'
#' @param border default = FALSE for panel.band.
#' @param font
#' @param fontface
#' @return The 'panel.bands', 'panel.fit', and 'panel.labels' functions are
#'         invoked for their graphical effect.
#' @author Georges Monette <georges@@yorku.ca>
#' @examples
#' \dontrun{
#'   library(yscs)
#'   library(latticeExtra)
#'   library(car)
#'
#'   fit <- lm(prestige ~ (income + I(income^2)) * type, Prestige,
#'       na.action = na.exclude)
#'   pred <- cbind(Prestige, predict(fit, newdata = Prestige, se = TRUE))
#'   head(pred)
#'   (p <- xyplot( prestige ~ income , pred,  groups = type,
#'                 subscripts = T,
#'                 fit = pred$fit,
#'                 lower = with(pred, fit - 2*se.fit),
#'                 upper = with(pred, fit + 2*se.fit)))
#'   p + glayer(gpanel.fit(...))
#'
#'   ###
#'   ### Use 'fillin' to add points in sparse regions of the predictor
#'   ### to produce smoother bands
#'   ###
#'
#'   fit <- lm( income ~
#'         (education+I(education^2)+I(education^3) +I(education^4))* type,
#'         Prestige, na.action = na.exclude)    # overfitting!
#'   # adding extra values of predictor to get smooth line
#'   Prestige$occupation <- rownames(Prestige)
#'   z <- fillin(Prestige, ~education + type, xpd = 1.1)
#'
#'   dim(z)
#'   dim(Prestige)
#'   z <- cbind(z, predict(fit, newdata = z, se = TRUE))
#'   head(z)
#'   gd(3,cex=2,lwd=2, alpha = .7)
#'   (p <-  xyplot( income ~ education, z, groups = type,
#'                  subscripts = T,
#'                  fit = z$fit,
#'                  lower = z$fit - z$se,
#'                  upper = z$fit + z$se,
#'                  auto.key = list(space='right', lines= T)))
#'   p + glayer( gpanel.fit(...))
#'   p + glayer( gpanel.fit(...,  alpha = .1))
#'   # Using spida:gd() to get a ggplot-like appearance
#'   gd(3,lty=1,lwd=2)
#'   p + glayer( gpanel.fit(...,alpha = .1))
#'
#'   ###
#'   ###  Using gpanel.fit with no groups
#'   ###
#'
#'   (p <-  xyplot( income ~ education| type, z,
#'                  subscripts = T,
#'                  fit = z$fit,
#'                  lower = z$fit - z$se,
#'                  upper = z$fit + z$se,
#'                  auto.key = list(space='right', lines= T)))
#'   p + layer( panel.fit(...))
#'
#'   # gd_(basecol = 'tomato4')  # Use 'gd_' to set parameters without groups
#'   gd_(base = 'tomato4')  # Use 'gd_' to set parameters without groups
#'   p + layer( panel.fit(...))
#'   p + layer( panel.fit(...,  col = 'grey10'))
#'
#'   ###
#'   ### With panels and groups
#'   ###
#'
#'   z <- Prestige
#'   z$gender <- with(z, cut( women, c(-1,15,50,101),labels = c("Male","Mixed","Female")))
#'   tab(z, ~ gender + type)
#'   z <- fillin( z, ~ education + type + gender, xpd = 1.1)
#'   fit <- lm( income ~ (education+I(education^2)+I(education^3) )* type * gender,
#'              z, na.action = na.exclude)    # overfitting!
#'   summary(fit)
#'   z <- cbind( z, predict(fit, newdata = z, se = TRUE))
#'   head(z)
#'   (p <-  xyplot( income ~ education| gender, z, groups = type,
#'                  subscripts = T,
#'                  fit = z$fit,
#'                  lower = z$fit - z$se,
#'                  upper = z$fit + z$se,
#'                  layout = c(1,3),
#'                  auto.key = list(space='right', lines= T, cex = 1.5)))
#'
#'   p + glayer(gpanel.fit(...))
#'   trellis.focus()
#'   panel.identify(labels= z$occupation)
#'   trellis.unfocus()
#'   z$type2 <- with( z, reorder(type,education, mean, na.rm=T))
#'   gd(3)
#'   (p <-  xyplot( income ~ education| type2, z, groups = gender,
#'                  subscripts = T,
#'                  fit = z$fit,
#'                  lower = z$fit - z$se,
#'                  upper = z$fit + z$se,
#'                  layout = c(1,3),
#'                  par.strip.text = list(cex = 2),
#'                  auto.key = list(space='right', lines= T, cex = 1.5)))
#'
#'   p + glayer( gpanel.fit(...))
#'   trellis.focus()
#'   panel.identify(labels= z$occupation)
#'   trellis.unfocus()
#'
#'   ###
#'   ### With panels^2
#'   ### need to remove 'col = col.line'
#'   ###
#'
#'   z <- Prestige
#'   z$occ <- rownames(Prestige)
#'   z$gender <- with(z, cut( women, c(-1,15,50,101),labels = c("Male","Mixed","Female")))
#'   z$type2 <- with( z, reorder(type,education, mean, na.rm=T))
#'   tab(z, ~ gender + type2)
#'   z <- fillin( z, ~ education + type + gender, xpd = 1.1)
#'   fit <- lm( income ~ (education+I(education^2)+I(education^3) )* type * gender,
#'              z, na.action = na.exclude)    # overfitting!
#'   summary(fit)
#'   z <- cbind( z, predict(fit, newdata = z, se = TRUE))
#'   head(z)
#'   (p <-  xyplot( income ~ education| gender*type, z,
#'                  subscripts = T,
#'                  fit = z$fit,
#'                  labels = z$occ,
#'                  lower = z$fit - z$se,
#'                  upper = z$fit + z$se,
#'                  auto.key = list(space='right', lines= T, cex = 1.5)))
#'
#'   p + layer( gpanel.fit(...))
#'   p + layer( gpanel.fit(..., col = 'black', alpha = .1)) + layer(gpanel.labels(...))
#' }
#' @export
gpanel.fit <-
  function(x, y, fit, lower, upper,
           subscripts, ..., type, group.number, alpha, col, col.line, col.symbol, border = F, font, fontface)
  {
    if( !missing(fit)) {
      if( missing(col) ) col <- 'blue'
      if( !missing(group.number)) {
        col <- col.line
      }
      if( !missing(subscripts) ) {
        fit <- fit[subscripts]
      }
      dd <- data.frame( x=x, fit = fit)
      dd <- dd[order(dd$x),]
      panel.lines( dd$x, dd$fit, ..., col = col, type = 'l')
    }
    if( !missing(lower)){
      if( missing(alpha) || alpha == 1) alpha <- .3
      if( missing(col) ) col <- 'blue'
      if( !missing(group.number)) {
        col <- col.symbol
      }
      if( !missing(subscripts) ) {
        upper <- upper[subscripts]
        lower <- lower[subscripts]
      }
      dd <- data.frame( x=x, lower = lower, upper = upper)
      dd <- dd[order(dd$x),]
      panel.polygon( c(dd$x, rev(dd$x)),c(dd$upper, rev(dd$lower)),
                     border = border, col = col, alpha = alpha,...)
    }
    #  panel.polygon(c(dd$x, rev(dd$x)), c(dd$upper, rev(dd$lower)), col = col, alpha = alpha, ...)
  }
#' @describeIn gpanel.fit to be used with 'layer' -- but, actually, identical to 'glayer'
#' @export
panel.fit <- gpanel.fit
#' @describeIn gpanel.fit identical to gpanel.fit but kept for backward compatibility
#' @export
gpanel.band <- gpanel.fit
#' gpanel.labels: shows all labels
#'
#' This is an experiment in writing a function that can be
#' called via layer or glayer without further complications
#' e.g. \code{xyplot(......,labels = rownames(data)) + layer(gpanel.labels(...))}
#' or  \code{xyplot(....., labels = rownames(data), subscripts = T) + glayer(gpanel.labels(...))}.
#' For selected labels see the examples with \code{\link{trellis.focus}} and \code{\link{panel.identify}}
#' @param x,y position of labels, usually supplied through panel call
#' @param labels default is rownames of data
#' @param ... NOTE: may specify anything you don't want passed through ...
#' @examples
#' \dontrun{
#' trellis.focus()
#' panel.identify(labels = rownames(data),rot=-15,col = col.symbol, etc.)
#' trellis.unfocus()
#' }
#' @param labels to display
#' @describeIn gpanel.fit
#'
# NOTE: Also include names of anything you DON'T want passed.
#
#' @export
gpanel.labels <-
  function (x, y, labels , subscripts, ...)
  {

    if(!missing(subscripts)) {
      labels <- labels[subscripts]
    }
    panel.text(x, y, labels, ...)
  }
#' Fill in values
#'
#' Fills in values in gaps between observed predictor values to help produce
#' a smooth graph of predicted values versus predictor values.
#'
#' @param data data frame with values of x that need filling in
#' @param n number of additional points over range of predictor (default 200)
#' @param form formula idenfying variable x to fill in and grouping variables, g1, g2, etc.
#'       using syntax: \code{~ x + g1 + g2} (the variable to fill in comes first)
#' @param xpd expansion beyond range of predictor (default 1.0, i.e. no expansion)
#' @describeIn gpanel.fit
#' @export
fillin <- function(data, form, n = 200, xpd = 1.0) {
  #       levels of g
  #
  # n:    number of values to fill in across the entire range
  # xpd:  proportional range expansion

  sel.mf <- model.frame(form, data, na.action = na.pass)
  ret <- sel.mf
  xmina <- min(sel.mf[[1]], na.rm = T)
  xmaxa <- max(sel.mf[[1]], na.rm = T)
  xmida <- (xmaxa+xmina)/2
  xrana <- (xmaxa-xmina)/2
  xmina <- xmida - xpd*xrana
  xmaxa <- xmida + xpd*xrana
  if (ncol(sel.mf) ==1 ) {
    ret$xmin <- min(sel.mf[[1]],na.rm=T)
    ret$xmax <- max(sel.mf[[1]],na.rm=T)
  } else {
    grps <- lapply( sel.mf[2:ncol(sel.mf)], function(x)
      as.numeric(as.factor(x)))
    grps <- factor( do.call(paste,grps))
    ret$xmin <- capply( sel.mf[[1]], grps, min, na.rm=T)
    ret$xmax <- capply( sel.mf[[1]], grps, max, na.rm=T)
  }
  summ <- up(ret, formula(ret)[-2])  # has every unique combination of g's
  added.vals <- summ[ rep(1:nrow(summ), n),]
  added.vals[[names(ret)[1]]] <- rep( seq(xmina, xmaxa, length.out = n), each = nrow(summ))
  x <- added.vals[[names(ret)[1]]]
  x.min <- added.vals[['xmin']]
  x.max <- added.vals[['xmax']]
  x.mid <- (x.min + x.max)/2
  x.ran <- (x.max - x.min)/2
  ok <- (x <= (x.mid + xpd*x.ran)) & (x >= (x.mid - xpd*x.ran))
  merge( data, added.vals[ok,], all = T)
}
# zd <- data.frame( x = 1:10,
#                   g = rep(c('a','b'), each = 5),
#                   y = 1:10,
#                   g2 = LETTERS[c(1,1,1,2,2,2,3,3,3,3)])
# zd
# (ff <- fillin(zd, ~ x + g, n = 10, xpd = 1.2))
# tab(ff, ~ g)
# barchart
# ?barchart
# library(latticeExtra)
# zz <- fillin(zd, ~ x + g+g2, n =100)
# xyplot( x ~ g|g2 ,zd, cex = 2) +
#   xyplot( x ~ g|g2, fillin(zd, ~ x + g+g2, n =100),col='red')
#
# as.data.frame(tab(zd, ~g))
# gpanel.fit ####
#   Adds error bands. Uses two additional arguments:
#   'lower' and 'upper' for the lower and upper values of the bands
#   There are, provisionally, two versions:
#   panel.band which is minimal. To get group colors, specify col = col.symbol
#   and panel.band2 based on panel.smoother in latticeExtra that does a lot
#   more work with colors, etc.
#   We'll see which approach works best.
#
#FUTURE: gpanel.model
# gpanel.model( x, y, model, se = T)
# will combine gpanel.fit and gpanel.band
