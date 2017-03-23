##
## panel functions:
##
## panel.fit incorporating panel.band
## panel.text
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
#' @param data data frame to be used to add additional values of numeric variable
#' @param form formula evaluated in data. The first term defines the variable with values to be filled in and the remaining terms define the variables to be used for grouping determining the minima and maxima within which values are added.
#' @param n the number of values to be added between the global mininum and maximum.
#'        Values falling outside conditional minima and maxima are culled. Default 200.
#' @param xpd expansion factor to add points beyond minima and maxima. Default 1.0.
#' @param y
#' @param fit fitted values of a model, generally passed through 'layer' from a call to 'xyplot': e.g. \code{xyplot( y ~ x, data, groups = g, fit = data$yhat, lower = with(data, yhat - 2*se), upper = with(data, yhat + 2*se), subscripts = T)}
#' @param lower,upper lower and upper limits of error bands, passed from main plotting function
#' @param subscripts subscripts, passed from main plotting function
#' @param dots passed from main plotting function
#' @param col passed from main plotting function
#' @param group.number , passed from main plotting function
#' @param alpha transparency, passed from main plotting function
#' @param col.symbol is used to control color when using 'groups'
#' @param border default = FALSE for panel.band.
#' @param font passed from main plotting function
#' @param fontface passed from main plotting function
#' @return The 'panel.bands', 'panel.fit', and 'panel.labels' functions are
#'         invoked for their graphical effect.
#' @author Georges Monette <georges@@yorku.ca>
#' @examples
#' \dontrun{
#'   library(spida2)
#'   library(latticeExtra)
#'   library(car)
#'   
#'   ### Show data, predicted values and confidence bands 
#'   ### using a 'predict' method that produces 'se'   
#'   ### -- showing predicted values only at points in the data
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
#'   p + glayer(panel.fit(...))
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
#'   z <- fillin(Prestige, ~education + type, xpd = 1.1) # fill in 'education' within 'type'
#'
#'   dim(z)
#'   dim(Prestige)
#'   z <- cbind(z, predict(fit, newdata = z, se = TRUE))
#'   head(z)
#'   gd(3, cex = 2, lwd = 2, alpha = .7)
#'   (p <-  xyplot( income ~ education, z, groups = type,
#'                  subscripts = T,
#'                  fit = z$fit,
#'                  lower = z$fit - z$se,
#'                  upper = z$fit + z$se,
#'                  auto.key = list(space='right', lines= T)))
#'   p + glayer( panel.fit(...))
#'   p + glayer( panel.fit(...,  alpha = .1))
#'   # Using spida2:gd() to get a ggplot-like appearance
#'   gd(3,lty=1,lwd=2)
#'   p + glayer( panel.fit(...,alpha = .1))
#'
#'   ###
#'   ###  Prediction with 'wald' for fits with predict
#'   ###  methods that don't produce se's, e.g. lme
#'   ###
#'   
#'   library(nlme)
#'   fit <- lme(mathach ~ (ses+I(ses^2)) * Sex * Sector, hs, random = ~ 1|school)
#'   summary(fit)
#'   pred <- expand.grid( ses = seq(-2,2,.1), Sex = levels(hs$Sex), Sector = levels(hs$Sector))
#'   head(pred)
#'   # merge pred with original data
#'   hs$from <- 'data'
#'   pred$from <- 'pred'
#'   dm <- merge(hs, pred, all = T)
#'   w <- wald(fit, getX(fit, data = dm)) # attaches data to wald.object so it can included in data frame
#'   w <- as.data.frame(w)
#'   head(w)
#'   library(latticeExtra)
#'   gd(pch = 1, alpha = .2)
#'   (p <- xyplot(mathach ~ ses | Sector, w, groups = Sex,
#'        auto.key = T, type = 'p',
#'        fit = w$coef,
#'        upper = with(w,coef+2*se),
#'        lower = with(w,coef-2*se),
#'        subscript = T))
#'   p + glayer( gpanel.fit(...,alpha=1))
#' 
#'   wald(fit, 'Sex')  # sig. overall effect of Sex
#'   wald( fit, ':Sex') # but no evidence of interaction with ses
#'   wald( fit, '\\^2') # nor of curvature
#' 
#'   # but we continue for the sake of illustration
#' 
#'   L <- Lform( fit, list( 0, 1, 2*ses, 0, Sex == 'Male', (Sex == 'Male')*2*ses), hs)
#'   L
#'   (ww <- wald ( fit, L ))
#'   wald.dd <- as.data.frame(ww, se = 2)
#'   head( wald.dd )
#' 
#'   xyplot(coef ~ ses | Sex, wald.dd, type = 'n',
#'       ylim = c(-5, 12),
#'       ylab = 'increase in mathach per unit increase in ses',
#'       fit = wald.dd$coef,
#'       upper = wald.dd$U2,
#'       lower = wald.dd$L2,
#'       subscripts = T) +
#'   layer(panel.fit(...))
#'
#'
#'   ###
#'   ###  Using panel.fit with no groups
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
#'   gd_(basecol = 'tomato4')  # Use 'gd_' to set parameters without groups
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
#'   p + glayer(panel.fit(...))
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
#'   p + glayer( panel.fit(...))
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
#'   p + layer( panel.fit(...))
#'   p + layer( panel.fit(..., col = 'black', alpha = .1)) + layer(panel.labels(...))
#'   
#'   
#'   
#' }
#' @export
panel.fit <-
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
#' @describeIn panel.fit to be used with 'layer' -- but, actually, identical to 'glayer'
#' @export
gpanel.fit <- panel.fit
#' @describeIn panel.fit identical to panel.fit but kept for backward compatibility
#' @export
panel.band <- panel.fit
#' panel.labels: shows all labels
#'
#' This is an experiment in writing a function that can be
#' called via layer or glayer without further complications
#' e.g. \code{xyplot(......,labels = rownames(data)) + layer(panel.labels(...))}
#' or  \code{xyplot(....., labels = rownames(data), subscripts = T) + glayer(panel.labels(...))}.
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
#' @describeIn panel.fit
#'
# NOTE: Also include names of anything you DON'T want passed.
#
#' @export
panel.labels <-
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
#' @describeIn panel.fit
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
# panel.fit ####
#   Adds error bands. Uses two additional arguments:
#   'lower' and 'upper' for the lower and upper values of the bands
#   There are, provisionally, two versions:
#   panel.band which is minimal. To get group colors, specify col = col.symbol
#   and panel.band2 based on panel.smoother in latticeExtra that does a lot
#   more work with colors, etc.
#   We'll see which approach works best.
#
#FUTURE: panel.model
# panel.model( x, y, model, se = T)
# will combine panel.fit and panel.band



##
## panel functions:
##
## panel.box
##
## designed to work easily with 'layer' and 'glayer' in latticeExtra
## Added: 2016 08 05
## Author: GM
##
#' Panel functions for predicted values and SE box
#'
#' With 'layer' and 'glayer' in 'latticeExtra', these functions can be used to easily generate fitted values and
#' error boxes that have a reasonable appearance whether a plot uses 'groups' or not.
#
#' @param data data frame to be used to add additional values of numeric variable
#' @param form formula evaluated in data. The first term defines the variable with values to be filled in and the remaining terms define the variables to be used for grouping determining the minima and maxima within which values are added.
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
#'   library(spida2)
#'   library(latticeExtra)
#'   library(car)
#'   Prestige$Education <- cut(Prestige$education,
#'          breaks = c(-Inf,8,12,+Inf), labels = c('< HS','some HS','some PS'))
#'   Prestige$Gender <- cut(Prestige$women, breaks = c(-Inf, 10, 50, +Inf),
#'          labels = c('Male','Mixed','Female'))
#'
#'   fit <- lm(income ~ Gender * Education, Prestige,
#'       na.action = na.exclude)
#'   pred <- with(Prestige, expand.grid(Education = levels(Education),
#'            Gender = levels(Gender)))
#'   pred <- cbind(pred, predict(fit, newdata = pred, se = TRUE))
#'   pred
#'   (p <- xyplot( fit ~ Education | Gender , pred,  
#'                 subscripts = T, width = .5,
#'                 fit = pred$fit,
#'                 labels = pred$Education,
#'                 lower = with(pred, fit - 2*se.fit),
#'                 upper = with(pred, fit + 2*se.fit)))
#'   gd(3,pch = 16)
#'   p + layer(panel.fit(...))
#'   (p <- update(p + layer(panel.box(...)), ylim = c(0, 15000), 
#'          auto.key = list(columns= 3)))
#'   ###### Need to fix handling of width with factor predictor   
#'   p + layer(panel.labels(...))
#'   
#' # NOTE: with groups use 'glayer' instead of 'layer' 
#' }
#' @export
panel.box <-
  function(x, y, se, lower = y - se, 
           upper = y + se, 
           width = if(is.numeric(x)) max(diff(sort(x)))/2 else 1/2,
           subscripts, ..., type, group.number, 
           alpha = .9, alpha.fit = alpha, alpha.box = .3,
           col, col.line, col.symbol, border = F, font, fontface)
  {
    circ <- function(x, y, width, n = 32) {
      th <- c(seq_len(n)*2*pi/n, NA)
      xyscale <- current.panel.limits()
      ytox <- diff(xyscale$ylim)/diff(xyscale$xlim)
      cbind(c(as.numeric(x) + width*sin(th)/2), c(y+width*ytox*cos(th)/2))
    }  
    if(missing(col)) col <- 'blue'
    if(!missing(group.number)) col <- col.symbol
    pmat <- do.call(rbind, lapply(seq_along(x), function(ii) circ(x[ii], y[ii], width = width)))
    panel.polygon(pmat, col = col, border = TRUE, alpha = alpha.fit, ...)
    # modify to allow width to vary 
    if(is.null(width)) width <- min(diff(x), na.rm = TRUE)/2
    if(!missing(lower)){
      if( missing(alpha) || alpha == 1) alpha <- .3
      if( !missing(subscripts) ) {
        upper <- upper[subscripts]
        lower <- lower[subscripts]
      }
      xx <- c(rbind(x,x,x,x,NA))
      xx <- xx + c(-1,-1,1,1,NA)*width/2  
      yy <- c(rbind(lower,upper,upper,lower,NA))
      panel.polygon(xx,yy,
                     border = border, col = col, alpha = alpha.box,...)
    }
  }

