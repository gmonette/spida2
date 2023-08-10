#' Paik-Agresti diagrams
#'
#' Construct a Paik-Agresti diagram adapted and enhanced from the
#' \code{asbio::paik} function.
#'
#' @param formula	A two sided formula, e.g. Y ~ X1 + X2, with
#'   cross-classified variables that will be treated as categorical variables. 
#'   The levels of the first explanatory variable, `X1` are displayed along
#'   the x-axis and the second
#'   explanatory variable, `X2`, is used as the conditioning variable
#'   whose levels are distinguished in the graph with different
#'   colors. Interactions and nested terms are not allowed.
#' @param data
#' 	 data.frame containing variables in formula. If rows of the the data frame should
#' 	 represent more than one observation, a variable named 'Freq' must
#' 	 record the frequency of each row, or the frequencies must be given through
#' 	 the 'counts' parameter below.  
#' @param counts A vector of counts for the associated
#'   categorical variables in formula. The variable 'Freq' is used
#'   if it exists in the data frame, \code{data}
#' @param resp.lvl The level in Y of primary interest. See
#'   example below.
#' @param circle.mult	
#'   Multiplier for circle radii in the diagram.
#' @param xlab	
#'   X-axis label. By default this is defined as the categories
#'   in the first explanatory variable, X1.
#' @param ylab	
#'   Y-axis label. By default these will be proportions with
#'   respect to the specified level of interest in the response.
#' @param leg.title	
#'   Legend title. By default the conditioning variable name.
#' @param leg.loc	
#'   Legend location. A legend location keyword; "bottomright",
#'   "bottom", "bottomleft", "left", "topleft", "top",
#'   "topright", "right" or "center".
#' @param show.mname	
#'   Logical, indicating whether or not the words "Overall proportion"
#'   should be printed in the graph above the dotted line indicating
#'   marginal proportions.
#' @param col
#'   list of colors for conditional levels. Default: the 8 colours
#'   of the 'Dark2' palette of RColorBrewer.
#' @param alpha transparency for circles expressed in hexadecimal, 
#'   e.g. 'AA' or 'FF' for no transparency. Default: '66'
#' @param cex, cex.cond, cex.marg cex for points, default: 1, default for
#'   cex.cond and cex.marg is cex
#' @param pch.marg pch for marginal points, default 15
#' @param pch.cond pch for conditional points, default 19
#' @param marginal logical, show marginal relationship, default: TRUE
#' @param lwd lwd for line segments, default: 2
#' @param lwd.marg lwd for marginal (overall) line segments. Default: lwd 
#' @param lwd.cond lwd for conditional line segments. Default: lwd
#' @param raise.prop proportion of vertical height of plot by which to 
#'   raise the "Overall Proportion" label. Default: 0.03
#' @param ...	
#'   Additional arguments from plot. Especially useful to 
#'   provide \code{ylim} if needed.
#' @author Ken Aho, modified by Georges Monette based on an idea by Fan Zhu
#' @references
#'   Agresti, A. (2012) Categorical Data Analysis, 3rd edition.
#'   New York. Wiley.
#'   Paik M. (1985) A graphical representation of a three-way
#'   contingency table: Simpson's paradox and correlation.
#'   American Statistician 39:53-54.
#' @examples
#' data(death.penalty) # from Agresti 2012 
#' print(death.penalty)
#' op <- par(mfrow=c(1,2), mar=c(4,4,0.1,0.1))
#' paik(verdict ~ d.race + v.race, 
#'    counts = death.penalty$count,
#'    data = death.penalty, 
#'    leg.title = "Victims' race", xlab = "Defendants' race", 
#'    ylab = "Proportion receiving death penalty")
#' par(mar=c(4,2,0,2))
#' paik(verdict ~ v.race + d.race, counts = death.penalty$count, data = death.penalty, 
#'      xlab = "Victims' race", leg.title = "Defendants' race",leg.loc="topleft", 
#'      ylab = "", yaxt = "n")
#' paik(am ~ gear + carb, mtcars)
#' par(op)
#' @importFrom gplots col2hex
#' @export 
paik <- function (formula, data, counts, resp.lvl = 2,  circle.mult = 1, 
    xlab = NULL, ylab = NULL, leg.title = NULL, leg.loc = NULL, 
    show.mname = TRUE,  
    col = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", 
            "#A6761D", "#666666"), 
    alpha = '66', 
    marginal = TRUE,
    cex = 1,
    cex.cond = cex,
    cex.marg = cex,
    pch.cond = 19,
    pch.marg = 15,
    lwd = 2,
    lwd.marg = lwd,
    lwd.cond = lwd,
    lwd.circle = lwd,
    raise.prop = .03,
    ...) 


{
  disp <- function(...) {}
  
  allvars <- all.vars(formula)
  if(!missing(counts)) data$Freq <- counts
  if(is.null(data$Freq) ) data$Freq <- 1
  form <- paste0("Freq ~ ", paste(allvars, collapse = '+'))
  disp(form)
  form <- as.formula(form)
  data <- as.data.frame(tab_(form, data))
  
  
  col <- gplots::col2hex(col)
  col <- paste0(col,alpha)
  draw.circle <- function (x, y, radius, nv = 100, border = NULL, col, lty = 1, 
                           density = NULL, angle = 45, lwd = 1, alpha = '66') 
  {  # copied from plotrix
    
    getYmult <- function () 
    {
      if (dev.cur() == 1) {
        warning("No graphics device open.")
        ymult <- 1
      }
      else {
        xyasp <- par("pin")
        xycr <- diff(par("usr"))[c(1, 3)]
        ymult <- xyasp[1]/xyasp[2] * xycr[2]/xycr[1]
      }
      return(ymult)
    }
    col <- ifelse(nchar(col) == 7, paste0(col,alpha), col)
    xylim <- par("usr")
    plotdim <- par("pin")
    ymult <- getYmult()
    angle.inc <- 2 * pi/nv
    angles <- seq(0, 2 * pi - angle.inc, by = angle.inc)
    if (length(col) < length(radius)) 
      col <- rep(col, length.out = length(radius))
    for (circle in 1:length(radius)) {
      xv <- cos(angles) * radius[circle] + x
      yv <- sin(angles) * radius[circle] * ymult + y
      polygon(xv, yv, border = border, col = col[circle], lty = lty, 
              density = density, angle = angle, lwd = lwd)
    }
    invisible(list(x = xv, y = yv))
  }
  
  # fix for R version 4: can't assume character variables are factors!
  data[] <- lapply(data, function(x) if(is.character(x)) factor(x) else x)
  
  vars <- as.character(attr(terms(formula), "variables")[-1])
  cond.var = vars[3]
  
  
  rv <- as.factor(data[, names(data) == vars[1]])  # response categories
  cv <- as.factor(data[, names(data) == cond.var])  # conditioning categories
  ov <- vars[vars != vars[1] & vars != cond.var]    # predictor categories(s)
  or <- data[, names(data) == ov]
  # new.formula <- formula(counts ~ rv + ov + cv) # not used?
  cl <- levels(data[, names(data) == cond.var])
  ol <- levels(data[, names(data) == ov])
  rl <- levels(data[, names(data) == vars[1]])
  # if(!missing(counts)) data$count <- counts
  # if(is.null(data$count)) {
  #   if(is.null(data$Freq)) {
  #     data$count <- 1 
  #   } else {
  #     data$count <- data$Freq
  #   }
  # }
  
  x <- xtabs(Freq ~ rv + or + cv, data = data)  # conditional
  xm <- xtabs(Freq ~ rv + or, data = data)      # marginal
  
  # marginal proportions
  m.prop <- apply(xm, 2, function(x) x[resp.lvl]/sum(x))
  
  r.prop <- matrix(nrow = length(cl), ncol = length(ol), dimnames = list(paste(cond.var, 
                                                                               cl, sep = "."), paste(ov, ol, sep = ".")))
  r.sum <- matrix(nrow = length(cl), ncol = length(ol), dimnames = list(paste(cond.var, 
                                                                              cl, sep = "."), paste(ov, ol, sep = ".")))
  for (i in 1:length(cl)) {
    tab <- x[, , cl = cl[i]]
    r.prop[i, ] <- apply(tab, 2, function(x) x[resp.lvl]/sum(x))
    r.sum[i, ] <- apply(tab, 2, sum)
  }
  b <- barplot(1:(length(ol) + 2), plot = FALSE)
  pts <- b[-c(1, (length(ol) + 2))]
  y.loc <- stack(as.data.frame(r.prop))[, 1]
  x.loc <- rep(pts, each = length(cl))
  temp.ylab <- bquote(paste("Proportion  ", .(vars[1]), 
                            "  =  ", .(rl[resp.lvl])))
  p <- plot(x.loc, y.loc, xlim = c(1.5, (length(ol) + 1.5)), 
            type = "n", xaxt = "n", xlab = ifelse(is.null(xlab), 
                                                  ov, xlab), ylab = ifelse(is.null(ylab), eval(temp.ylab), 
                                                                           ylab), ...)
  grid(p)
  axis(1, at = pts, labels = ol)
  tprop <- stack(as.data.frame(t(r.prop)))[, 1]
  tx <- rep(pts, length(cl))
  tpropm <- matrix(ncol = length(pts), nrow = length(cl), data = tprop, 
                   byrow = TRUE)
  txm <- matrix(ncol = length(pts), nrow = length(cl), data = tx, 
                byrow = TRUE)
  circle.col <- rep(col[1:length(cl)], length(ol))
  radii <- r.sum/sum(r.sum)
  radii <- stack(as.data.frame(radii))[, 1] * circle.mult * 0.4
  for (i in 1:length(radii)) draw.circle(x.loc[i], y.loc[i], 
                                         radii[i], col = circle.col[i], alpha = alpha,
                                         lwd = lwd.circle)
  if (length(pts) == 2) {
    for (i in 1:length(cl)) {
      segments(txm[i, ][1], tpropm[i, ][1], txm[i, ][2], 
               tpropm[i, ][2], col = circle.col[i], lwd = lwd.cond)
    }
    if(marginal)
      segments(pts[1], m.prop[1], pts[2], m.prop[2], lty = 2, 
               lwd = lwd.marg)
  }
  if (length(pts) == 3) {
    for (i in 1:length(cl)) {
      segments(txm[i, ][1], tpropm[i, ][1], txm[i, ][2], 
               tpropm[i, ][2], col = circle.col[i], lwd = lwd.cond)
      segments(txm[i, ][2], tpropm[i, ][2], txm[i, ][3], 
               tpropm[i, ][3], col = circle.col[i], lwd = lwd.cond)
    }
    if(marginal){
      segments(pts[1], m.prop[1], pts[2], m.prop[2], lty = 2, 
               lwd = lwd.marg)
      segments(pts[2], m.prop[2], pts[3], m.prop[3], lty = 2, 
               lwd = lwd.marg)
    }
  }
  if (length(pts) == 4) {
    for (i in 1:length(cl)) {
      segments(txm[i, ][1], tpropm[i, ][1], txm[i, ][2], 
               tpropm[i, ][2], col = circle.col[i], lwd = lwd.cond)
      segments(txm[i, ][2], tpropm[i, ][2], txm[i, ][3], 
               tpropm[i, ][3], col = circle.col[i], lwd = lwd.cond)
      segments(txm[i, ][3], tpropm[i, ][3], txm[i, ][4], 
               tpropm[i, ][4], col = circle.col[i], lwd = lwd.cond)
    }
    if(marginal) {
      
      segments(pts[1], m.prop[1], pts[2], m.prop[2], lty = 2, 
               lwd = lwd.marg)
      segments(pts[2], m.prop[2], pts[3], m.prop[3], lty = 2, 
               lwd = lwd.marg)
      segments(pts[3], m.prop[3], pts[4], m.prop[4], lty = 2, 
               lwd = lwd.marg)
    }  
  }
  if (length(pts) == 5) 
    stop("Number of rows in table must be less than 5")
  points(x.loc, y.loc, pch = pch.cond, cex = cex.cond)
  if(marginal) {
    points(pts, m.prop, pch = pch.marg, cex = cex.marg )
  }
  if(!marginal){
    
    legend(ifelse(is.null(leg.loc), "topright", leg.loc), 
           pch = rep(21, length(cl)), pt.bg = col, bg = "white", 
           pt.cex = 1.5, title = ifelse(is.null(leg.title), cond.var, 
                                        leg.title), legend = cl)
  } else {
    cl <- c(cl, 'Overall')
    legend(ifelse(is.null(leg.loc), "topright", leg.loc), 
           pch = c(rep(21, length(cl)-1),pch.marg) , pt.bg = c(col,'black'), bg = "white", 
           pt.cex = 1.5, title = ifelse(is.null(leg.title), cond.var, 
                                        leg.title), legend = cl)
    
  }
  pin <- par('pin')
  
  degree <- diff(m.prop)/diff(pts)/(pin[2]/pin[1]) * 360
  degree <- diff(m.prop)/diff(pts)/(pin[1]/pin[2]) * 360
  disp(m.prop)
  disp(pts)
  degree <- diff(m.prop)/diff(pts) * 360
  degree <- atan2(diff(m.prop)*pin[2], diff(pts)*pin[1]) * (180/pi)
  degree <- atan2(diff(m.prop), diff(pts)) * (180/pi)
  usr <- par('usr')
  ud <- c(usr[2] - usr[1], usr[4] - usr[3])
  degree <- atan2(diff(m.prop)*ud[2]*pin[2], diff(pts)*ud[1]*pin[1]) * (180/pi)
  degree <- 45
  degree <- atan2((diff(m.prop)/ud[2])*pin[2], (diff(pts)/ud[1])*pin[1]) * (180/pi)
  disp(degree)
  if (marginal & show.mname == TRUE & length(pts) ==2) 
    text(mean(pts), mean(m.prop) + raise.prop * (usr[4] - usr[3]), srt = degree, 
         "Overall proportion")
  res <- invisible(list(marginal.prop = r.prop, group.prop = r.sum/sum(r.sum)))
  invisible(res)
}

if(FALSE){
  library(spida2)
  
  paik(verdict ~ d.race + v.race,
       counts = death.penalty$count,
       data = death.penalty,
       leg.title = "Victims race", xlab = "Defendants race",
       show.mname = TRUE,
       marginal = TRUE,
       lwd = 3, raise.prop = .03,
       ylab = "Proportion receiving death penalty")
  
  paik(verdict ~ d.race + v.race,
       circle.mult = .5,
       counts = death.penalty$count,
       data = death.penalty,
       leg.title = "Victims race", xlab = "Defendants race",
       ylab = "Proportion receiving death penalty")
  
  paik(cyl ~ vs + carb, mtcars)
  paik(am ~ vs + carb, mtcars)
 
  par(mar=c(4,2,0,2))
  paik(verdict ~ v.race + d.race, counts = death.penalty$count, data = death.penalty,
       xlab = "Victims race", leg.title = "Defendants race",leg.loc="topleft",
       ylab = "", yaxt = "n")
  paik(verdict ~ v.race + d.race, data = death.penalty,
       xlab = "Victims race", leg.title = "Defendants race",leg.loc="topleft",
       ylab = "", yaxt = "n")
  
  z <- death.penalty
  z$Freq <- z$count
  z$count <- NULL
  paik(verdict ~ v.race + d.race, data = death.penalty,
       xlab = "Victims race", leg.title = "Defendants race",leg.loc="topleft",
       ylab = "", yaxt = "n")
}
#' @describeIn paik version for numerical Y value
#' @examples
#' dd <- death.penalty  
#' dd$Freq <- dd$count
#' dd$ver <- as.numeric(dd$verdict == 'Y')
#' paik2(ver ~ d.race + v.race, dd, cex = 2, circle.factor = .2)  
#' 
#' d <- mtcars
#' paik2(mpg ~ cyl + gear, d, cex = 3)  
#' @export
paik2 <- function (formula, data, counts, resp.lvl = 2,  circle.mult = 1, 
                  xlab = NULL, ylab = NULL, leg.title = NULL, leg.loc = NULL, 
                  mtext = "Overall proportion",  
                  col = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", 
                          "#A6761D", "#666666"), 
                  alpha = '66', 
                  marginal = TRUE,
                  cex = 1,
                  cex.cond = cex,
                  cex.marg = cex,
                  pch.cond = 19,
                  pch.marg = 15,
                  lwd = 2,
                  lwd.marg = lwd,
                  lwd.cond = lwd,
                  lwd.circle = lwd,
                  raise.prop = .03,
                  circle.factor = .2,
                  ylim,
                  ...) 


{
  disp <- function(...) {}

  allvars <- all.vars(formula)
  yvname <- allvars[1]
  xvname <- allvars[2]
  cvname <- allvars[3]
  
  # Weighted averages
  
  if(is.null(data$Freq)) data$Freq <- 1
  if(!missing(counts)) data$Freq <- counts
  
  dd <- data
  dd$yvar <- data[[yvname]]
  if(is.factor(dd$yvar)) stop("Use paik if y is a factor")
  dd$xvar <- as.factor(data[[xvname]])
  dd$cvar <- as.factor(data[[cvname]])
  
  dd$ymean <- capply(dd, ~ xvar+cvar, with, weighted.mean(yvar, Freq) )
  
  dd$ymarg <- capply(dd, ~ xvar, with, weighted.mean(yvar, Freq) )
  
   
  col <- gplots::col2hex(col)
  col <- paste0(col,alpha)
  draw.circle <- function (x, y, radius, nv = 100, border = NULL, col, lty = 1, 
                           density = NULL, angle = 45, lwd = 1, alpha = '66') 
  {  # copied from plotrix
    
    getYmult <- function () 
    {
      if (dev.cur() == 1) {
        warning("No graphics device open.")
        ymult <- 1
      }
      else {
        xyasp <- par("pin")
        xycr <- diff(par("usr"))[c(1, 3)]
        ymult <- xyasp[1]/xyasp[2] * xycr[2]/xycr[1]
      }
      return(ymult)
    }
    col <- ifelse(nchar(col) == 7, paste0(col,alpha), col)
    xylim <- par("usr")
    plotdim <- par("pin")
    ymult <- getYmult()
    angle.inc <- 2 * pi/nv
    angles <- seq(0, 2 * pi - angle.inc, by = angle.inc)
    if (length(col) < length(radius)) 
      col <- rep(col, length.out = length(radius))
    for (circle in 1:length(radius)) {
      xv <- cos(angles) * radius[circle] /ymult + x
      yv <- sin(angles) * radius[circle]  + y
      polygon(xv, yv, border = border, col = col[circle], lty = lty, 
              density = density, angle = angle, lwd = lwd)
    }
    invisible(list(x = xv, y = yv))
  }
  
  xlevs <- levels(dd$xvar)
  clevs <- levels(dd$cvar)
  
  b <- barplot(1:(length(xlevs) + 2), plot = FALSE)
  x.pts <- b[-c(1, (length(xlevs) + 2))]
  
  dd$xvals <- x.pts[dd$xvar]
  
  ddcond <- up(dd, ~ xvar + cvar)
  ddmarg <- up(dd, ~ xvar)
  
  radii <- ddcond$Freq/sum(ddcond$Freq)
  ylims <- range(ddcond$ymean)
  radii <- radii * circle.factor * diff(ylims) / max(radii)
  ddcond$radii <- radii
  
  limits <- if(missing(ylim)) range(c(ddcond$ymean + radii,ddcond$ymean-radii )) else ylim
  # disp(limits)
  p <- plot(ddcond$xvals, ddcond$ymean, xlim = c(1.5, (length(xlevs) + 1.5)), 
            type = "n", xaxt = "n", 
            ylim = limits,
            xlab = ifelse(is.null(xlab), xvname, xlab), 
            ylab = ifelse(is.null(ylab), yvname, ylab), ...)
  grid(p)
  axis(1, at = x.pts, labels = xlevs)
  
  with(ddcond,points(xvals,ymean, pch = 16))
  if(marginal) with(ddmarg, points(xvals, ymarg, pch = 18))
  
  for(i in seq_along(clevs)) {
    with(ddcond[ddcond$cvar==clevs[i],], 
         lines(xvals,ymean, col = col[i], lwd = lwd.cond ))
  }
  
  if(marginal) with(ddmarg, lines(xvals, ymarg, col = 'black', lwd = lwd.marg))
  
  for (i in 1:nrow(ddcond)) {
    with(ddcond[i,],
         draw.circle(
           xvals, ymean, radii, col = col[cvar],
           alpha = alpha, lwd = lwd.circle))
  }
  
  if(!marginal){
    legend(ifelse(is.null(leg.loc), "topright", leg.loc), 
           pch = rep(21, length(clevs)), pt.bg = col, bg = "white", 
           pt.cex = 1.5, title = ifelse(is.null(leg.title), cvname, 
                                        leg.title), legend = clevs)
  } else {
    clevs <- c(clevs, 'Overall')
    legend(ifelse(is.null(leg.loc), "topright", leg.loc), 
           pch = c(rep(21, length(clevs)-1),pch.marg) , pt.bg = c(col,'black'), bg = "white", 
           pt.cex = 1.5, title = ifelse(is.null(leg.title), cvname, 
                                        leg.title), legend = clevs)
    
  }
  
  
  pin <- par('pin')
  usr <- par('usr')
  ud <- c(usr[2] - usr[1], usr[4] - usr[3])
  degree <- atan2((diff(ddmarg$ymarg[1:2])/ud[2])*pin[2], (diff(ddmarg$xvals[1:2])/ud[1])*pin[1]) * (180/pi)
  if (marginal & !is.null(mtext)) 
    text(mean(ddmarg$xvals[1:2]), mean(ddmarg$ymarg[1:2]) + raise.prop * (usr[4] - usr[3]), srt = degree, 
         mtext)
  res <- list(data = data, cond = ddcond, marginal = ddmarg)
  invisible(res)
}
 