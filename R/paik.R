#' Paik-Agresti diagrams
#'
#' Construct a Paik-Agresti diagram adapted and enhanced from the
#' \code{asbio::paik} function.
#'
#' @param formula	A two sided formula, e.g. Y ~ X1 + X2, with
#'   cross-classified categorical variables. The second
#'   explanatory variable, i.e. X2, is used as the trace variable
#'   whose levels are distinguished in the graph with different
#'   colors. Interactions and nested terms are not allowed.
#' @param counts A vector of counts for the associated
#'   categorical variables in formula. The variable 'Freq' is used
#'   if it exists in the data frame, data
#' @param resp.lvl The level in Y of primary interest. See
#'   example below.
#' @param data
#' 	 Dataframe containing variables in formula.
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
#'   Logical, indicating whether or not the words "Marginal prop"
#'   should printed in the graph above the dotted line indicating
#'   marginal proportions.
#' @param alpha transparency for circles expressed in hexadecimal, 
#'   e.g. 'AA' or 'FF' for no transparency. Default: '66'
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
#' 
#' data(death.penalty) # from Agresti 2012 
#' 
#' op <- par(mfrow=c(1,2), mar=c(4,4,0.1,0.1))
#' paik(verdict ~ d.race + v.race, counts = count, data = death.penalty, 
#' leg.title = "Victims race", xlab = "Defendants race", 
#' ylab = "Proportion receiving death penalty")
#' par(mar=c(4,2,0,2))
#' paik(verdict ~ v.race + d.race, counts = count, data = death.penalty, 
#' xlab = "Victims race", leg.title = "Defendants race",leg.loc="topleft", 
#' ylab = "", yaxt = "n")
#' par(op)
#' 
#' @export 
paik <- function (formula, counts, resp.lvl = 2, data, circle.mult = 1, 
    xlab = NULL, ylab = NULL, leg.title = NULL, leg.loc = NULL, 
    show.mname = FALSE,  
    col = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", 
            "#A6761D", "#666666"), 
    alpha = '66', ...) 
{
    draw.circle <- function (x, y, radius, nv = 100, border = NULL, col, lty = 1, 
              density = NULL, angle = 45, lwd = 1, alpha = '66') 
    {
        # copied from plotrix
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
  
    # fix for R version 4
    data[] <- lapply(data, function(x) if(is.character(x)) factor(x) else x)
    vars <- as.character(attr(terms(formula), "variables")[-1])
    cond.var = vars[3]
    rv <- data[, names(data) == vars[1]]
    cv <- data[, names(data) == cond.var]
    ov <- vars[vars != vars[1] & vars != cond.var]
    or <- data[, names(data) == ov]
    new.formula <- formula(counts ~ rv + ov + cv) # not used?
    cl <- levels(data[, names(data) == cond.var])
    ol <- levels(data[, names(data) == ov])
    rl <- levels(data[, names(data) == vars[1]])
    if(is.null(data$count)) data$count <- data$Freq
    x <- xtabs(count ~ rv + or + cv, data = data)
    xm <- xtabs(count ~ rv + or, data = data)
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
        radii[i], col = circle.col[i], alpha = alpha)
    if (length(pts) == 2) {
        for (i in 1:length(cl)) {
            segments(txm[i, ][1], tpropm[i, ][1], txm[i, ][2], 
                tpropm[i, ][2], col = circle.col[i])
        }
        segments(pts[1], m.prop[1], pts[2], m.prop[2], lty = 2, 
            lwd = 2)
    }
    if (length(pts) == 3) {
        for (i in 1:length(cl)) {
            segments(txm[i, ][1], tpropm[i, ][1], txm[i, ][2], 
                tpropm[i, ][2], col = circle.col[i])
            segments(txm[i, ][2], tpropm[i, ][2], txm[i, ][3], 
                tpropm[i, ][3], col = circle.col[i])
        }
        segments(pts[1], m.prop[1], pts[2], m.prop[2], lty = 2, 
            lwd = 2)
        segments(pts[2], m.prop[2], pts[3], m.prop[3], lty = 2, 
            lwd = 2)
    }
    if (length(pts) == 4) {
        for (i in 1:length(cl)) {
            segments(txm[i, ][1], tpropm[i, ][1], txm[i, ][2], 
                tpropm[i, ][2], col = circle.col[i])
            segments(txm[i, ][2], tpropm[i, ][2], txm[i, ][3], 
                tpropm[i, ][3], col = circle.col[i])
            segments(txm[i, ][3], tpropm[i, ][3], txm[i, ][4], 
                tpropm[i, ][4], col = circle.col[i])
        }
        segments(pts[1], m.prop[1], pts[2], m.prop[2], lty = 2, 
            lwd = 2)
        segments(pts[2], m.prop[2], pts[3], m.prop[3], lty = 2, 
            lwd = 2)
        segments(pts[3], m.prop[3], pts[4], m.prop[4], lty = 2, 
            lwd = 2)
    }
    if (length(pts) == 5) 
        stop("Number of rows in table must be less than 5")
    points(x.loc, y.loc, pch = 19, cex = 0.6)
    legend(ifelse(is.null(leg.loc), "topright", leg.loc), 
        pch = rep(21, length(cl)), pt.bg = col, bg = "white", 
        pt.cex = 1.5, title = ifelse(is.null(leg.title), cond.var, 
            leg.title), legend = cl)
    degree <- diff(m.prop)/diff(pts) * 360
    if (show.mname == TRUE) 
        text(mean(pts), mean(m.prop) + 0.03 * max(y.loc), srt = degree, 
            "Marginal prop.")
    res <- invisible(list(marginal.prop = r.prop, group.prop = r.sum/sum(r.sum)))
    invisible(res)
}
