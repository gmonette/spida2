
###
###  td
###

#' @name td
NULL


#' Set lattice parameters for multiple groups
#'
#' Easier alternative to using trellis.par.set and trellis.par.get to change lattice colors, line types, etc.
#'
#' Designed to easily set lattice parameters for multiple groups. Setting
#' parameters before calling the lattice function allows parameters to be used
#' consistently in the group key.
#'
#' 'td' calls 'trellis.device' and sets graphical parameters for
#' 'superpose.line' and 'superpose.symbol'. 'td' also initializes a new trellis
#' device with a white background if new = TRUE.
#'
#' 'gd' is similar to 'td' except that it uses a theme that resembles that of 'ggplot'
#'
#' @param n in 'gd' specifies the number of distinct colours to generate to
#' distinguish groups. 'gd' uses 'latticeExtra' to set defaults for a
#' ggplot2-like appearance. Default is n = 4
#' @param new If new = TRUE, open a new window, otherwise modify the existing
#' active window, if there is one.
#' @param col FIXME
#' @param lty FIXME
#' @param lwd FIXME
#' @param pch FIXME
#' @param cex FIXME
#' for each level of the groups variable
#' @param font FIXME
#' @param fill FIXME
#' @param col.line FIXME
#' not given
#' @param col.symbol FIXME
#' 'groups' not given
#' @param alpha FIXME
#' @param alpha.line FIXME
#' @param alpha.symbol graphical parameters for superpose.line and
#' superpose.symbol
#' @param len extend the length of parameters by recycling to length 'len'
#' @param long if TRUE generate a default combination of col, lty and pch with
#' length 42.
#' @param record If TRUE, set history to 'recording'. Caution: this can use a
#' lot of memory in .GlobalEnv.  Consider adding selected graphs to memory with
#' the 'Insert' key instead.
#' @param basecol FIXME
#' @param colsets FIXME
#' @param \dots FIXME
#' parameter: e.g. \code{plot.symbol=list(cex=2,col='red')}. Particular useful
#' for the cex, col, lty, lwd, alpha, pch parameters in plot.line and
#' plot.symbol.
#' @author Georges Monette
#' @concept lattice
#' @examples
#' td( lty = 1:7)   # sets line types for 7 groups
#' gd(7)            # sets line types for 7 groups using colors from RColorBrewer
#' td( plot.symbol = list(col = 'red', pch = 17))
#' gd_(col='blue')  # set main color to 'blue'
#' @export
td <- function(...) gd(..., gglike = FALSE)
# td <- function(
#     new = FALSE,
#     col=c("#0080ff",   "#ff00ff",   "darkgreen", "#ff0000" ,  "orange" ,   "#00ff00",   "brown" ),
#     lty=1:7, lwd=1,
# 	pch = 1:7, cex = 0.8, font = 1,
# 	long = FALSE,
# #    record = FALSE, # not supported in RStudio
#     basecol = NULL,
# 	  colsets = c('plot.symbol','plot.line','dot.symbol',
# 			'dot.line','cloud.3d','box.dot'),
#     ...) {
# 
#         # Modified for R: Oct. 10, 2004
# 	#
# 	# reset superpose.symbol and superpose.line so they have consistent length
# 	# equal to the max of input parameters:
# 	#			sps: cex, col, font, pch
# 	#			spl: lty, col, lwd
# 	# or to len
# 	#  This allows distinctive line styles, for example for 12 objects by using
# 	#  good lty's:     1,4,8,6,3,5,7
# 	#  and good col's: 1,8,6,3,4,5,2
# 	#
#   require(lattice)
#   aargs <- list(...)
#   if ( long ) {
#     col <- c(3,5,4,6,8,2)   # drop yellow
#     len <- 42 # generates 42 unique combinations of pch/lty and col
#   }
#   if(new) trellis.device(theme = col.whitebg, record = record, new = new)
#                                         # NOTE: fixed panel.superpose so lty argument
#                                         # is passed to points for type = 'b'
#   len <- max(len,length(col),length(lty),length(lwd),length(pch),length(cex),length(font))
#   spl <- trellis.par.get("superpose.line")
#   spl$alpha <- rep(alpha.line, length = len)
#   spl$lty <- rep(lty,length=len)
#   spl$col <- rep(col.line,length=len)
#   spl$lwd <- rep(lwd,length=len)
#   trellis.par.set("superpose.line",spl)
#   sps <- trellis.par.get("superpose.symbol")
#   sps$alpha <- rep( alpha.symbol, length = len)
#   sps$pch <- rep(pch, length=len)
#   sps$col <- rep(col.symbol, length=len)
#   sps$cex <- rep(cex, length=len)
#   sps$font <- rep(font, length=len)
#   sps$fill <- rep(fill, length=len)
# 
#   trellis.par.set("superpose.symbol",sps)
#   list(superpose.symbol = sps, superpose.line = spl)
#   if ( !is.null(basecol)) {
#     for ( ii in colsets ) {
#       tt <- trellis.par.get(ii)
#       tt$col <- basecol
#       trellis.par.set(ii,tt)
#     }
#   }
#   if ( length(aargs)){
#     tpg <- trellis.par.get()
#     for ( nn in names(aargs)){
#       for(mm in names(aargs[[nn]])){
#         tpg[[nn]][[mm]] <- aargs[[nn]][[mm]]
#       }
#     }
#     trellis.par.set(theme = tpg)
#   }
#   ret <- trellis.par.get()
#   invisible(ret[grep('superpose',names(ret))])
# }
#' Set lattice parameters for plotting by group
#'  
#' @param n number of groups for which to set colors, line types, etc. using RColorBrewer.
#' @examples
#' #   - setting colors for groups, i.e. 'superpose.symbol' in trellis.par.get():
#' gd(5)  # where 5 is the number of groups
#' gd(5, lwd = 2, lty = 1)
#' gd(5, col = brewer.pal(5,"Dark2"),cex = 1.5)
#'
#' # To set colors when not using groups (superpose = FALSE)
#' gd_(col='tomato4')
#' # changing the default color for lines and symbols
#' gd(plot.line=list(col='red',lwd=2),
#'            plot.symbol=list(col='blue', cex = 1.3))
#' # OR using superpose = FALSE
#' gd(superpose = FALSE, col = 'red', lwd = 2)
#' # OR using the utility function 'gd_':
#' gd_(col = 'red', lwd = 2)
#' #
#' # To set colors for lattice::barchart:
#' library(lattice)
#' gd(superpose.polygon = list(col=brewer.pal(4,'Paired'), border='black'))
#' barchart(Titanic, 
#'     auto.key=list(title = 'survived',
#'                   space = 'right',
#'                   reverse.rows = T), 
#'     horizontal = F)
#' #  For a complete list of elements that can be changed:
#' names(trellis.par.get())
#' # For a list of colors
#' colors()
#' grepv('pink',colors()) # types of pink
#' # Using magrittr
#' library(magrittr)
#' colors()  %>%  grepv('blue', .)  %>%
#'   pal  %>%
#'   as.data.frame %>%
#'   sortdf( ~ red)  %>%
#'   as.matrix  %>%
#'   divide_by(255)  %>%
#'   rgb  %>%
#'   pal
#' @concept lattice
#' @describeIn td with ggplot2-like theme
#' @export
gd <- function (n=8, pal = "Dark2",
            col = brewer.pal(n, pal), lty = 1:n, lwd = 1,
            pch = 19, cex = 1.4, font = 1, fill = "transparent",
            col.line = col, col.symbol = col,
            alpha = 1, alpha.line = alpha, alpha.symbol = alpha,
            len = n,
            # arguments to ggplot2like:
            h = c(0,360) + 15, l = 65, c = 100, h.start = 0, direction = 1,
            low = "#3B4FB8", high = "#B71B1A", space = "rgb",
            # trellis par set parameters for basecol:
            basecol = NULL,
            colsets = c("plot.symbol","plot.line", "dot.symbol",
                        "dot.line", "cloud.3d", "box.dot"),
            superpose = TRUE,
            # set ggplot2 like environment even if not first call
            gginit = FALSE,
            gglike = TRUE,
            # other arguments of form:
            #      plot.symbol = list( col = 'red', pch = 4)
            ...)
{
    # gd makes it easy to set graphical parameters,
    # i.e. col, lwd, lty, fill, font, cex, pch, alpha
    # for 'superpose.symbol' and 'superpose.line' used
    # for different groups in lattice
    # Note: fill works with pch 21:25
    #
    # gd can also be used to set other graphical parameters
    # by specifying the list in which they are set

    library(lattice)
    library(latticeExtra)
    library(RColorBrewer)
    
  
  if(FALSE){ # feature turned off GM:2020_02_13
    # If there are arguments other than n
    # AND all are of length no greater than 1, superpose is FALSE
    arglist <- as.list(match.call())[-1]
    arglist$n <- NULL
    if(length(arglist) > 0 && 
       pmax(sapply(arglist,length)) == 1 &&
       missing(superpose) &&
       missing(n)) superpose <- FALSE
  }
  aargs <- list(...)
    
    # ggplot2
    if(gglike) {
      if(is.null(lattice.options('gginit')[[1]]) | gginit == TRUE){
        lattice.options(gginit=TRUE)
        trellis.par.set(ggplot2like(n = n,h = h,l = l,c = c,
                                    h.start = h.start, direction = direction,
                                    low = low , high = high , space = space))
        lattice.options(ggplot2like.opts())
        aargs$superpose.polygon <- list(col=col,border='black')
        aargs$plot.polygon <- list(col=col[1],border='black')
      }
    }
    len <- max(len, length(col), length(lty), length(lwd), length(pch),
               length(cex), length(font))
    if (superpose) {
      spl <- trellis.par.get("superpose.line")
      spl$alpha <- rep(alpha.line, length = len)
      spl$lty <- rep(lty, length = len)
      spl$col <- rep(col.line, length = len)
      spl$lwd <- rep(lwd, length = len)
      trellis.par.set("superpose.line", spl)
      sps <- trellis.par.get("superpose.symbol")
      sps$alpha <- rep(alpha.symbol, length = len)
      sps$pch <- rep(pch, length = len)
      sps$col <- rep(col.symbol, length = len)
      sps$cex <- rep(cex, length = len)
      sps$font <- rep(font, length = len)
      sps$fill <- rep(fill, length = len)
      trellis.par.set("superpose.symbol", sps)
    } else { # use to set non-panel setting
      tt <- trellis.par.get()
      if ( !missing(col) ) {
        tt$plot.symbol$col <- col.symbol
        tt$plot.line$col <- col.line
      }
      if ( !missing(col.line)) {
        tt$plot.line$col <- col.line
      }
      if ( !missing(col.symbol)) {
        tt$plot.symbol$col <- col.symbol
      }
      if ( !missing(alpha) ) {
        tt$plot.symbol$alpha <- alpha.symbol
        tt$plot.line$alpha <- alpha.line
      }
      if ( !missing(alpha.line)) {
        tt$plot.line$alpha <- alpha.line
      }
      if ( !missing(alpha.symbol)) {
        tt$plot.symbol$alpha <- alpha.symbol
      }
      if ( !missing(lty)) {
        tt$plot.line$lty <- lty
      }
      if ( !missing(lwd)) {
        tt$plot.line$lwd <- lwd
      }
      if ( !missing(pch)) {
        tt$plot.symbol$pch <- pch
      }
      if ( !missing(cex)) {
        tt$plot.symbol$cex <- cex
      }
      if ( !missing(fill)) {
        tt$plot.symbol$fill <- fill
      }
      trellis.par.set(theme = tt)
    }
    if (!is.null(basecol)) {
      for (ii in colsets) {
        tt <- trellis.par.get(ii)
        tt$col <- basecol
        trellis.par.set(ii, tt)
      }
    }
    if (length(aargs)) {
      tpg <- trellis.par.get()
      for (nn in names(aargs)) {
        for (mm in names(aargs[[nn]])) {
          tpg[[nn]][[mm]] <- aargs[[nn]][[mm]]
        }
      }
      trellis.par.set(theme = tpg)
    }
    ret <- trellis.par.get()
    invisible(ret[grep("superpose", names(ret))])
  }
#' @describeIn td gd to set non-group parameters
#' @export
gd_ <- function(...) {
  gd(superpose = FALSE, ...)
}
#' @describeIn td trellis par set for superpose parameters
#' @export
tps <- function (...)
{
  args <- list(...)
  # disp(args)
  theme <- trellis.par.get()
  nams <- names(args)
  for(i in seq_along(args)) {
    # disp(nams[i])
    switch(nams[i],
         'lty' = {
           theme$superpose.line$lty <- args[[i]]
         },
         'col' = {
           theme$superpose.line$col <- args[[i]]
           theme$superpose.symbol$col <- args[[i]]
         }, 
          'symbol_col' = {
           
             theme$superpose.symbol$col <- args[[i]]
             
           },
         'line_col' = {
           theme$superpose.line$col <- args[[i]]
         },
         'lwd' = {
           theme$superpose.line$lwd <- args[[i]]
           # disp('lwd')
         },
         'alpha' = {
           theme$superpose.line$alpha <- args[[i]]
           theme$superpose.symbol$alpha <- args[[i]]
           # disp('alpha')
         },
         'symbol_alpha' = {
           theme$superpose.symbol$alpha <- args[[i]]
           # disp('did this too')
           
         },
         'line_alpha' = {
           theme$superpose.line$alpha <- args[[i]]
           
         },
         'font' = {
           theme$superpose.symbol$font <- args[[i]]
           
         },
         'fill' = {
           theme$superpose.symbol$fill <- args[[i]]
           
         },
         'cex' = {
           theme$superpose.symbol$cex <- args[[i]]
           
         },
         'pch' = {
           theme$superpose.symbol$pch <- args[[i]]
         },
         {
             lnames <- strsplit(nams[i], split = '_', fixed = TRUE)[[1]]
             larg <- paste(lnames, collapse ='"]][["')
             expr <- paste0('theme[["',larg,'"]] <- args[[i]]')
             eval(str2lang(expr))
         }
    )
  }
  invisible(trellis.par.set(theme=theme))
}
#' @describeIn td trellis par set for all parameters, special treatment for col, lwd, lty, pch
#' @export
tps_ <- function (...)
{
  args <- list(...)
  # disp(args)
  theme <- trellis.par.get()
  nams <- names(args)
  for(i in seq_along(args)) {
    # disp(nams[i])
    switch(nams[i],
           'lty' = {
             theme$plot.line$lty <- args[[i]]
           },
           'col' = {
             theme$plot.line$col <- args[[i]]
             theme$plot.symbol$col <- args[[i]]
           }, 
           'symbol_col' = {
             
             theme$plot.symbol$col <- args[[i]]
             
           },
           'line_col' = {
             theme$plot.line$col <- args[[i]]
           },
           'lwd' = {
             theme$plot.line$lwd <- args[[i]]
             # disp('lwd')
           },
           'alpha' = {
             theme$plot.line$alpha <- args[[i]]
             theme$plot.symbol$alpha <- args[[i]]
             # disp('alpha')
           },
           'symbol_alpha' = {
             theme$plot.symbol$alpha <- args[[i]]
             # disp('did this too')
             
           },
           'line_alpha' = {
             theme$plot.line$alpha <- args[[i]]
             
           },
           'font' = {
             theme$plot.symbol$font <- args[[i]]
             
           },
           'fill' = {
             theme$plot.symbol$fill <- args[[i]]
             
           },
           'cex' = {
             theme$plot.symbol$cex <- args[[i]]
             
           },
           'pch' = {
             theme$plot.symbol$pch <- args[[i]]
           },
           {
             lnames <- strsplit(nams[i], split = '_', fixed = TRUE)[[1]]
             larg <- paste(lnames, collapse ='"]][["')
             expr <- paste0('theme[["',larg,'"]] <- args[[i]]')
             eval(str2lang(expr))
           }
    )
  }
  invisible(trellis.par.set(theme=theme))
}
### tpg and tpg_
#' @describeIn td trellis par set for superpose parameters
#' @export
tpg <- function (...)
{
  args <- list(...)
  # disp(args)
  theme <- trellis.par.get()
  if(length(args) == 0) return(theme)
  # nams <- names(args)
  ret <- list()
  for(i in seq_along(args)) {
    # disp(nams[i])
    ret[[args[[i]]]] <- switch(args[[i]],
           'lty' = {
             theme$superpose.line$lty
           },
           'col' = {
             if(identical(theme$superpose.line$col,theme$superpose.symbol$col))
               theme$superpose.line$col
             else
               list(superpose.line = theme$superpose.line$col, superpose.symbol = theme$superpose.symbol$col)
           }, 
           'symbol_col' = {
             
             theme$superpose.symbol$col
             
           },
           'line_col' = {
             theme$superpose.line$col
           },
           'lwd' = {
             theme$superpose.line$lwd
             # disp('lwd')
           },
           'alpha' = {
             list(theme$superpose.line$alpha,
             theme$superpose.symbol$alpha)
             if(identical(theme$superpose.line$alpha,theme$superpose.symbol$alpha))
               theme$superpose.line$alpha
             else
               list(superpose.line = theme$superpose.line$alpha, superpose.symbol = theme$superpose.symbol$alpha)
             
                         
           },
           'symbol_alpha' = {
             theme$superpose.symbol$alpha
         
             
           },
           'line_alpha' = {
             theme$superpose.line$alpha 
             
           },
           'font' = {
             theme$superpose.symbol$font 
             
           },
           'fill' = {
             theme$superpose.symbol$fill 
             
           },
           'cex' = {
             theme$superpose.symbol$cex
             
           },
           'pch' = {
             theme$superpose.symbol$pch
           },
           {
             lnames <- strsplit(args[[i]], split = '_', fixed = TRUE)[[1]]
             larg <- paste(lnames, collapse ='"]][["')
             expr <- paste0('theme[["',larg,'"]]')
             eval(str2lang(expr))
           }
    )
  }
  ret

}

###
###  xqplot
###
#' Extended Quantile Plots
#'
#' An easy way to see a dataset's variables at a glance. Shows uniform quantile
#' plot for numerical varibles and barcharts for factors. Quantile plots also
#' show a horizontal line at the position of the mean and at mean plus or minus
#' one standard deviation.
#'
#' @param x a data frame or list of variables to plot
#' @param ptype "quantile" (default) or "normal": kind of quantile to plot on x
#' axis.
#' @param labels names for each plot
#' @param \dots additional arguments passed to 'plot' command
#' @param mfrow number of rows and columns per page. If missing, an attempt is
#' made to choose a reasonable number.
#' @param ask if TRUE pause after each page, default: FALSE
#' @param mcex character expansion factor for marginal text
#' \code{mcex}
#' @param maxlab maximum number of categories to label in barcharts
#' @param debug if TRUE, print additional information
#' @param mar size of margins
#' @param text.cex.factor character expansion factor for barchart labels
#' @param left.labs determines placement of barchart labels
#' @param class show class of object, default: TRUE
#' @param xlab.pos position of xlab, default: 2
#' @param xlab.cex cex for xlab, default: .7
#' @param sublab.pos position for sublab, default: xlab.pos + xlab.cex * 1.15
#' @param xaxs style of x axis, default: 'i'
#' @param maxvarnamelength maximum length of variable name without splitting on
#' two lines.
#' @note Bugs:
#' 'mfrow' should take the total number of variables into account if they will
#' fill more than one page so the last page is close to being full.
#'
#' The current version of the function could be made much simpler and more
#' transparent. Some code is redundant.
#' @examples
#' require(car)
#' xqplot(Prestige)
#' xqplot(Prestige,"n") # normal quantiles
#' @export
xqplot <- function(x,
                   ptype = "quantile",
                   labels = dimnames(x)[[2]],
                   ...,
                   mfrow = findmfrow (ncol(x)),
                   # ask = prod(mfrow) <
                   # ncol(x) && dev.interactive(),
                   ask = FALSE,
                   mcex = 0.8,
                   maxlab = 12 ,
                   debug = F,
                   mar = c(4, 2.5, 2, 1),
                   xlab.pos = 2,
                   # new param
                   xlab.cex = .7,
                   # new param
                   sublab.pos = xlab.pos + xlab.cex * 1.15,
                   # new param
                   text.cex.factor = 1 ,
                   left.labs = F,
                   class = TRUE,
                   xaxs = 'i',
                   maxvarnamelength = 20)
{
  ## Adapted from myplot.data.frame for R by G. Monette, Oct. 25, 2004
  ##    maxlab is maximum number of labels
  # Turn matrices into variables:
  x <- as.data.frame(x)
  if (any (sapply(x, class) == 'matrix')) {
    zz <- list()
    for (ii in seq_along(x)) {
      if (is.matrix(x[[ii]])) {
        if (is.null (colnames(x[[ii]]))) {
          cnames <- paste(names(x)[ii], 1:ncol(x[[ii]]), sep = '.')
        } else {
          cnames <- paste(names(x)[ii], colnames(x[[ii]]), sep = '.')
        }
        for (jj in seq_len(ncol (x[[ii]]))) {
          zz[[cnames[jj]]] <- x[[ii]][, jj]
        }
        
      } else {
        zz[[names(x)[[ii]]]] <- x[[ii]]
      }
    }
    x <- as.data.frame(zz)
    #disp( x )
  }
  
  
  left.labs <- rep(left.labs, length = length(x))
  findmfrow <- function(x) {
    if (x > 9)
      c(3, 4)
    else
      cbind(
        '1' = c(1, 1),
        '2' = c(1, 2),
        '3' = c(2, 2),
        '4' = c(2, 2),
        '5' = c(2, 3),
        '6' = c(2, 3),
        '7' = c(3, 3),
        '8' = c(3, 3),
        '9' = c(3, 3)
      ) [, x]
  }
  
  opt <- par(mfrow = mfrow,
             ask = ask ,
             mar = mar + 0.1)
  on.exit(par(opt))
  if (debug) {
    cat("opt:\n")
    print(opt)
  }
  
  iscat <- function(x)
    is.factor(x) || is.character(x)
  
  Levels <- function(x) {
    if (is.factor(x))
      levels(x)
    else
      unique(x)
  }
  
  
  compute.cex <- function(x) {
    ll <- length(x)
    cex <- 2 * ifelse(ll < 5, 2,
                      ifelse(ll < 10, 1,
                             ifelse(ll < 20, .7, .5))) / mfrow[1]
  }
  for (ii in 1:dim(x)[2]) {
    vv <- x[[ii]]
    nam <- labels[[ii]]
    Nmiss <- sum(is.na(vv))
    N <- length(vv)
    if (iscat(vv)) {
      tt <- table(vv)
      
      xlab <- paste("N =", N)
      if (Nmiss > 0) {
        tt <- c("<NA>" = sum(is.na(vv)), tt)
        xlab <- paste(xlab, "  Nmiss =", Nmiss)
      }
      ll <- names(tt)
      nn <- length(ll)
      if (left.labs[ii]) {
        barplot(
          tt,
          horiz = TRUE,
          xlab = '',
          cex.names = text.cex.factor * compute.cex(nn)
        )
        mtext(xlab, 1, xlab.pos, cex = xlab.cex)
      }
      else {
        zm <- barplot(tt,
                      names = rep("", nn),
                      horiz = TRUE,
                      xlab = '')
        mtext(xlab, 1, xlab.pos, cex = xlab.cex)
        
        ## If nn > maxlab drop labels for smaller frequencies
        sel <- rep(T, length(tt))
        tt.sorted <- rev(sort(tt))
        if (nn > maxlab)
          sel <- tt > tt.sorted[maxlab]
        if (debug) {
          disp(sel)
          disp(nam)
          disp(tt)
          disp(tt.sorted)
          disp(maxlab)
          disp(tt.sorted[maxlab])
          disp(sel)
          disp(zm[sel])
          disp(rep(max(tt), nn)[sel])
          disp(ll[sel])
        }
        if (any(sel))
          text(rep(max(tt), nn)[sel]  ,
               zm[sel],
               ll[sel],
               adj = 1,
               cex = text.cex.factor * compute.cex(nn))
      }
    } # end of iscat(vv)
    else {
      sublab <- ""
      N <- length(vv)
      Ninfinite <- 0
      if (any(is.infinite (vv))) {
        n.pi <- sum(vv == Inf , na.rm = TRUE)
        n.ni <- sum(vv == -Inf, na.rm = TRUE)
        Ninfinite <- n.pi + n.ni
        vv <- vv[!is.infinite(vv)]
        sublab <- paste(sublab, "-Inf:", n.ni, "+Inf:", n.pi)
      }
      Nmiss <- 0
      if (any (is.na(vv))) {
        Nmiss <- sum(is.na(vv))
        vv  <- vv[!is.na(vv)]
        sublab <- paste(sublab, "NA:", Nmiss)
      }
      Nok <- N - Nmiss - Ninfinite
      if (pmatch(ptype, 'normal', nomatch = 0) == 1) {
        xxvar <- qnorm(ppoints(length(vv)))
        xlab <- paste("Normal quantile for", Nok, "obs.")
      }
      else {
        xxvar <- ppoints(length(vv))
        xlab <- paste("Fraction of", Nok, "obs.")
      }
      
      ## Plot continuous
      if (Nok == 0) {
        xxvar <- 1
        vv <- 1
        if (sublab == "") {
          plot(xxvar,
               vv,
               xlab = '',
               ylab = "",
               type = 'n')
          # mtext(xlab,1, xlab.pos)
        } else {
          plot(xxvar,
               vv,
               xlab = '',
               ylab = "",
               type = 'n')
          mtext(sublab, 1, xlab.pos, cex = xlab.cex)
          #mtext(xlab,1, xlab.pos)
        }
        text(1, 1, "NA")
      }
      else {
        if (sublab == "") {
          plot(xxvar,
               sort(vv),
               xlab = '',
               ylab = "Data",
               ...)
          mtext(xlab, 1, xlab.pos, cex = xlab.cex)
        } else {
          plot(xxvar,
               sort(vv),
               xlab = '',
               ylab = "Data",
               ...)
          mtext(xlab, 1, xlab.pos, cex = xlab.cex)
          mtext(sublab, 1, sublab.pos , cex = xlab.cex)
        }
        xm <- mean(vv)
        xs <- sqrt(var(vv))
        abline(h = xm, lty = 1)
        abline(h = c(xm - xs, xm + xs), lty = 2)
      }
    }
    ## Titles for all plots
    vlab <- labels[ii]
    line.offset <- 1.0
    if (nchar(vlab) > maxvarnamelength) {
      vlab <-
        paste(
          substring(vlab, 1, maxvarnamelength),
          "\n",
          substring(vlab, maxvarnamelength + 1)
        )
      line.offset <- 0.2
    }
    mtext(vlab, 3, line.offset , cex = mcex)
    if (class)
      mtext(paste(class(vv), collapse = ', '), 3, 0.2, cex = .7 * mcex)
  }
  # par(opt)
  if (debug) {
    disp(par())
  }
  invisible(0)
}

#' Show available characters, colours, etc.
#'
#' @param n FIXME
#' @param all default: FALSE
#' @concept lattice
#' @seealso \code{\link[lattice]{lattice::show.settings}}
#' @export
sampler <-
    function( n=24, all = FALSE ) {
    # sample of lines and symbols
     old.par <- par(ask=T)
     on.exit( par(old.par))
      require(lattice)

     y <- 0:n
     x <- 0:n
if(all) {
     print(xyplot( y ~ x, type = 'n', xlab = 'lty', ylab = 'col',
      panel = function(x,y,...) {
      for ( i in x) {
       panel.xyplot(c(i,i),range(y),type='l',lty=i,col=1,lwd = 3)
      }
      for ( i in y) {
       for ( j in seq(0,.9, by = .1)) {
        panel.xyplot(c(min(x)+ j*(max(x)-min(x)),min(x)+ (j+.1)*(max(x)-min(x))),c(i,i),type='l',lty=1,col=i, lwd = 3)
       }
      }
     }))

     # print(z$x, z$y, ylim=c(0,7))
     spl <-trellis.par.get('superpose.line')
     z <- expand.grid( y = 1:length(spl$lty), x = 0:2)
     print(xyplot( y ~ x , z, ylim =c(0,length(spl$lty)),groups = y, type='b',
            main="superpose.line and .symbol"))
}
     y <- 10*(0:25)
     x <- 0:9
     print(xyplot( y ~ x, type = 'n', main = 'pch',
        xlab = expression( ~ alpha + beta + gamma + delta[epsilon] + zeta^eta + theta + iota+kappa),
        ylab = expression( ~ lambda + mu + nu + xi + omicron + pi + rho + sigma + tau + upsilon + phi + chi +psi + omega),
      panel = function(x,y,...) {
      for ( i in x) {
       for ( j in y ) {
        panel.xyplot(i,j,pch=i+j,cex = 2)
       }
      }
     }))
     print(show.settings())

     invisible(0)
}

#' Generate a palette of colours -- possibly superseded
#'
#' @param col colors to show
#' @param border (default 'light gray')
#' @param \dots FIXME
#' @export
pal <- function(col=c('blue','pink'), border = "light gray", ...) {
     n <- length(col)
     plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
             xlab = "", ylab = "", ...)
     rect(0, 0:(n - 1)/n, .6, 1:n/n,  col = col, border = border)
     ret <- col2rgb(col)
     dimnames(ret)[[2]] <- col
     ret <- t(ret)
     txt <- paste( as.character(col), "(",
        apply( ret, 1, paste, collapse=" "), ")")
     text( rep(.6, n), (0:(n-1)+.5)/n, txt, adj = 0)
     ret <- col2rgb(col)
     dimnames(ret)[[2]] <- col
     t(ret)
}
#' Display selected colors n at a time
#'
#' @param pattern match to select colors
#' @param ignore match to ignore colors (default 'grey|gray') 
#' @param n maximum number to show in one page (default 30)
#' @param ask (default FALSE) prompt for new page
#' @export
pals <- function(pattern = '', ignore = '^gray[0-9]|^grey[0-9]', 
                 n = 30, ask = FALSE){
    if(is.numeric(pattern)) {
      n <- pattern
      pattern <- ''
    }
    cc <- grep(pattern, colors(), value = TRUE)
    if( !is.null(ignore) && !(ignore == '')) {
      cc <- grep(ignore, cc, value = TRUE, invert = TRUE)
    }
    N <- length(cc)
    ii <- 1
    while( ii < N ){
        pal(cc[ii:min(ii+n,N)], ask = ask)
        ii <- ii + n + 1
    }
}

## brace() now moved to brace.R

#' Replace elements of x with correspondingly named elements of ll
#'
#' @param x FIXME
#' @param ll FIXME
#' @export
change <- function(x,ll) {
  #
  # Modifies elements in a list
  # Ideal for changing ... arguments in calls to panel.groups etc.
  #
  nams <- names(ll)
  for ( ii in  seq_along(ll) ) {
    x[[nams[ii]]] <- ll[[ii]]
  }
  x
}
#' Replace elements of x with correspondingly named elements of ll
#'
#' @param x list or vector some of whose elements will be replaced or appended to
#' @param ... elements for replacement or appending
#' @export
change <- function(x,...) {
  #
  # Modifies elements in a list
  # Ideal for changing ... arguments in calls to panel.groups etc.
  #
  ll <- list(...)
  nams <- names(ll)
  if(is.null(nams)) nams <- rep('',length(ll))
  for ( ii in  seq_along(ll) ) {
    x[[nams[ii]]] <- ll[[ii]]
  }
  x
}

#' Panel function to display subgroups within groups within panels
#'
#' This function is designed to be used as the argument to \code{panel.groups}
#' in \code{xyplot}. It effectively adds another level of subgrouping to that
#' implemented by the \code{groups} argument in \code{xyplot}.  Useful mainly
#' to display data and fitted lines in groups within panels.
#'
#' This function is designed to be used as the argument to 'panel.groups' in
#' 'xyplot'. It allows the plotting of points versus lines within subgroups of
#' the data identified by the 'groups' argument in xyplot.  It requires a
#' variable to identify the subgroups. Points or lines are used within
#' subgroups depending on 'subgroups.type' where the order is that of the
#' levels of the 'subgroups' argument coerced as a factor, if necessary.  See
#' the examples below.
#'
#' @param x,y coordinates of the points to be displayed
#' @param subscripts subscripts giving indices in original data frame
#' @param subgroups a subgrouping variable. Use a full reference, e.g.
#' data$subvar
#' @param subgroups.type plotting type, typically 'p' or 'l', for each level of
#' the variable passed through the \code{subgroups} argument
#' @param type FIXME
#' @param panel.subgroups function use to plot data within in each group
#' referring to the levels of the variable passed by \code{subgroups}.  Define
#' a \code{panel.subgroups} argument in the call to \code{xyplot} and it will
#' be used to plot the groups. See the examples below.
#' @param \dots any other arguments to be passed to the panel plotting function
#' @seealso \code{link[lattice]{panel.superpose}},
#' \code{link[lattice]{panel.superpose.2}}, \code{link[lattice]{panel.xyplot}}
#' @examples
#' \dontrun{
#' library(car)
#' data(Prestige)
#' fit <- lm( prestige ~ (education +I(education^2)) * type, Prestige, na.action = na.omit)
#' pred <- expand.grid( education = seq( 6, 18, .5), type = levels( Prestige$type))
#' pred$prestige <- predict( fit, newdata = pred )
#'
#' Prestige$what <- 'data'
#' pred$what <- 'fit'         # this works because 'fit' follows 'data' lexicographically
#'
#' combined <- merge( Prestige, pred, all = T)
#'
#' xyplot( prestige ~ education, combined,
#'           groups = type,
#'           subgroups = combined$what,  # note that a full reference to the variable is needed
#'           panel = panel.superpose,    # might not be necessary in future version of lattice
#'           panel.groups = panel.subgroups)  # uses the default of points for the first level of 'what'
#'                                            # and lines for the second level
#'
#' ## Using the argument 'panel.subgroups' instead of the default 'panel.xyplot'
#' ## Note that panel.subgroups is a function (this one) and also an argument that
#' ## is a function passed to the function. The argument defines the action to
#' ## be taken within each level of 'what'
#'
#' xyplot( prestige ~ education, combined,
#'         groups = type,
#'         subgroups = combined$what,  # note that a full reference to the variable is needed
#'         panel = panel.superpose,    # might not be necessary in future version of lattice
#'         panel.groups = panel.subgroups,
#'         panel.subgroups = function( x, y, subgroup, type, ...) {
#'              # note that you need to include 'type' among the arguments
#'              # if you need to prevent it from being passed through '...'
#'              # When called, this function will be passed arguments
#'              # subgroup, subgroup.number, subscripts, and type from
#'              # subgroups.type.
#'        if ( subgroup == 'data' ) {
#'               panel.xyplot( x, y, ...)
#'               panel.lines( dell(x,y), ...)
#'        } else {
#'               panel.lines( x,y, ...)
#'        }
#'    })
#' }
#' @export
panel.subgroups <- function( x, y, subscripts,
     subgroups, subgroups.type = c('p','l'),type,
     panel.subgroups = panel.xyplot, ...) {
help = "Use help: ?panel.subgroups"
         subgroups <- as.factor(subgroups)
         levs <- levels(subgroups)
         subgroups.type <- rep( subgroups.type, length.out = length(levs))

         subgroups = subgroups[subscripts]
         for ( i in seq_along( levs) ) {
             sel <- subgroups == levs[i]
             if ( any( sel )) {
                panel.subgroups( x[sel], y[sel], type = subgroups.type[i],
                     subscripts = subscripts[sel], subgroup.number = i,
                     subgroup = levs[i], ...)
             }
         }
}
#' Display RColorBrewer palette
#' 
#' Function name that follows usual conventions and is
#' easier to remember.
#' 
#' @export
brewer.pal.show <- RColorBrewer::display.brewer.all