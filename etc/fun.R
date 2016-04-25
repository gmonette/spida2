# "
# LOG OF CHANGES:
# 2014:
#    Sep 10: Added gd from fun.R
#    Sep 9: Modified 'tolong' to handle ragged wide files (i.e. not the same
#           times for each varying variable) and 'reshape' bug that allocates
#           times incorrectly if all varying variables do not have times in the
#           same order.
# 2013:
#    Aug 15: added ability to other arguments via trellis.par.set.
# 2011:
#    Jun 16: added osculant and osculant.default to generate locus
#            of osculation between two families of ellipses
#    Mar 16: added capply.formula based on aggregate.formula
# 2010:
#    Aug 13: removed obsolete version of Lmat in fun.R superceded by version in
#            wald.R
#            Ldiff moved to wald.R
#            TODO: clean up other wald and eventually spline stuff
#    Jul 21: renames Lag and LagI, cLag and cLagI respectively avoid conflicts
#            with Hmisc::Lag.  Lag and and LagI are kept as aliases.
#    Jun 16: Added non-ordered version of reorder.factor
#    Jan 30: Added cvars
#
# 2009:
#    May 8:  Added naresid to model.matrix.lme
#
#
#
# '''fun.R:  A collection of utility functions and functions for multilevel modeling'''
# '''        kept by Georges Monette.  These file can be downloaded or sourced in R'''
# '''        from http://www.math.yorku.ca/~georges/R/fun.R'''
#
# '''NOTE 1: Please send any suggestions or problems to georges@yorku.ca.
#
# '''NOTE 2: There is a copy of this file at //wiki.math.yorku.ca/index.php/R:fun.R'''
# '''        which can be edited. Changes will be incorporated regularly into the'''
# '''        the downloadable file.'''
#
# '''NOTE 3: THIS FILE USES QUOTES AND TAGS SO IT CAN BE RENDERED AS A WIKI FILE AS WELL AS R SOURCE'''
#
# '''BE SURE TO LEAVE WIKI TEXT IN DOUBLE QUOTES AND AVOID DOUBLE QUOTES IN WIKI TEXT'''
#
# Last uploaded to http://www.math.yorku.ca/~georges/R/fun.R : Auguest 23, 2006
#
# :Decribe modifications here
# :March 2009
# :: panel.subgroups  allows different plotting 'types' within subgroups
# :July 13, 2008
# :: gsp: general spline program
# :: smsp: smoothing spline using random effects model
# :New: Aug 10, 2007
# :: pchisq.mix: mixture of chisquares for different dfs for testing hypotheses
# :::  concering null random effects in lme. Use simulate.lme to verify whether correct
# :New: May 27, 2007
# :: sasin      - read a SAS ODS CSV file and extract individual tables into a list
# :: brace
# :New: April 9, 2007
# :: vif.lme    - variance inflation factors for lme fits
# :: Rbind      - combine data and prediction data frame to plot
# :::          results together
# :New: November 9, 2006
# :: oplot - plots number of observations that overplot
# :New: August 23, 2006
# :: cvar( x, id ) ; dvar( x, id )
# ::: cvar creates a contextual variable that is the group mean of 'x' within each level of the factor 'id'. If 'x' is a factor, 'cvar' returns a suitably labelled matrix that is the group mean of the coding variables for the factor 'x'.
# ::: dvar is x - cvar( x, id ) where x is turned into its coding matrix if x is a factor.
# :New: July 10, 2006
# :: xmerge( x , y , by )
# ::: merge that tests consistency of common names not in the by argument. Values are combined with priority given to 'y'. If variables are inconsistent, the '.x' and '.y' versions are also left in the merged data frame.
# ::: merge two data.frames with diagnostics
# :: long(data, varying)
# ::: reshapes data from wide to long format using variable names given in list 'varying'. Vectors in the list are old names, names of vectors are new names. Names alone serve are roots.
# :Modified: June 10, 2006
# :: added:
# ::: up(dd, form)  - creates a higher level data set with one row per case defined by 'form'
# :Modified:  --[[User:Georges|Georges Monette]] 14:30, 28 May 2006 (EST)
# :: added functions from funRDC
# :Modified:  --[[User:Georges|Georges Monette]] 05:56, 22 Feb 2006 (EST)
# :: added Lmat
# :Modified:  --[[User:Georges|Georges Monette]] 09:45, 2 Nov 2005 (EST)
# :: added class 'cat' to coursefun
# :: defined print.cat to print with cat
# :Modified:
# :: plot3d, identify3d by John Fox so they work with matrices, data.frames or variable arguments
# :: ell  - corrected error when using radius
# ::
#
#
#
# == General description ==
# <pre>
# "
# ##
# ##
# ##  Some R functions for PSYC 6140 and MATH 6630
# ##  2005-2006
# ##
# ## Last update: October 27, 2005
# ## Summary:
# ##
# ##
# ## Tables:
# ##
# ##    atotal: border an array with sums
# ##    abind : glue two arrays together on a selected dimension
# ##
# ## General Linear Hypothesis
# ##
# #     glh   : glh( fit, L )
# ##    Lmat  : generates L matrix for anova in library(lmer) or lht in library(car)
# ##
# ## Graphics:
# ##
# ##    td    : easy front end to trellis.device and related functions
# ##    xqplot: extended quantile plots
# ##
# ## 3D graphics by John Fox:
# ##
# ##    scatter3d
# ##    identify3d
# ##    ellipsoid
# ##    plot3d     - wrapper for scatter3d
# ##
# ## Inference
# ##    cell   - a modified version of car::confidence.ellipse.lm that
# ##             creates a confidence ellipse to plot with lines(...)
# ##    dell   - data ellipse
# ##
# ##
# # Note that some functions that were transferred and improved in coursefun.R
# # have been 'disabled' by appending '.rdc' to function name
# #
# #  Splus functions written in RDC
# #  2005:
# #  April 19       Modified for R
# #  May 3          copied atotal, abind
# #                 wrote acond
# #  May 10         new function: cap1 to capitalize initial letters and turn underscores to blanks
# #  May 12         new function adapted from gm: td
# #  May 13         new function: write.sas to write data frame to SAS without truncating variable names
# #  May 16         added Poly and centering to splines
# #  June 13        added constant to check if values are constant within levels of id
# #                 added varLevel( data.frame, ~lev1/lev2) to report level of each variable
# #                 added Lmat to generate L matrix with 1's for effects containing a string
# #                 modified anova.lme to use min DFs in 'denominator' DFs
# #  August 3       modified Lag to accept 'at' argument
# #  August 5       changed 'arep' to 'apct' in order to parallel atotal and acond
# #                 changed acond to aprop
# #  August 15      getFix, glh and and print.glh, Q, Vcov, Vcor
# #  August 25      Contrasts
# #  2006
# #  May 19         fill and capply
# #  October 2      cvar:  create contextual variable
#
#
# ##
# ##  Crude predict methods for 'mer'  Dec 6, 2008
# ##
#
#
#
# #' A tentative version of predict for mer objects
# #'
# #'
# #'
# #'
# #'
# #' @param model
# #' @param data
# #' @param form
# #' @param verbose
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (model, data = model.matrix(model), form, verbose = FALSE)
# #' {
# #'     help = "\n    This is a crude predict method for class 'mer'\n    Limitations:\n    1) When invoked on the model it only returns the linear predictor for the\n       complete data. Note that a complete data frame is easily obtained with\n       model.frame( model )\n    2) When the 'data' argument is provided for 'new' data, it is also\n       necessary to provide the correct rhs of the fixed effects formula.\n\n    "
# #'     if (missing(form) & !missing(data)) {
# #'         cat(help)
# #'         stop("Need 'form' if 'data' given")
# #'     }
# #'     if (!missing(data)) {
# #'         data = model.matrix(form, data)
# #'         cnames = colnames(data)
# #'         if (verbose)
# #'             print(cnames)
# #'         fnames = names(fixef(model))
# #'         if (verbose)
# #'             print(fnames)
# #'         if (any(cnames != fnames)) {
# #'             cat("\nMatrix names:\n")
# #'             print(cnames)
# #'             cat("\nCoeff names:\n")
# #'             print(fnames)
# #'             warning("matrix and coeff names not the same")
# #'         }
# #'     }
# #'     data %*% fixef(model)
# #'   }
# #'
# #' @export
# predict.mer <- function( model, data = model.matrix(model), form , verbose = FALSE) {
#
# help   = "
#     This is a crude predict method for class 'mer'
#     Limitations:
#     1) When invoked on the model it only returns the linear predictor for the
#        complete data. Note that a complete data frame is easily obtained with
#        model.frame( model )
#     2) When the 'data' argument is provided for 'new' data, it is also
#        necessary to provide the correct rhs of the fixed effects formula.
#
#     "
#         if (missing(form) & ! missing(data)) {
#             cat( help)
#             stop( "Need 'form' if 'data' given")
#         }
#         if( ! missing( data )){
#             data = model.matrix(form, data)
#             cnames = colnames(data)
#             if( verbose ) print( cnames )
#             fnames = names( fixef( model ))
#             if (verbose) print( fnames)
#             if ( any( cnames != fnames)) {
#                 cat("\nMatrix names:\n")
#                 print( cnames )
#                 cat("\nCoeff names:\n")
#                 print( fnames)
#                 warning("matrix and coeff names not the same")
#             }
#         }
#         data %*% fixef( model)
#
# }
#
#
#
#
# ##
# ## Linear algebra
# ##
#
#
#
#
#
#

# #' Indicate number of points overplotted
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param y
# #' @param \dots
# #' @param verbose
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, y, ..., verbose = T)
# #' {
# #'     pat <- paste(x, y, sep = ",")
# #'     keep <- !duplicated(pat)
# #'     ns <- table(pat)
# #'     ns <- ns[pat[keep]]
# #'     nps <- as.character(ns)
# #'     x <- x[keep]
# #'     y <- y[keep]
# #'     if (verbose) {
# #'         print(pat)
# #'         print(table(pat))
# #'         print(keep)
# #'         print(ns)
# #'     }
# #'     plot(x, y, pch = "o", cex = 5, ...)
# #'     text(x, y, nps)
# #'   }
# #'
# #' @export
# oplot <- function( x, y, ..., verbose = TRUE) {
#     pat <- paste( x, y, sep = ",")
#     keep <- !duplicated(pat)
#
#     ns <- table(pat)
#     ns <- ns[pat[keep]]      # to order ns so it matches pat[keep]
#     nps <- as.character(ns)
#     #nps [ns>9] <- "*"
#     x <- x[keep]
#     y <- y[keep]
#
#
#     if (verbose) {
#         print(pat)
#         print(table(pat))
#         print(keep)
#         print(ns)
#     }
#     plot( x, y,pch = "o",cex = 5,...)
#     text( x, y, nps)
#
# }
# #plot( c(1,2,3,2,1,2,3), c(1,2,3,2,1,2,3), pch = 'AB')
# #oplot( c(rep(1,10),2,3,2,1,2,3), c(rep(1,10),2,3,2,1,2,3), type = 'b')
#
#
#
#
# #' Left square root of X'X
# #'
# #'
# #'
# #'
# #'
# #' @param x
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x)
# #' {
# #'     xx <- svd(x)
# #'     ret <- t(xx$v) * sqrt(pmax(xx$d, 0))
# #'     ret
# #'   }
# #'
# #' @export
# fac <- function(x) {
#    xx <- svd(x)
#    ret <- t(xx$v) * sqrt(pmax( xx$d,0))
#    ret  #ret [ nrow(ret):1,]
# }
#
#
#
#
#
# #' Display the name and value of an object
# #'
# #' Display the name and value of an object, useful for debugging
# #' concise (1-5 lines) description of what the function does. ~~
# #'
# #'
# #'
# #' @param x
# #' @param head
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, head = deparse(substitute(x)))
# #' {
# #'     cat("\n::: ", head, " :::\n")
# #'     print(x)
# #'   }
# #'
# #' @export
# disp <- function( x , head = deparse(substitute(x))) {
#     # for debugging
#     cat("\n::: ", head , " :::\n")
#     print(x)
# }
#
#
# "
# </pre>
#
# == General Linear Hypothesis ==
# <pre>
# "
#
#
#
#
# #' Attempt at streamlining expand.grid for prediction
# #'
# #'
# #'
# #'
# #'
# #' @param df
# #' @param by
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (df, by, ...)
# #' {
# #'     dots = list(...)
# #'     by = model.frame(by, df, na.action = na.include)
# #'     byn = names(by)
# #'     names(byn) = byn
# #'     args = lapply(byn, function(x) {
# #'         vv = df[[x]]
# #'         if (is.factor(vv))
# #'             levels(vv)
# #'         else unique(vv)
# #'     })
# #'     args = c(args, dots)
# #'     do.call("expand.grid", args)
# #'   }
# #'
# eg <-
# function( df, by, ...) {
#     # a quicker version of expand.grid
#     # should work with fits
#         dots = list(...)
#         by = model.frame( by, df, na.action = na.include)
#         byn = names(by)  # will this stay
#         names(byn) = byn
#         args = lapply( byn, function(x){
#                 vv = df[[x]]
#                 if ( is.factor(vv) ) levels(vv) else unique(vv)
#         })
#         args = c( args,dots)
#         do.call('expand.grid', args)
# }
#
#
#
#
#
# #' Transform NAs to FALSE
# #'
# #'
# #'
# #'
# #'
# #' @param x
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x)
# #' {
# #'     x[is.na(x)] <- F
# #'     x
# #'   }
# #'
# #' @export
# na2f <- function(x)  {
#     x[is.na(x)] <- FALSE
#     x
# }
#
#
# # fit <- lmer( Yield ~ Location *   Family  + (1|Block), data = Genetics)
# # getFix(fit)
# # fit <- lme( Yield ~ Location *   Family  , data = Genetics, random = ~1|Block)
#
# # L <- rbind( c(1,1,1,1), c(0,1,0,0), c(1,0,1,1))
#
#
#
#
# #' Utility function
# #'
# #'
# #'
# #'
# #'
# #' @param x
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x)
# #' {
# #'     ret <- NULL
# #'     if (is.character(x))
# #'         ret <- "it's character"
# #'     else {
# #'         if (is.numeric(x)) {
# #'             if (is.null(dim(x))) {
# #'                 if (length(x) != 4)
# #'                   ret <- diag(4)[x, ]
# #'                 else ret <- rbind(x)
# #'             }
# #'         }
# #'     }
# #'     ret
# #'   }
# #'
# #' @export
# tfun <- function( x) {
#       ret <- NULL
#       if ( is.character(x)) ret <- "it's character"
#       else {
#          if ( is.numeric(x)) {
#             if ( is.null(dim(x))) {
#                 if ( length(x) != 4 ) ret <- diag(4)[x,]
#                 else ret <- rbind( x )
#             }
#          }
#       }
#       ret
# }
#
# ##
# ##  Extension of avp from car
# ##
#
#
#
#
#
# #' Create data frame for added variable plot
# #'
# #' av.frame( model, variable) returns a data frame with model.frame(model)
# #' augmented by y.res and x.res, the residuals for an added variable plot
# #'
# #' The purpose of this function is to facilitate OLS av.plots for mixed models.
# #'
# #'
# #'
# #' @param model
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'        library(nlme)
# #'        library(lattice)
# #'        hs <- read.csv( 'http://www.math.yorku.ca/~georges/Data/hs.csv')
# #'
# #'        # Mixed model where ses and Sex are Level 1 and Sector is Level 2
# #'
# #'        fit.mm <- lme( mathach ~ ses * Sex * Sector, hs, random = ~ 1+ses| school)
# #'
# #'        # for diagnostics fit an OLS model using only level 1 variables interacting
# #'        # with the id variable
# #'
# #'        fit.ols <- lm( mathach ~ (ses * Sex ) * factor(school), hs)
# #'        xyplot( y.res ~ x.res | factor(school), cbind(av.frame(fit.ols, 'ses:Sex'),hs), sub = 'ses:Sex')
# #'        xyplot( y.res ~ x.res | factor(school), cbind(av.frame(fit.ols, '^Sex'),hs), sub = 'Sex')
# #'        xyplot( y.res ~ x.res | factor(school), cbind(av.frame(fit.ols, '^ses$|^ses:f'),hs), sub = 'ses')
# #'
# #'
# #'        Note : y.res is the residual from fitting the response on
# #'               the model matrix for fit.ols omitting any column
# #'               whose names is matched (as a regular expression)
# #'               by 'effect'
# #'               x.res is the residual of the first column of the
# #'               model matrix that is matched by 'effect' on the
# #'               same matrix used for y.res.
# #'        Caution: To make sure that the correct columns were
# #'               matched, the list of matched columns that are omitted
# #'               is printed.
# #'
# #' ## The function is currently defined as
# #' function (model, ...)
# #' {
# #'     UseMethod("av.frame")
# #'   }
# #'
# #' @export
# av.frame <- function( model, ..., help = FALSE) {
# if(help) {
#  cat("
#        av.frame( model, variable)
#        returns a data frame with model.frame(model) augmented
#        by y.res and x.res, the residuals for an added variable
#        plot
#
#        The purpose of this function is to facilitate OLS av.plots
#        for mixed models.
#
#
#
#        Example:
#
#        library(nlme)
#        library(lattice)
#        hs <- read.csv( 'http://www.math.yorku.ca/~georges/Data/hs.csv')
#
#        # Mixed model where ses and Sex are Level 1 and Sector is Level 2
#
#        fit.mm <- lme( mathach ~ ses * Sex * Sector, hs, random = ~ 1+ses| school)
#
#        # for diagnostics fit an OLS model using only level 1 variables interacting
#        # with the id variable
#
#        fit.ols <- lm( mathach ~ (ses * Sex ) * factor(school), hs)
#        xyplot( y.res ~ x.res | factor(school), cbind(av.frame(fit.ols, 'ses:Sex'),hs), sub = 'ses:Sex')
#        xyplot( y.res ~ x.res | factor(school), cbind(av.frame(fit.ols, '^Sex'),hs), sub = 'Sex')
#        xyplot( y.res ~ x.res | factor(school), cbind(av.frame(fit.ols, '^ses$|^ses:f'),hs), sub = 'ses')
#
#
#        Note : y.res is the residual from fitting the response on
#               the model matrix for fit.ols omitting any column
#               whose names is matched (as a regular expression)
#               by 'effect'
#               x.res is the residual of the first column of the
#               model matrix that is matched by 'effect' on the
#               same matrix used for y.res.
#        Caution: To make sure that the correct columns were
#               matched, the list of matched columns that are omitted
#               is printed.
#
# ")
#   return( invisible(0))
# }
#          UseMethod("av.frame")
# }
#
#
#
#
#
# #' lm method for av.frame
# #'
# #'
# #'
# #'
# #'
# #' @param model
# #' @param variable
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (model, variable, ...)
# #' {
# #'     mod.mat <- model.matrix(model)
# #'     var.names <- colnames(mod.mat)
# #'     omit <- grep(variable, var.names)
# #'     if (0 == length(omit))
# #'         stop(paste(variable, "is not matched among columns of the model matrix."))
# #'     cat("x.var =", var.names[omit[1]], "\n", "omitted vars =",
# #'         var.names[omit[-1]], "\n")
# #'     response <- response(model)
# #'     x.var <- mod.mat[, omit[1]]
# #'     Xpred <- mod.mat[, -omit]
# #'     preds <- predict(update(model, na.action = na.exclude))
# #'     responseName <- responseName(model)
# #'     if (is.null(weights(model)))
# #'         wt <- rep(1, length(response))
# #'     else wt <- weights(model)
# #'     res <- lsfit(mod.mat[, -omit], cbind(mod.mat[, omit[1]],
# #'         response), wt = wt, intercept = FALSE)$residuals
# #'     ret <- matrix(NA, nrow = length(preds), ncol = 2)
# #'     ret[!is.na(preds), ] <- res
# #'     data.frame(x.res = ret[, 1], y.res = ret[, 2])
# #'   }
# #'
# #' @export
# av.frame.lm <- function (model, variable,...){
# # code borrowed from 'car' by J. Fox, function 'av.plot'
# # labels = names(residuals(model)[!is.na(residuals(model))]),
# #    identify.points = TRUE, las = par("las"), col = palette()[2],
# #    pch = 1, lwd = 2, main = "Added-Variable Plot", ...)
#
#
#     mod.mat <- model.matrix(model)
#     var.names <- colnames(mod.mat)
#     omit <- grep( variable, var.names)
#     if (0 == length(omit))
#         stop(paste(variable, "is not matched among columns of the model matrix."))
#
#     cat( "x.var =", var.names[ omit[1] ], "\n",
#     "omitted vars =", var.names[omit[-1]], "\n")
#
#     response <- response(model)
#     x.var <- mod.mat[,omit[1]]
#     Xpred <- mod.mat[, - omit ]
#     preds <- predict( update( model, na.action = na.exclude))
#
#     responseName <- responseName(model)
#     if (is.null(weights(model)))
#         wt <- rep(1, length(response))
#     else wt <- weights(model)
#     res <- lsfit(mod.mat[, -omit], cbind(mod.mat[, omit[1]], response),
#         wt = wt, intercept = FALSE)$residuals
#     ret <- matrix(NA, nrow = length( preds), ncol = 2)
#     ret[ !is.na(preds),] <- res
#     data.frame( x.res = ret[,1], y.res = ret[,2])
#
# }
#
#
#
#
#
# #' Variance Inflation Factors for Mixed Models
# #'
# #' Calculates versions of the variance-inflation and generalized
# #' variance-inflation factors for mixed models.
# #'
# #' The concept of Variance Inflation in linear models can be applied to mixed
# #' models in a number of ways since the variance-covariance matrix of the
# #' estimated fixed coefficients is not simply proportional to the inverse of
# #' the cross-product matrix for the data. This method for the generic function
# #' \code{vif} in the \code{car} package, implements the version based on the
# #' variance-covariance funtion of the estimated fixed coefficients.  Since the
# #' \code{vcov} method is available for \code{lm} objects and \code{lme}
# #' objects, uses the code for the \code{lm} method for \code{lme} objects.
# #'
# #' @param mod
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (mod)
# #' {
# #'     if (any(is.na(fixef(mod))))
# #'         stop("there are aliased coefficients in the model")
# #'     v <- vcov(mod)
# #'     mm <- model.matrix(formula(mod), mod$data)
# #'     assign <- attributes(mm)$assign
# #'     if (names(fixef(mod)[1]) == "(Intercept)") {
# #'         v <- v[-1, -1]
# #'         assign <- assign[-1]
# #'     }
# #'     else warning("No intercept: vifs may not be sensible.")
# #'     terms <- labels(terms(mod))
# #'     n.terms <- length(terms)
# #'     if (n.terms < 2)
# #'         stop("model contains fewer than 2 terms")
# #'     R <- cov2cor(v)
# #'     detR <- det(R)
# #'     result <- matrix(0, n.terms, 3)
# #'     rownames(result) <- terms
# #'     colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
# #'     for (term in 1:n.terms) {
# #'         subs <- which(assign == term)
# #'         result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs,
# #'             -subs]))/detR
# #'         result[term, 2] <- length(subs)
# #'     }
# #'     if (all(result[, 2] == 1))
# #'         result <- result[, 1]
# #'     else result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
# #'     result
# #'   }
# #'
# "vif.lme" <-
# function (mod)
# {
#     if (any(is.na(fixef(mod))))
#         stop("there are aliased coefficients in the model")
#     v <- vcov(mod)  # vcov.lme is in library(stats)
#     mm <- model.matrix( formula(mod), mod$data)
#     assign <- attributes(mm)$assign
#     if (names(fixef(mod)[1]) == "(Intercept)") {
#         v <- v[-1, -1]
#         assign <- assign[-1]
#     }
#     else warning("No intercept: vifs may not be sensible.")
#     terms <- labels(terms(mod))
#     n.terms <- length(terms)
#     if (n.terms < 2)
#         stop("model contains fewer than 2 terms")
#     R <- cov2cor(v)
#     detR <- det(R)
#     result <- matrix(0, n.terms, 3)
#     rownames(result) <- terms
#     colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
#     for (term in 1:n.terms) {
#         subs <- which(assign == term)
#         result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs,
#             -subs]))/detR
#         result[term, 2] <- length(subs)
#     }
#     if (all(result[, 2] == 1))
#         result <- result[, 1]
#     else result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
#     result
# }
#
#
#
#
#
#
#
# "</pre>
#
# == Trellis graphics ==
# <pre>
# "
#
#
#

#
#
# "</pre>
# == Contingency tables ==
# <pre>"
#
#
#
#

#
# ###
# ###   tab
# ###
#
#
#
#
# #' otab
# #'
# #'
# #'
# #'
# #'
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (...)
# #' {
# #'     aa <- list(...)
# #'     if (length(aa) == 1 && is.list(aa[[1]])) {
# #'         return(do.call("tab", aa[[1]]))
# #'     }
# #'     for (ii in 1:length(aa)) aa[[ii]] <- factor(aa[[ii]], exclude = NULL)
# #'     ret <- do.call("table", aa)
# #'     ret
# #'   }
# #'
# #' @export
# otab <- function(...) {
# help <- "
# abind                fun.R     for PSYC 6140/MATH 6630 05/06
#
# Cross Tabulation and Table Creation Including Missing Values
#
# Description:
#
#      'tab' does the same thing as 'table' except that it includes
#      missing values for factors.  The argument 'exclude = NULL' to 'table'
#      results in the inclusion of missing values for numeric variable but
#      excludes missing values for factors. 'tab' is intended to remedy this
#      deficiency of 'table'.
#
# Usage:
#
#      tab(...)
#
# Arguments:
#      ...: objects which can be interpreted as factors (including
#           character strings), or a list (or data frame) whose
#           components can be so interpreted.
#
# Details:
#
# Value:
#
# Notes:
#
# References:
#
# Bugs:
#       Does not use argument name as a dimension name, in contrast with 'table'.
#
# Contributed by:  G. Monette  2005-10-10
#
# Modifications:
#
# "
#   args <- list(...)
#   if( is.list(args[[1]])) args <- args[[1]]
#   # for ( ii in 1:length(a)) if ( is.factor( a[[ii]])) a[[ii]] <- factor(a[[ii]],exclude = NULL)
#   for ( ii in 1:length(args))  args[[ii]] <- factor(args[[ii]],exclude = NULL)
#   do.call("table", args)
# }
#
# "</pre>
#
# == Ellipses: data and confidence ==
# <pre>"
#
# #' @export
# ellplus <- function ( center = rep(0,2), shape = diag(2), radius = 1, n = 100,
#                angles = (0:n)*2*pi/n,
#                fac = chol ,
#                ellipse = all,
#                diameters = all,
#                box = all,
#                all = FALSE) {
#         help <- "
#         ellplus can produce, in addition to the points of an ellipse, the
#         conjugate axes corresponding to a chol or other decomposition
#         and the surrounding parallelogram.
#         "
#         rbindna <- function(x,...) {
#             if ( nargs() == 0) return(NULL)
#             if ( nargs() == 1) return(x)
#             rbind( x, NA, rbindna(...))
#         }
#         if( missing(ellipse) && missing(diameters) && missing(box)) all <- TRUE
#         circle <- function( angle) cbind( cos(angle), sin(angle))
#         Tr <- fac(shape)
#         ret <- list (
#             t( c(center) + t( radius * circle( angles) %*% Tr)),
#             t( c(center) + t( radius * circle( c(0,pi)) %*% Tr)),
#             t( c(center) + t( radius * circle( c(pi/2,3*pi/2)) %*% Tr)),
#             t( c(center) + t( radius * rbind( c(1,1), c(-1,1),c(-1,-1), c(1,-1),c(1,1)) %*% Tr)))
#         do.call( 'rbindna', ret[c(ellipse, diameters, diameters, box)])
#     }
#
#
# #' @export
# dellplus <- function( x, y,  ...) {
#     if ( (is.matrix(x) && (ncol(x) > 1))|| is.data.frame(x)) mat <- as.matrix(x[,1:2])
#     else if (is.list(x)) mat <- cbind(x$x, x$y)
#     else mat <- cbind( x,y)
#     ellplus( apply(mat,2,mean), var(mat), ...)
#
# }
#
# # Replaced with version below from p3d
# # ell <- function(center = rep(0,2) , shape = diag(2) , radius = 1, n = 100,
# #         angles = (0:n)*2*pi/n) {
# #        circle <- radius * cbind( cos(angles), sin(angles))
# #        t( c(center) + t( circle %*% fac(shape)))
# # }
#
#
#
#
# #' Old confidence ellipse
# #'
# #'
# #'
# #'
# #'
# #' @param model
# #' @param which.coef
# #' @param levels
# #' @param Scheffe
# #' @param dfn
# #' @param center.pch
# #' @param center.cex
# #' @param segments
# #' @param xlab
# #' @param ylab
# #' @param las
# #' @param col
# #' @param lwd
# #' @param lty
# #' @param add
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (model, which.coef, levels = 0.95, Scheffe = FALSE,
# #'     dfn = 2, center.pch = 19, center.cex = 1.5, segments = 51,
# #'     xlab, ylab, las = par("las"), col = palette()[2], lwd = 2,
# #'     lty = 1, add = FALSE, ...)
# #' {
# #'     help <- "\nSee help for car::confidence.ellipse.lm\nexcept that 'cell' returns the points to form the ellipse\nwhich must be plotted with plot(...,type='l') or lines(...)\n-- Use dfn to determine Sheffe dimension, i.e. dfn = 1 to generate ordinary CIs, dfn = 2 for 2-dim CE, etc.\n"
# #'     require(car)
# #'     which.coef <- if (length(coefficients(model)) == 2)
# #'         c(1, 2)
# #'     else {
# #'         if (missing(which.coef)) {
# #'             if (has.intercept(model))
# #'                 c(2, 3)
# #'             else c(1, 2)
# #'         }
# #'         else which.coef
# #'     }
# #'     coef <- coefficients(model)[which.coef]
# #'     xlab <- if (missing(xlab))
# #'         paste(names(coef)[1], "coefficient")
# #'     ylab <- if (missing(ylab))
# #'         paste(names(coef)[2], "coefficient")
# #'     if (missing(dfn)) {
# #'         if (Scheffe)
# #'             dfn <- sum(df.terms(model))
# #'         else 2
# #'     }
# #'     dfd <- df.residual(model)
# #'     shape <- vcov(model)[which.coef, which.coef]
# #'     ret <- numeric(0)
# #'     for (level in rev(sort(levels))) {
# #'         radius <- sqrt(dfn * qf(level, dfn, dfd))
# #'         ret <- rbind(ret, c(NA, NA), ell(coef, shape, radius))
# #'     }
# #'     colnames(ret) <- c(xlab, ylab)
# #'     ret
# #'   }
# #'
# old.cell <-
# function (model, which.coef, levels = 0.95, Scheffe = FALSE, dfn = 2,
#     center.pch = 19, center.cex = 1.5, segments = 51, xlab, ylab,
#     las = par("las"), col = palette()[2], lwd = 2, lty = 1,
#     add = FALSE, ...)
# {
# help <- "
# See help for car::confidence.ellipse.lm
# except that 'cell' returns the points to form the ellipse
# which must be plotted with plot(...,type='l') or lines(...)
# -- Use dfn to determine Sheffe dimension, i.e. dfn = 1 to generate ordinary CIs, dfn = 2 for 2-dim CE, etc.
# "
#     require(car)
#     which.coef <- if (length(coefficients(model)) == 2)
#         c(1, 2)
#     else {
#         if (missing(which.coef)) {
#             if (has.intercept(model))
#                 c(2, 3)
#             else c(1, 2)
#         }
#         else which.coef
#     }
#     coef <- coefficients(model)[which.coef]
#     xlab <- if (missing(xlab))
#         paste(names(coef)[1], "coefficient")
#     ylab <- if (missing(ylab))
#         paste(names(coef)[2], "coefficient")
#     if(missing(dfn)) {
#         if (Scheffe) dfn <- sum(df.terms(model))
#         else 2
#     }
#     dfd <- df.residual(model)
#     shape <- vcov(model)[which.coef, which.coef]
#     ret <- numeric(0)
#     for (level in rev(sort(levels))) {
#         radius <- sqrt(dfn * qf(level, dfn, dfd))
#         ret <- rbind(ret, c(NA,NA), ell( coef, shape, radius) )
#     }
#     colnames(ret) <- c(xlab, ylab)
#     ret
# }
#
# # from Plot3d.R
#
#
#
#
#
# #' Calculate coordinates of a data ellipse
# #'
# #'
# #' \code{dell} to calculates the coordinates of a 2D data ellipse
# #' (concentration ellipse) from (X, Y) variables.
# #'
# #' \code{dellplus} can produce, in addition to the points of an ellipse, the
# #' conjugate axes corresponding to a \code{chol} or other decomposition and the
# #' surrounding parallelogram defined by these axes.
# #'
# #' These functions simply calculate the mean vector and covariance matrix and
# #' call \code{ell} or \code{ellplus}.
# #'
# #' @aliases dell dellplus
# #' @param x,y Either a two-column matrix or numeric vectors of the same length
# #' @param radius Radius of the ellipse-generating unit circle.  The default,
# #' \code{radius=1} corresponds to a "standard" ellipse.
# #' @param \dots Other arguments passed down to \code{ell} or \code{ellplus}.
# #' @return Returns a 2-column matrix of (X,Y) coordinates suitable for drawing
# #' with \code{lines()}.
# #'
# #' For \code{dellplus}, when more than one of the options \code{ellipse},
# #' \code{diameters}, and \code{box} is \code{TRUE}, the different parts are
# #' separated by a row of \code{NA}.
# #' @author Georges Monette
# #' @seealso \code{\link{cell}}, \code{\link{ell}}, \code{\link{ellplus}},
# #' @references Monette, G. (1990). Geometry of Multiple Regression and
# #' Interactive 3-D Graphics. In Fox, J. & Long, S. (ed.)  \emph{Modern Methods
# #' of Data Analysis}, Sage Publications, 209-256.
# #' @keywords dplot aplot
# #' @examples
# #'
# #' data(Prestige)   # from car
# #' attach(Prestige)
# #' fit.simple <- lm( prestige ~ education, Prestige)
# #'
# #' plot(prestige ~ education, type='p')
# #' lines(dell(education, prestige), col="blue", lwd=3)
# #' lines(bbox <- dellplus(education, prestige, box=TRUE))
# #' lines(dellplus(education, prestige, diameter=TRUE, radius=2), col="gray")
# #' detach(Prestige)
# #'
# #'
# #' @export
#     dell <- function( x, y, radius = 1, ...) {
#         if ( (is.matrix(x) && (ncol(x) > 1))|| is.data.frame(x)) mat <- as.matrix(x[,1:2])
#         else if (is.list(x)) mat <- cbind(x$x, x$y)
#         else mat <- cbind( x,y)
#         ell( apply(mat,2,mean), var(mat), radius = radius, ...)
#
#     }
#
#
#
#
#
# #' Calculate coordinates of an ellipse
# #'
# #' \code{ell} is a utility function used to calculate the (X, Y) coordinates of
# #' a 2D ellipse for the purpose of drawing statistical diagrams and plots.
# #'
# #' \code{ellplus} can produce, in addition to the points of an ellipse, the
# #' conjugate axes corresponding to a \code{chol} or other decomposition and the
# #' surrounding parallelogram defined by these axes.
# #'
# #'
# #' @aliases ell ellplus
# #' @param center X,Y location of the center of the ellipse
# #' @param shape A 2x2 matrix, typically a covariance matrix of data (for a data
# #' ellipse), or a covariance matrix of estimated parameters in a model (for a
# #' confidence ellipse).
# #' @param radius Radius of the ellipse-generating unit circle.  The default,
# #' \code{radius=1} corresponds to a "standard" ellipse.
# #' @param n Number of points on the unit circle used to calculate the ellipse
# #' @param angles Angles around the unit circle used to calculate the ellipse
# #' @param fac A function defining the conjugate axes used to transform the unit
# #' circle into an ellipse.  The default, \code{chol}, uses the right Cholesky
# #' factor of \code{shape}.
# #' @param ellipse Logical to indicate if the points on the ellipse should be
# #' returned
# #' @param diameters Logical to indicate if the points defining the ends of the
# #' conjugate axes of the ellipse should be returned
# #' @param box Logical to indicate if the points on the conjugate-axes bounding
# #' box should be returned
# #' @param all Logical to request all of \code{ellipse}, \code{diameters} and
# #' \code{box}. If \code{FALSE}, only the components specified separately by
# #' \code{ellipse}, \code{diameters} and \code{box} are returned.
# #' @return Returns a 2-column matrix of (X,Y) coordinates suitable for drawing
# #' with \code{lines()}.
# #'
# #' For \code{ellplus}, when more than one of the options \code{ellipse},
# #' \code{diameters}, and \code{box} is \code{TRUE}, the different parts are
# #' separated by a row of \code{NA}.
# #' @author Georges Monette
# #' @seealso \code{\link{cell}}, \code{\link{dell}}, \code{\link{dellplus}},
# #' @keywords dplot aplot
# #' @examples
# #'
# #' plot( x=0,y=0, xlim = c(-3,3), ylim = c(-3,3),
# #'       xlab = '', ylab = '', type = 'n', asp=1)
# #' abline( v=0, col="gray")
# #' abline( h=0, col="gray")
# #' A <- cbind( c(1,2), c(1.5,1))
# #' W <- A %*% t(A)
# #'
# #' lines( ell(center=c(0,0), shape = W ), col = 'blue', lwd=3)
# #' lines( ellplus(center=c(0,0), shape = W, box=TRUE, diameters=TRUE ), col = 'red')
# #'
# #' # show conjugate axes for PCA factorization
# #' pca.fac <- function(x) {
# #'     xx <- svd(x)
# #'     ret <- t(xx$v) * sqrt(pmax( xx$d,0))
# #'     ret
# #' }
# #'
# #' plot( x=0,y=0, xlim = c(-3,3), ylim = c(-3,3),
# #'       xlab = '', ylab = '', type = 'n', asp=1)
# #' abline( v=0, col="gray")
# #' abline( h=0, col="gray")
# #' lines( ell(center=c(0,0), shape = W ), col = 'blue', lwd=3)
# #' lines( ellplus(center=c(0,0), shape = W, box=TRUE, diameters=TRUE, fac=pca.fac ), col = 'red')
# #'
# #'
# #' @export
#     ell <- function( center = rep(0,2) , shape = diag(2), radius  = 1, n =100) {
#           fac <- function( x )  {
#               # fac(M) is a 'right factor' of a positive semidefinite M
#               # i.e. M = t( fac(M) ) %*% fac(M)
#               # similar to chol(M) but does not require M to be PD.
#               xx <- svd(x,nu=0)
#           t(xx$v) * sqrt(pmax( xx$d,0))
#           }
#          angles = (0:n) * 2 * pi/n
#          if ( length(radius) > 1) {
#             ret <- lapply( radius, function(r) rbind(r*cbind( cos(angles), sin(angles)),NA))
#             circle <- do.call( rbind, ret)
#          }
#          else circle = radius * cbind( cos(angles), sin(angles))
#          ret <- t( c(center) + t( circle %*% fac(shape)))
#          attr(ret,"parms") <- list ( center = rbind( center), shape = shape, radius = radius)
#          class(ret) <- "ell"
#          ret
#     }
#
#
#
#
#
#
# #' Centre of an object
# #'
# #'
# #'
# #'
# #'
# #' @param obj
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (obj, ...)
# #' UseMethod("center")
# #'
# #' @export
#     center <- function( obj, ... ) UseMethod("center")
#
#
#
#
#
# #' Center of an ellipse
# #'
# #'
# #'
# #'
# #'
# #' @param obj
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (obj, ...)
# #' attr(obj, "parms")$center
# #'
# #' @export
#     center.ell <- function( obj, ...) attr(obj, 'parms') $ center
#
#
# # line between two ellipsoids ( would like to extend to cylinders)
#
#
#
#
# #'
# #' ellipses
# #'
# #'
# #' Generates points on the curve of osculation between the centers of two
# #' families of ellipses
# #'
# #'
# #'
# #' @aliases osculant.default osculant
# #' @param center1
# #' ellipses.
# #' @param shape1
# #' first family of ellipses.
# #' @param center2 of second family.
# #' @param shape2
# #' @param n n + 1 is the number of points to generate along the locus. \code{n
# #' = 1} generates the two centers, \code{n=0} generates to point on the first
# #' ellipse that lines on the locus of osculation provided the centre of the
# #' second family lines outside the first ellipse.
# #' @param range of values of \code{u} to use to generate points. (See the
# #' algorithm in the code)
#
# #' @author Georges Monette (georges@@yorku.ca)
#
# #' \code{\link{cell}}, \code{\link{dell}},
# #' \code{\link{ellplus}},\code{\link{dellplus}},
#
# #' @keywords ellipse ellipse geometry
# #' @examples
# #'
# #' v1 <- 36*diag(2)
# #' v2 <- .5 * (diag(2) - .4)
# #' v2[2,2] <- 2
# #' plot( 0:10,0:10,type = 'n')
# #' lines( ell( c(2,2), v1))
# #' lines( ell( c(4,4), v2), col = 'red')
# #' osculant(  c(2,2), v1, c(4,4), v2, n = 3)
# #' osculant(  c(2,2), v1, c(4,4), v2, n = 1)
# #' lines( osculant( c(2,2), v1, c(4,4), v2), col = 'red')
# #'
# #' lines( ell( c(8,8), v2), col = 'blue')
# #' lines( osculant( c(2,2), v1, c(8,8), v2), col = 'blue')
# #' points( osculant( c(2,2), v1, c(8,8), v2, n=1),pch = 16, col = 'blue')
# #' points( osculant( c(2,2), v1, c(8,8), v2, n=0),pch = 16, col = 'blue')
# #' points( osculant( c(8,8), v2, c(2,2), v1,  n=0),pch = 16, col = 'blue')
# #'
# #'
# #' @export
# osculant <- function(x, ...) UseMethod("osculant")
#
# osculant.default <- function( center1, shape1, center2, shape2, n = 100, range =c(0,1), maxu = 100) {
# # Use solution from Lagrangean:
# # p = ( shape1^-1 + lam * shape1^-1)^-1 %*% lam2 shape2^-1 delta
#     pt <- function(u)  u* solve( u*diag(p) + (1-u)* shape, delta)
#
#
# #' p - quick paste with sep = ''
# #'
# #' Works like \code{paste}, using an empty separator string.
# #'
# #'
# #' @param \dots one or more R objects, to be converted to character vectors.
# #' @return A character vector of the concatenated values.
# #' @author Georges Monette
# #' @seealso \code{\link[base]{paste}}
# #' @keywords manip
# #' @examples
# #'
# #' p(letters[1:5], 1:5)
# #'
#     p <- nrow(shape1)
#     delta <- center2 - center1
#     shape <- t(solve(shape1,shape2))
#     shape <- shape/mean(diag(shape))   # attempt to equalize intervals over range
#     if( n == 0) {
#       norm1 <- function( u ) {
#           pp <- pt(u)
#           sqrt( sum( pp  * solve( shape1, pp))) -1
#       }
#       if ( norm1(1) < 0 ) {
#           warning( "Center of second ellipse inside first ellipse")
#           return( NULL)
#       }
#       u <- uniroot( norm1, c(0,1))$root
#       rbind( pt(u) + center1)
#
#     } else {
#       vec <- sapply( seq(range[1],range[2],diff(range)/n), pt)
#       t( vec + center1 )
#     }
# }
#
# #
# #
# #
# # v1 <- 36*diag(2)
# # v2 <- .5 * (diag(2) - .4)
# # v2[2,2] <- 2
# # plot( 0:10,0:10,type = 'n')
# # lines( ell( c(2,2), v1))
# # lines( ell( c(4,4), v2), col = 'red')
# # osculant(  c(2,2), v1, c(4,4), v2, n = 3)
# # osculant(  c(2,2), v1, c(4,4), v2, n = 1)
# # lines( osculant( c(2,2), v1, c(4,4), v2), col = 'red')
# #
# # lines( ell( c(8,8), v2), col = 'blue')
# # lines( osculant( c(2,2), v1, c(8,8), v2), col = 'blue')
# # points( osculant( c(2,2), v1, c(8,8), v2, n=1),pch = 16, col = 'blue')
# # points( osculant( c(2,2), v1, c(8,8), v2, n=0),pch = 16, col = 'blue')
# # points( osculant( c(8,8), v2, c(2,2), v1,  n=0),pch = 16, col = 'blue')
# #
#
#
#
#
#
#
# #' Confidence ellipse
# #'
# #'
# #'
# #'
# #'
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (...)
# #' {
# #'     UseMethod("cell")
# #'     help <- "\n        See help for car::confidence.ellipse.lm\n        except that 'cell' returns the points to form the ellipse\n        which must be plotted with plot(...,type='l') or lines(...)\n        -- Use dfn to determine Sheffe dimension, i.e. dfn = 1 to generate ordinary CIs, dfn = 2 for 2-dim CE, etc.\n        -- TODO: extend to 3 dimensions if which.coef has length 3\n        "
# #'   }
# #'
# #' @export
# cell <- function( ... )  {
#         UseMethod("cell")
#         help <- "
#         See help for car::confidence.ellipse.lm
#         except that 'cell' returns the points to form the ellipse
#         which must be plotted with plot(...,type='l') or lines(...)
#         -- Use dfn to determine Sheffe dimension, i.e. dfn = 1 to generate ordinary CIs, dfn = 2 for 2-dim CE, etc.
#         -- TODO: extend to 3 dimensions if which.coef has length 3
#         "
# }
#
#
# cell.wald <-
#  function (obj, which.coef = 1:2, levels = 0.95, Scheffe = FALSE, dfn = 2,
#     center.pch = 19, center.cex = 1.5, segments = 51, xlab, ylab,
#     las = par("las"), col = palette()[2], lwd = 2, lty = 1,
#     add = FALSE, ...)
# {
#
# # BUGS: works only on first element of glh list
# # glh should be restructured to have two classes: waldList and wald
#
#     obj <- obj[[1]]
#     coef <- obj$coef[which.coef]
#     xlab <- if (missing(xlab))
#         paste(names(coef)[1], "coefficient")
#     ylab <- if (missing(ylab))
#         paste(names(coef)[2], "coefficient")
#
#     dfd <- obj$anova$denDF
#     shape <- obj$vcov[which.coef, which.coef]
#     ret <- ell( coef, shape , sqrt( dfn * qf( levels, dfn, dfd)))
#     colnames(ret) <- c(xlab, ylab)
#     ret
# }
#
#
#
#
#
#
# #' Confidence ellipse
# #'
# #'
# #'
# #'
# #'
# #' @param model
# #' @param which.coef
# #' @param levels
# #' @param Scheffe
# #' @param dfn
# #' @param center.pch
# #' @param center.cex
# #' @param segments
# #' @param xlab
# #' @param ylab
# #' @param las
# #' @param col
# #' @param lwd
# #' @param lty
# #' @param add
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (model, which.coef, levels = 0.95, Scheffe = FALSE,
# #'     dfn = 2, center.pch = 19, center.cex = 1.5, segments = 51,
# #'     xlab, ylab, las = par("las"), col = palette()[2], lwd = 2,
# #'     lty = 1, add = FALSE, ...)
# #' {
# #'     which.coef <- if (length(coefficients(model)) == 2)
# #'         c(1, 2)
# #'     else {
# #'         if (missing(which.coef)) {
# #'             if (any(names(coefficients(model)) == "(Intercept)"))
# #'                 c(2, 3)
# #'             else c(1, 2)
# #'         }
# #'         else which.coef
# #'     }
# #'     coef <- coefficients(model)[which.coef]
# #'     xlab <- if (missing(xlab))
# #'         paste(names(coef)[1], "coefficient")
# #'     ylab <- if (missing(ylab))
# #'         paste(names(coef)[2], "coefficient")
# #'     if (missing(dfn)) {
# #'         if (Scheffe)
# #'             dfn <- sum(df.terms(model))
# #'         else 2
# #'     }
# #'     dfd <- df.residual(model)
# #'     shape <- vcov(model)[which.coef, which.coef]
# #'     ret <- numeric(0)
# #'     ret <- ell(coef, shape, sqrt(dfn * qf(levels, dfn, dfd)))
# #'     colnames(ret) <- c(xlab, ylab)
# #'     ret
# #'   }
# #'
# cell.default <-
# function (model, which.coef, levels = 0.95, Scheffe = FALSE, dfn = 2,
#     center.pch = 19, center.cex = 1.5, segments = 51, xlab, ylab,
#     las = par("las"), col = palette()[2], lwd = 2, lty = 1,
#     add = FALSE, ...)
# {
#
#     #require(car)
#     which.coef <- if (length(coefficients(model)) == 2)
#         c(1, 2)
#     else {
#         if (missing(which.coef)) {
#             if (any(names(coefficients(model)) == "(Intercept)"))
#                 c(2, 3)
#             else c(1, 2)
#         }
#         else which.coef
#     }
#     coef <- coefficients(model)[which.coef]
#     xlab <- if (missing(xlab))
#         paste(names(coef)[1], "coefficient")
#     ylab <- if (missing(ylab))
#         paste(names(coef)[2], "coefficient")
#     if(missing(dfn)) {
#         dfn <- if (Scheffe) sum(df.terms(model)) else 2
#     }
#     dfd <- df.residual(model)
#     shape <- vcov(model)[which.coef, which.coef]
#     ret <- numeric(0)
#
#     ret <- ell( coef, shape,sqrt(dfn * qf(levels, dfn, dfd)))
#     colnames(ret) <- c(xlab, ylab)
#     ret
# }
#
#
#
#
#
#
#
#
#
#
#
#
#
# "</pre>
#
# == Diagnostics ==
# <pre>"
#
#
#
#
# #' Generic diagnostics
# #'
# #' Generic diagnostics
# #'
# #'
# #'
# #' @param x
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, ...)
# #' UseMethod("diags")
# #'
# #' @export
# diags <- function(x, ...) UseMethod("diags")
#
#
#
#
# #' Standard diagnostics for lm objects
# #'
# #' Standard diagnostics for lm objects
# #'
# #'
# #'
# #' @param x
# #' @param y
# #' @param \dots
# #' @param ask
# #' @param labels
# #' @param showlabs
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, y, ..., ask, labels = names(residuals(x)), showlabs = text)
# #' {
# #'     if (!missing(ask)) {
# #'         op <- par(ask = ask)
# #'         on.exit(par(op))
# #'     }
# #'     form <- formula(x)
# #'     f <- predict(x)
# #'     r <- residuals(x)
# #'     nams <- names(r)
# #'     if (!missing(labels)) {
# #'         nams <- names(residuals(x))
# #'         if (length(nams) != length(labels))
# #'             labels <- labels[nams]
# #'     }
# #'     ret <- NULL
# #'     if (missing(y)) {
# #'         y <- f + r
# #'         yname <- deparse(form[[2]])
# #'     }
# #'     else yname <- deparse(substitute(y))
# #'     fname <- paste("Fitted:", deparse(form[[3]]), collapse = " ")
# #'     plot(f, y, xlab = fname, ylab = yname, main = "Dependent var. vs. Predicted",
# #'         ...)
# #'     abline(0, 1, lty = 1)
# #'     lines(supsmu(f, y))
# #'     showlabs(f, y, labels, ...)
# #'     lmi <- lm.influence(x)
# #'     hat <- lmi$hat
# #'     sigma <- lmi$sigma
# #'     mm <- scale(model.matrix(x), scale = F)
# #'     mp <- predict(x, type = "terms")
# #'     comp.res <- mp + r
# #'     plot(f, abs(r), xlab = fname, ylab = deparse(substitute(abs(resid(x)))),
# #'         main = "Absolute Residual vs. Predicted", ...)
# #'     showlabs(f, abs(r), labels, ...)
# #'     zq <- qqnorm(r, main = "Normal Quantile Plot", ylab = "Residual",
# #'         sub = fname)
# #'     qqline(r)
# #'     showlabs(zq, labels, ...)
# #'     n <- length(r)
# #'     r.o <- sort(r)
# #'     half <- (n + 1)/2
# #'     if (n%%2 == 1) {
# #'         med <- r.o[half]
# #'         below <- med - r.o[half:1]
# #'         above <- r.o[half:n] - med
# #'     }
# #'     else {
# #'         med <- sum(r.o[c(half, half + 1)])/2
# #'         below <- med - r.o[(n/2):1]
# #'         above <- r.o[(n/2 + 1):n] - med
# #'     }
# #'     opt <- par(pty = "s")
# #'     ran <- range(c(below, above))
# #'     plot(below, above, main = "Symmetry plot of residuals", xlab = "Distance below median",
# #'         ylab = "Distance above median", xlim = ran, ylim = ran)
# #'     abline(0, 1, lty = 2)
# #'     par(opt)
# #'     std.r <- r/(sigma * sqrt(1 - hat))
# #'     plot(hat, std.r, xlab = "Leverage (hat)", ylab = yname, sub = fname,
# #'         main = "Studentized residual vs. Leverage", ...)
# #'     showlabs(hat, std.r, labels, ...)
# #'     nams <- dimnames(lmi$coefficients)[[1]]
# #'     pairs(lmi$coefficients)
# #'     pairs(lmi$coefficients, panel = function(x, y, nams) {
# #'         points(x, y)
# #'         text(x, y, nams)
# #'     }, nams = nams)
# #'     invisible(0)
# #'   }
# #'
# #' @export
# diags.lm <- function(x, y, ..., ask, labels = names(residuals(x)), showlabs = text)
#     {
#     # diags.lm
#     # graphical diagnostics for lm, locally first-order for glm
#     # enlarged version of plot.lm with emphasis on diagnostics
#     # G. Monette, Dec. 94
#     # modified Nov. 97, May 98
#     # Slight modification to pairs adding labels, Jan 03
# 	if(!missing(ask)) {
# 		op <- par(ask = ask)
# 		on.exit(par(op))
# 	}
# 	form <- formula(x)
# 	f <- predict(x)
# 	r <- residuals(x)
# 	nams <- names(r)
# 	if(!missing(labels)) {
# 		nams <- names(residuals(x))	#
#     # if labels not same length as residuals assume it's a vector
#     # of len == original data and select elements included in residuals
# 		if(length(nams) != length(labels))
# 			labels <- labels[nams]
# 	}
# 	ret <- NULL
# 	if(missing(y)) {
# 		y <- f + r
# 		yname <- deparse(form[[2]])
# 	}
# 	else yname <- deparse(substitute(y))
# 	fname <- paste("Fitted:", deparse(form[[3]]), collapse = " ")
# 	plot(f, y, xlab = fname, ylab = yname, main = "Dependent var. vs. Predicted",
# 		...)
# 	abline(0, 1, lty = 1)
# 	lines(supsmu(f,y))
# 	showlabs(f, y, labels,...)
#     #
#     # get influence diags and model matrix while looking at first plot
#     #
# 	lmi <- lm.influence(x)
# 	hat <- lmi$hat
# 	sigma <- lmi$sigma	# drop 1 sigma
# 	mm <- scale(model.matrix(x), scale = F)	# centres each column
# 	mp <- predict(x, type = "terms")
# 	comp.res <- mp + r	# effect + residual
#     #
#     # Absolute residual vs. predicted
#     #
# 	plot(f, abs(r), xlab = fname, ylab = deparse(substitute(abs(resid(x)))),
# 		main = "Absolute Residual vs. Predicted", ...)
# 	showlabs(f, abs(r), labels, ...)	#
#     #
#     # Normal quantile plot
#     #
# 	zq <- qqnorm(r, main = "Normal Quantile Plot", ylab = "Residual", sub
# 		 = fname)
# 	qqline(r)
# 	showlabs(zq, labels,...)	#
#     #
#     # Symmetry plot of residuals (Lawrence C. Hamilton, Regression with
#     #       Graphics, Duxbury, 1992)
# 	n <- length(r)
# 	r.o <- sort(r)
# 	half <- (n + 1)/2
# 	if (n%%2 == 1) {    # n is odd
# 		med <- r.o[half]
# 		below <- med - r.o[half:1]
# 		above <- r.o[half:n] - med
# 	}
# 	else {
#     # n is even
# 		med <- sum(r.o[c(half, half + 1)])/2
# 		below <- med - r.o[(n/2):1]
# 		above <- r.o[(n/2 + 1):n] - med
# 	}
# 	opt <- par(pty = "s")
# 	ran <- range(c(below, above))
# 	plot(below, above, main = "Symmetry plot of residuals", xlab =
# 		"Distance below median", ylab = "Distance above median", xlim
# 		 = ran, ylim = ran)
# 	abline(0, 1, lty = 2)
# 	par(opt)	#
#
#     #
#     # Studentized residual vs. leverage
#     #
#
# 	std.r <- r/(sigma * sqrt(1 - hat))
# 	plot(hat, std.r, xlab = "Leverage (hat)", ylab = yname, sub = fname,
# 		main = "Studentized residual vs. Leverage", ...)
# 	showlabs(hat, std.r, labels,...)	#	plot(lmi$sig, std.r)	#
#
#     #
#     # effect of dropping one observation DFBETA
#     #
#
# 	nams <- dimnames(lmi$coefficients)[[1]]
# 	pairs(lmi$coefficients)
# 	pairs(lmi$coefficients, panel = function(x,y,nams){
# 		points(x,y)
# 		text(x,y,nams)
# 	}, nams = nams)
#
# 	# main = "Effect of dropping one case", sub = fname)
# 	invisible(0)
# }
#
# #' @export
# model.matrix.lme <- function( fit , data = fit$data,
#                               na.action = fit$na.action,
#                               ...){
#   mCall <- fit$call
#   fixed <- eval(eval(mCall$fixed)[-2])
#   data <- model.frame(fixed, data = data)
#   naresid( na.action, model.matrix(fixed, data = data))
# }
# # model.matrix(fit, data=pred) %>% dim
# # model.frame(fit, data=pred) %>% dim
# # model.matrix(fit, na.action = na.pass) %>% dim
# # model.frame(fit, na.action = na.pass) %>% dim
#
#
# # BUG: not working as it should for na.action=na.exclude
# # model.frame.lme <- function( fit , data = fit$data,
# #                              na.action = fit$call$na.action,...)
# #   model.frame(formula(fit), data = data, na.action = na.action)
#
# #' @export
# model.frame.lme <- function (object, data =object$data, na.action = object$na.action,
#                              ...)
# {
#   # adapted from portions of predict.lme
#   mCall <- object$call
#   fixed <- eval(eval(mCall$fixed)[-2])
#   Terms <- object$terms
#   data <- as.data.frame(data)
#   mfArgs <- list(formula = fixed,
#                  data = data, na.action = na.action, drop.unused.levels = TRUE)
#   dataMix <- do.call("model.frame", mfArgs)
#   dataMix
# }
#
#
#
#
#
# #' Standard diagnostics for lme objects
# #'
# #'
# #'
# #'
# #'
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (...)
# #' cat("Being implemented")
# #'
# #' @export
# diags.lme <- function( ... ) cat("Being implemented")
#
#
#
# "
# </pre>
# == vplot -- a plot function for matrix algebra ==
# <pre>
# "
#
# ##
# ##  vplot  plots columns of a 2 x n matrix
# ##  Transferred to coursfun: Nov. 15, 2005
#
#
#
# #' Collection of functions to help teach matrix geometry in 2 dimensions
# #'
# #'
# #'
# #' vplot - plots the columns of a 2 x n matrix or a vector of length 2 - vplot
# #' adds to the current plot resizing it to include all plotted objects in a
# #' 'euclidean' frame - to start a new plot, use 'new = T' - to remove the last
# #' element added use 'vplot(pop=1)' Associated functions: - vell( mean, var)
# #' generates an ellipse, default = unit circle - vbox() generates a box -
# #' vobj() generates a circle in a box - orthog(theta) generates an orthog
# #' matrix rotating through angle theta - orthog.proj generates the matrix of an
# #' orthog. projection into span (x) - vmat( .... ) generates a 2 by n matrix
# #' Examples: vplot( new = T ) vplot( vell(), 'l' ) vplot( cbind(c(3,1),c(1,4))
# #' \%*\% vell()) vplot( pop = 1) vplot( cbind(c(3,1),c(1,4)) \%*\% vell(), type
# #' = 'l', col = 'red')
# #' above ~~
# #'
# #' @param mat
# #' @param type
# #' @param new
# #' @param pch
# #' @param pop
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (mat, type = "p", new = F, pch = 16, pop = 0, ...)
# #' {
# #'
# #'
# #'     if (new || !exists(".vplot"))
# #'         assign(".vplot", list(list(x = 0, y = 0, type = "n")),
# #'             pos = 1)
# #'     a <- .vplot
# #'     if (!missing(mat)) {
# #'         mat <- cbind(mat)
# #'         if (type == "v") {
# #'             zz <- rbind(0 * mat, mat, mat/0)
# #'             mat <- matrix(c(zz), nrow = 2)
# #'             type = "b"
# #'         }
# #'         d <- dim(mat)
# #'         if (d[1] != 2 && d[2] == 2) {
# #'             mat <- t(mat)
# #'             warning("mat is n x 2 and has been transposed")
# #'         }
# #'         a <- c(a, list(list(x = mat[1, ], y = mat[2, ], type = type,
# #'             pch = pch, ...)))
# #'     }
# #'     dat <- NULL
# #'     for (i in seq(along = a)) {
# #'         dat <- c(dat, a[[i]]$x, a[[i]]$y)
# #'     }
# #'     par(pty = "s")
# #'     plot(range(na.omit(dat)), range(na.omit(dat)), type = "n",
# #'         xlab = "", ylab = "")
# #'     if (pop > 0) {
# #'         keep <- 1:max(1, (length(a) - (pop + 1)))
# #'         a <- a[keep]
# #'     }
# #'     abline(h = 0, v = 0)
# #'     for (i in seq(along = a)) do.call("points", a[[i]])
# #'     assign(".vplot", a, pos = 1)
# #'     invisible(a)
# #'   }
# #'
# #' @export
# vplot <- function( mat , type = 'p', new = F,  pch = 16, pop = 0, ...) {
# help <- "
# vplot    - plots the columns of a 2 x n matrix or a vector of length 2
#          - vplot adds to the current plot resizing it to include all plotted
#            objects in a 'euclidean' frame
#          - to start a new plot, use 'new = TRUE'
#          - to remove the last element added use 'vplot(pop=1)'
#          Associated functions:
#          - vell( mean, var) generates an ellipse, default = unit circle
#          - vbox() generates a box
#          - vobj() generates a circle in a box
#          - orthog(theta) generates an orthog matrix rotating through angle theta
#          - orthog.proj generates the matrix of an orthog. projection into span (x)
#          - vmat( .... ) generates a 2 by n matrix
#          Examples:
#            vplot( new = TRUE )
#            vplot( vell(), 'l' )
#            vplot( cbind(c(3,1),c(1,4)) %*% vell())
#            vplot( pop = 1)
#            vplot( cbind(c(3,1),c(1,4)) %*% vell(), type = 'l', col = 'red')
# "
#      if (  new || !exists(".vplot")) assign(".vplot", list(list(x=0,y=0,type='n')),pos=1)
#      a <- .vplot
#      if ( ! missing(mat) ) {
#         mat <- cbind(mat)
#         if ( type == 'v' ) {
#            zz <- rbind( 0*mat, mat, mat/0)
#            mat <- matrix( c(zz), nrow = 2)
#            type = 'b'
#         }
#         d <- dim(mat)
#         if ( d[1] != 2 && d[2] == 2){
#            mat <- t(mat)
#            warning("mat is n x 2 and has been transposed")
#         }
#         a <- c(a,list( list(x=mat[1,],y = mat[2,],type=type, pch = pch, ...)))
#      }
#      dat <- NULL
#      for ( i in seq( along = a )) {
#          dat <- c( dat, a[[i]]$x, a[[i]]$y)
#      }
#      # print(a)
#      par ( pty = 's')
#      plot( range(na.omit(dat)), range(na.omit(dat)), type = 'n', xlab = '', ylab ='')
#      if ( pop > 0 ) {
#         keep <- 1:max(1,(length(a)-(pop+1)))
#         a <- a[keep]
#      }
#      abline( h = 0, v = 0)
#      for ( i in seq( along = a)) do.call('points', a[[i]])
#      assign(".vplot", a, pos = 1)
#      invisible(a)
# }
#
#
#
# #' Vector around an ellipse
# #'
# #'
# #'
# #'
# #'
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (...)
# #' t(ell(...))
# #'
# #' @export
# vell <- function(...) t( ell(...))
#
#
# #' Unit box
# #'
# #'
# #'
# #'
# #'
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (...)
# #' cbind(c(-1, -1), c(-1, 1), c(1, 1), c(1, -1), c(-1, -1))
# #'
# #' @export
# vbox <- function(...) cbind( c(-1,-1), c(-1,1), c(1,1), c(1,-1), c(-1,-1))
#
#
# #' Combine ellipse with subtending box
# #'
# #'
# #'
# #'
# #'
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (...)
# #' {
# #'     cbind(vell(), NA, vbox(), NA, c(0, -1), c(0, 1), NA, c(-1,
# #'         0), c(1, 0))
# #'   }
# #'
# #' @export
# vobj <- function(...) {
#      cbind( vell(), NA, vbox(), NA, c(0,-1),c(0,1), NA, c(-1,0), c(1,0))
# }
#
#
# #' Square in 2 dimensions
# #'
# #'
# #'
# #'
# #'
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (...)
# #' vmat(0, 0, 0, 1, 1, 1, 1, 0, 0, 0)
# #'
# #' @export
# vsquare <- function(...) vmat( 0,0,0,1,1,1,1,0,0,0)
#
#
#
# #' Create a matrix entering vectors column by column
# #'
# #'
# #'
# #'
# #'
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (...)
# #' {
# #'     help <- "\nvmat creates a matrix entering data column by column\n"
# #'     aa <- list(...)
# #'     aa <- do.call("c", aa)
# #'     matrix(aa, nrow = 2)
# #'   }
# #'
# #' @export
# vmat <- function(...) {
# help <- "
# vmat creates a matrix entering data column by column
# "
#      aa <- list(...)
#      aa <- do.call('c', aa)
#      matrix(aa, nrow = 2)
# }
#
#
#
# #' Two by two matrix of orthogonal rotation
# #'
# #'
# #'
# #'
# #'
# #' @param theta
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (theta)
# #' cbind(c(cos(theta), sin(theta)), c(-sin(theta), cos(theta)))
# #'
# #' @export
# orthog <- function( theta ) cbind( c( cos(theta), sin(theta)), c( - sin(theta), cos(theta)))
#
#
# #' Orthogonal projection matrix -- check possibly poor numerical behaviour
# #'
# #'
# #'
# #'
# #'
# #' @param x
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x)
# #' {
# #'     x <- cbind(x)
# #'     x %*% solve(t(x) %*% x, x)
# #'   }
# #'
# #' @export
# orthog.proj <- function ( x ) {
#        x <- cbind(x)
#        x %*% solve(t(x) %*% x , x)
# }
#
#
# "
# </pre>
# == Read.spss and Read.dta ==
# <pre>
# "
#
# ###
# ###  trim
# ###
#
#
#
# #' Trim trailing blanks
# #'
# #' Generic function to trim leading and trailing blanks from character vectors
# #' and factors.
# #'
# #' The main application is in reading SPSS files that often have leading or
# #' trailing blanks in character and factor values. These blanks are often
# #' inconsistent so that values will appear to differ even though they are
# #' equal.  The trim function is called in \code{Read.spss} to remove leading
# #' and trailing blanks from all factors.
# #'
# #' @param x a data frame, factor, character or numeric vector
# #' \code{x} here~~
# #' @return A character vector or factor with leading and trailing blanks
# #' removed.
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x)
# #' {
# #' trim in fun.R
# #'     UseMethod("trim")
# #'   }
# #'
# #' @export
#   trim <- function(x) {
#        help <- "
# trim in fun.R
#   removes trailing blanks from character variables or from
#   factor levels
#   Use mainly to trim .dta files produced with SPSS
# "
#        UseMethod("trim")
# }
#
#
#
# #' Trim trailing blanks from all variables in a data frame
# #'
# #'
# #'
# #'
# #'
# #' @param x
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x)
# #' {
# #'     for (nn in names(x)) x[[nn]] <- trim(x[[nn]])
# #'     x
# #'   }
# #'
# #' @export
#   trim.data.frame <- function(x) {
#       for ( nn  in names(x)) x[[nn]] <- trim(x[[nn]])
#       x
#   }
#
#
# #' Trim trailing blanks from a factor object
# #'
# #'
# #'
# #'
# #'
# #' @param x
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x)
# #' {
# #'     levels(x) <- trim(levels(x))
# #'     x
# #'   }
# #'
# #' @export
#   trim.factor <- function( x ) {
#       levels(x) <- trim(levels(x))
#       x
#   }
#
#
# #' Trim trailing blanks from character vector
# #'
# #'
# #'
# #'
# #'
# #' @param x
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x)
# #' {
# #'     x[] <- sub(" +$", "", x)
# #'     x
# #'   }
# #'
# #' @export
#   trim.character <- function( x ) {
#       x[] <- sub(" +$", "", x )
#       x[] <- sub("^ +", "", x )
#       x
#   }
#
#
# #' Trim default -- identity
# #'
# #'
# #'
# #'
# #'
# #' @param x
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x)
# #' x
# #'
# #' @export
#   trim.default <- function(x) x
#
#
#
# #' Read an SPSS file
# #'
# #'
# #'
# #'
# #'
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (...)
# #' {
# #'     require("Hmisc")
# #'     dd <- spss.get(...)
# #'     trim(dd)
# #'   }
# #'
# #' @export
#   Read.spss <- function( ... ) {
#             require("Hmisc")
#             dd <- spss.get ( ... )
#             trim( dd )
#   }
#
#
#
# #' read a STATA .dta file
# #'
# #'
# #'
# #'
# #'
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (...)
# #' {
# #'     help <- "\n  Read.dta reads Stata files using 'read.dta' in 'library(foreign)'\n  This appears to be an ideal way of importing spss files in order\n  to keep full variable names. Direct use of 'read.spss' on a SPSS\n  '.sav' file abbreviates variable names to 8 characters.\n  Note: missing.type = T produces warnings.\n"
# #'     require("foreign")
# #'     dd <- read.dta(...)
# #'     cls <- sapply(dd, class)
# #'     ch.nams <- names(dd)[cls == "character"]
# #'     for (nn in ch.nams) dd[[nn]] <- factor(trim(dd[[nn]]))
# #'     dd
# #'   }
# #'
# #' @export
#   Read.dta <- function ( ... ) {
# help <- "
#   Read.dta reads Stata files using 'read.dta' in 'library(foreign)'
#   This appears to be an ideal way of importing spss files in order
#   to keep full variable names. Direct use of 'read.spss' on a SPSS
#   '.sav' file abbreviates variable names to 8 characters.
#   Note: missing.type = TRUE produces warnings.
# "
#            require("foreign")
#            #  dd <- read.dta(... , missing.type = TRUE)  # Note: missing.type = T produces warnings.
#            dd <- read.dta(...)
#            cls <- sapply(dd,class)
#            ch.nams <- names(dd) [ cls == "character" ]
#            for ( nn in ch.nams ) dd[[nn]] <- factor(trim(dd[[nn]]) )
#            dd
#   }
#
#
#
#
# #' Write a STATA file
# #'
# #'
# #'
# #'
# #'
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (...)
# #' {
# #'     require("foreign")
# #'     write.dta(..., version = 7)
# #'   }
# #'
# #' @export
#   Write.dta <- function( ... ) {
#        require("foreign")
#        write.dta( ..., version = 7)
#   }
#
#
#
# #' Write and SPSS file
# #'
# #'
# #'
# #'
# #'
# #' @param dataframe
# #' @param file
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (dataframe, file, ...)
# #' {
# #'     require(foreign)
# #'     dname <- deparse(substitute(dataframe))
# #'     disp(dname)
# #'     cls <- sapply(dataframe, class)
# #'     for (nn in names(dataframe)[cls == "Date"]) {
# #'         dataframe[[nn]] <- as.character(dataframe[[nn]], "%Y/%m/%d")
# #'     }
# #'     if (any(cls == "Date")) {
# #'         cat("\nOpen .dta file in SPSS and convert following variables to dates\nwith yyyy/mm/dd format:\n")
# #'         for (nn in names(dataframe)[cls == "Date"]) cat("       ",
# #'             nn, "\n")
# #'     }
# #'     for (nn in names(dataframe)[cls == "factor"]) {
# #'         dataframe[[nn]] <- as.character(dataframe[[nn]])
# #'     }
# #'     if (missing(file))
# #'         file <- paste(dname, ".dta", sep = "")
# #'     else file <- sub("\.dta$|\.DTA$|$", ".dta", file)
# #'     cat(paste("\nData saved in", file, "\n"))
# #'     write.dta(dataframe, file, version = 7, ...)
# #'   }
# #'
# #' @export
#   Write.spss <- function( dataframe, file, ... ) {
#         require(foreign)
#         dname <- deparse(substitute(dataframe))
#         disp(dname)
#         cls <- sapply( dataframe, class )
#         for ( nn in names(dataframe)[cls == "Date"] ){
#             dataframe[[nn]] <- as.character( dataframe[[nn]], "%Y/%m/%d")
#         }
#         if ( any ( cls == "Date")) {
#             cat("\nOpen .dta file in SPSS and convert following variables to dates\nwith yyyy/mm/dd format:\n")
#             for ( nn in names(dataframe) [ cls == "Date" ] ) cat("       ", nn, "\n")
#         }
#         for ( nn in names(dataframe)[cls == "factor"]) {
#             dataframe[[nn]] <- as.character( dataframe[[nn]])
#         }
#         if ( missing(file) )  file <- paste(dname,".dta", sep ="")
#         else file <- sub("\\.dta$|\\.DTA$|$",".dta",file)
#
#         cat(paste("\nData saved in", file,"\n"))
#         write.dta( dataframe, file, version = 7, ...)
#   }
#
#    # zd <- data.frame( x=1:10, a = LETTERS[1:10], d=seq(as.Date("2000-01-01"), by ="4 months", length.out = 10))
#    # zd
#    # Write.spss(zd)
#
#
#
#
#
# "
# </pre>
# == RDC functions -- needs to be organized ==
# <pre>
# "
#
#
#
#
#

# #' @export
# oldcapply.default <- function ( x, by, FUN , ...) {
#   FUN <- match.fun(FUN)
#   if (inherits(by,'formula')) by <- model.frame( by , x , na.action = na.include)
#   ret <- unsplit ( lapply ( split ( x , by ), FUN, ...), by )
#   if ( !is.null( dim(ret)) && length(dim(ret)) ==1) ret <- c(ret)
#   ret
# }
#
#
# #' @export
# xapply <- function(x, ...) UseMethod("xapply")
#
#
#
# #' Apply a function over a ragged array and return a data frame
# #'
# #' Splits the data into subsets, computes summary statistics for each, and
# #' returns the result in a convenient form.
# #' description of what the function does. ~~
# #'
# #' \code{xapply} works like \code{\link{aggregate}} except that the result is
# #' returned as a vector in a data frame with the levels of the \code{by}
# #' variable replicated according to the length of the result.
# #'
# #' The intention in writing \code{xapply} was to facilitate the creation of
# #' 'prediction data frames' extending \code{expand.grid} so that the values of
# #' a variable can easily depend on the value of a factor. For example,
# #' predicting weight from height it might be desired to have a range of heights
# #' for each sex that is consistent with the conditional distribution.
# #'
# #' Suppose a data frame \code{hw} contains variables Sex, height and weight.
# #' Instead of \preformatted{ > fit <- lm ( weight ~ height * Sex, hw) > pred <-
# #' expand.grid( Sex = levels(hw$Sex), height = quantile( hw$height,
# #' c(0,5,25,50,75,95,100)/100)) > pred$weight <- predict( fit, pred) > xyplot(
# #' weight ~ height, pred, groups = Sex, type = 'b') } we can have:
# #' \preformatted{ > fit <- lm ( weight ~ height * Sex, hw) > pred <-
# #' expand.grid( Sex = levels(hw$Sex) ) > pred <- merge( pred, xapply( height ~
# #' Sex, hw, quantile, c(0,5,25,50,75,95,100)/100)) > pred$weight <- predict(
# #' fit, pred) > xyplot( weight ~ height, pred, groups = Sex, type = 'b') }
# #'
# #' @aliases xapply.formula xapply xpandlist
# #' @param formula a two-sided formula: the lhs identifies the variable(s) that
# #' are the first argument of FUN (if there is more than one as in \code{cbind(
# #' x, y)}), then \code{FUN} is applied to each in turn. Note that this
# #' behaviour is different from that of \code{capply}, where \code{FUN} is
# #' applied to the resulting data frame. The rhs identifies the variables to be
# #' used to split the data into subsets.
# #' @param data a data frame.
# #' @param FUN a function.
# #' @param \dots additional arguments to the function
# #' @param subset a condition to subset the data frame.
# #' @param na.action determines whether to omit or include NAs in the grouping
# #' factors. Use \code{na.include} so NAs will form a distinct level. %%
# #' ~~Describe \code{na.action} here~~
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #'
# #' @export
# xapply.formula <- function(formula, data, FUN, ..., subset, na.action = na.omit, debug = FALSE) {
# # the first portion of this code is from stats:::aggregate.formula
# # to avoid evaluating formula
#       FUN <- match.fun(FUN)
#     if (missing(formula) || !inherits(formula, "formula"))
#         stop("'formula' missing or incorrect")
#     if (length(formula) != 3L)
#         stop("'formula' must have both left and right hand sides")
#     m <- match.call(expand.dots = FALSE)
#     if (is.matrix(eval(m$data, parent.frame())))
#         m$data <- as.data.frame(data)
#     m$... <- m$FUN <- NULL
#     m[[1L]] <- as.name("model.frame")
#     if (as.character(formula[[2L]] == ".")) {
#         rhs <- unlist(strsplit(deparse(formula[[3L]]), " *[:+] *"))
#         lhs <- sprintf("cbind(%s)", paste(setdiff(names(data),
#             rhs), collapse = ","))
#         m[[2L]][[2L]] <- parse(text = lhs)[[1L]]
#     }
#     mf <- eval(m, parent.frame())
#     if( debug ) disp(str(mf))
#     if (is.matrix(mf[[1L]])) {
#         lhs <- as.data.frame(mf[[1L]])
#         names(lhs) <- as.character(m[[2L]][[2L]])[-1L]
#         if( debug ) disp( str(lhs) )
#         ret <- aggregate.data.frame(lhs, mf[-1L], FUN = FUN, ..., simplify = FALSE)
#     }
#     else ret <- aggregate.data.frame(mf[1L], mf[-1L], FUN = FUN, ..., simplify = FALSE)
#     xpandlist(ret)
# }
#
#
# #' @export
# xapply.default <- function( x, ...) {
#     ret <- aggregate( x, ..., simplify = FALSE)
#     xpandlist(ret)
# }
#
#
#
# #' @export
# xpandlist <- function ( x ) {
# # make lists long in a data frame
# # internal function to spida to turn the ouput of  aggregate with simplify = FALSE
# #
#   # expand lists:
#     list.vars <- names(x)[sapply(x, is.list)]
#     if ( length( list.vars ) == 0) return(x)
#     nn <- list.vars[1]
#     v <- x[[nn]]
#     ns <- sapply( v, length)
#     val <- unlist(v)
#     ret <- x[ rep( 1:nrow(x), ns),]
#     ret[[nn]] <- val
#     if (length( list.vars) > 1){
#       for ( nn in list.vars[-1]){
#             v <- x[[nn]]
#             ns2 <- sapply( v, length)
#             if ( ! all.equal(ns,ns2)) stop("Length of results for different variables do not have compatible lengths")
#             val <- unlist(v)
#             ret[[nn]] <- val
#       }
#     }
#     ret
# }
#
#
# ###
# ###   pchisq.mix
# ###
#
#
#
# #' P-value for a mixed Chi-Square
# #'
# #'
# #'
# #'
# #'
# #' @param q
# #' @param df
# #' @param mix
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (q, df, mix = rep(1, length(df))/length(df))
# #' {
# #'     pc <- function(df) if (df == 0)
# #'         1 * (q >= 0)
# #'     else pchisq(q, df)
# #'     sum(sapply(df, pc) * mix)
# #'   }
# #'
#    pchisq.mix <- function( q, df , mix = rep(1,length(df))/length(df) ) {
#          # returns cdf for mixture of chi-squares. Usefule for testing
#          # random effects in mixed models
#            pc <- function( df ) if ( df == 0 ) 1*(q >= 0) else pchisq( q, df )
#            sum ( sapply( df, pc) * mix)
#    }
#
#
#
#
#
#
# #' as.character
# #'
# #'
# #'
# #'
# #'
# #' @param x
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x)
# #' as.character(x)
# #'
#    ch <- as.character
#
#
#
# #' Find the first non-missing element of a vector and check whether there are
# #' other inconsistent non-missing values
# #'
# #' Find the first non-missing element of a vector and check whether there are
# #' other inconsistent non-missing values
# #' description of what the function does. ~~
# #'
# #'
# #'
# #' @param x
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, ...)
# #' {
# #'     ret <- unique(x[!is.na(x)])
# #'     if (length(ret) > 1 && is.character(ret) && ("" %in% ret))
# #'         ret <- ret[ret != ""]
# #'     if (length(ret) > 1) {
# #'         cat("\n\n============= Multiple values in select.first ==========\n")
# #'         for (i in 1:length(ret)) cat("\nValue", i, ": <<", ret[i],
# #'             ">>\n")
# #'     }
# #'     ret[1]
# #'   }
# #'
# #' @export
#    select.first <- function( x , ... ) {
#        ret <- unique( x [ !is.na(x) ])
#        if ( length(ret) > 1 && is.character( ret ) && ( "" %in% ret)) ret <- ret [ ret != ""]
#        if ( length(ret) > 1 ) {
#         cat("\n\n============= Multiple values in select.first ==========\n")
#         for ( i in 1: length(ret)) cat("\nValue",i,": <<",ret[i],">>\n")
#        }
#        ret [1]
#    }
#
# #' @export
#    fill <- function( x, ... ) UseMethod("fill")
#
#
#
# #' Create time invariant variable by filling in missing values
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param by
# #' @param FUN
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, by, FUN = select.first, ...)
# #' {
# #'     levs <- levels(x)
# #'     ret <- capply(ch(x), by, FUN, ...)
# #'     factor(ret, levels = levs)
# #'   }
# #'
# #' @export
#    fill.factor <- function( x, by, FUN = select.first, ...) {
#       levs <- levels(x)
#       ret <- capply( ch( x), by, FUN , ... )
#       factor( ret, levels = levs)
#    }
#
#
#
# #' Create a time-invariant variable by filling in NAs
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param by
# #' @param FUN
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, by, FUN = select.first, ...)
# #' {
# #'     as.Date(capply(ch(x), by, FUN, ...))
# #'   }
# #'
# #' @export
#    fill.Date <- function( x, by, FUN = select.first, ...) {
#       as.Date( capply(ch(x), by, FUN, ...))
#    }
#
#
#
# #' Default method to fill in values of time-invariant variable
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param by
# #' @param FUN
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, by, FUN = select.first, ...)
# #' {
# #'     capply(x, by, FUN, ...)
# #'   }
# #'
# #' @export
#    fill.default <- function( x , by , FUN = select.first, ...) {
#       capply( x, by, FUN, ...)
#    }
#
#
#
# #' Fill in missing values to create a time-invariant variable
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param by
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, by, ...)
# #' {
# #'     ret <- list()
# #'     for (nn in names(x)) {
# #'         ret[[nn]] <- fill(x[[nn]], by)
# #'     }
# #'     as.data.frame(ret)
# #'   }
# #'
# #' @export
#    fill.data.frame <- function ( x, by ,... ) {
#           ret <- list()
#           for ( nn in names(x)) {
#               ret[[nn]] <- fill( x[[nn]], by)
#           }
#           as.data.frame(ret)
#    }
#
#
#     ##
#     ##   cvar: V0.1 August 15, 2006
#     ##   Creating contextual variables for categorical variables
#     ##
#     ##   cvar is designed to create contextual variables
#     ##   for factors, as well as for numerical variables.
#     ##   If a factor has g levels, convar will create a
#     ##   matrix with g-1 columns each of which is the within group
#     ##   mean value of the correponding column of the "contrast"
#     ##   matrix for the factor.
#     ##
#
#
# #' @export
#     cvar <- function( x, id ,... ) {
#         help = "
#         cvar: creates contextual group mean variables within levels of 'id'.\n
#               If 'x' is a factor, 'cvar' returns a matrix labelled so that it\n
#               is consistent with labels generated for coding variables for 'x'.\n
#               Example:\n
#                \n
#                 dd <- data.frame(x= 1:100, id = rep( LETTERS[1:10], each = 10))\n
#                 dd$a <- factor(sample( c('a','b','c'), 100, replace = T))\n
#                 dd$y <- dd$x + rep(rnorm(10), each = 10) + rnorm(100) + as.numeric(dd$a)\n
#                 library(nlme)\n
#                 fit <- lme( y ~ x + cvar(x,id), dd, random = ~ 1 + dvar(x,id) | id)\n
#                 anova( fit , type = 'm')\n
#                                         \n
#               The output of 'anova' can be used to test whether a contextual effect\n
#               needs to be included in the model.\n
#                                                 \n
#               See also: dvar for group-mean centering: x - cvar(x, id)\n
#         "
#         UseMethod("cvar")
#     }
#
# #' @export
#     cvar.factor <- function( x, id, ... ) {
#         mat <- contrasts( x) [ x,]
#         ret <- cvar(mat, id, ...)
#         colnames(ret) <- colnames(mat)
#         ret
#     }
#
# #' @export
#     cvar.default <- function( x, id, ... ) {
#         if ( is.matrix (x) ) {
#             if ( dim(x)[2] == 1) return( cvar( x[,], id, ...))
#             else {
#                 ret <-  cbind( cvar(x[,1], id, ...), cvar(x[,-1],id,...))
#                 colnames(ret) <- colnames(x)
#                 return( ret )
#             }
#         } else {
#             capply( x, id, mean, na.rm = T)
#         }
#     }
#
#
# #' @export
#     dvar <- function( x, id ,... ) {
#         help = "
#         dvar: produces group mean centering: x - cvar(x, id)
#         See 'cvar'
#         "
#         UseMethod("dvar")
#     }
#
# #' @export
#     dvar.factor <- function( x, id, ... ) {
#         mat <- contrasts( x) [ x,]
#         ret <- mat - cvar(mat, id, ...)
#         colnames(ret) <- colnames(mat)
#         ret
#     }
#
# #' @export
#     dvar.default <- function( x, id, ... ) {
#         if ( is.matrix (x) ) {
#             if ( dim(x)[2] == 1) return( dvar( x[,], id, ...))
#             else {
#                 ret <-  cbind( dvar(x[,1], id, ...), dvar(x[,-1],id,...))
#                 colnames(ret) <- colnames(x)
#                 return( ret )
#             }
#         } else {
#             x - capply( x, id, mean, na.rm = T)
#         }
#     }
#
# ##
# ##  sum
# ##
#
# #' @export
# cvars <- function(  x, by, ...) {
#       if ( length(x) == 1 && x == 1) {
#             n <- nrow(as.data.frame(by))
#             capply( rep(1,n), by, sum)
#       } else {
#             capply( x, by, sum, ...)
#       }
# }
#
#
#
#
#
# #' Transform NAs to 0
# #'
# #'
# #'
# #'
# #'
# #' @param x
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x)
# #' {
# #'     x[is.na(x)] <- 0
# #'     x
# #'   }
# #'
# #' @export
# na20 <- function(x) {
#      x[is.na(x)] <- 0
#      x
# }
#
#
#
# #' Describe a vector
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param descript
# #' @param exclude.missing
# #' @param digits
# #' @param weights
# #' @param normwt
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, descript, exclude.missing = TRUE, digits = 4, weights = NULL,
# #'     normwt = FALSE, ...)
# #' {
# #'     oldopt <- options(digits = digits)
# #'     on.exit(options(oldopt))
# #'     if (length(weights) == 0)
# #'         weights <- rep(1, length(x))
# #'     special.codes <- attr(x, "special.miss")$codes
# #'     labx <- attr(x, "label")
# #'     if (missing(descript))
# #'         descript <- as.character(sys.call())[2]
# #'     if (length(labx) && labx != descript)
# #'         descript <- paste(descript, ":", labx)
# #'     un <- attr(x, "units")
# #'     if (length(un) && un == "")
# #'         un <- NULL
# #'     fmt <- attr(x, "format")
# #'     if (length(fmt) && (is.function(fmt) || fmt == ""))
# #'         fmt <- NULL
# #'     if (length(fmt) > 1)
# #'         fmt <- paste(as.character(fmt[[1]]), as.character(fmt[[2]]))
# #'     present <- if (all(is.na(x)))
# #'         rep(FALSE, length(x))
# #'     else if (is.character(x))
# #'         (if (.R.)
# #'             x != "" & x != " " & !is.na(x)
# #'         else x != "" & x != " ")
# #'     else !is.na(x)
# #'     present <- present & !is.na(weights)
# #'     if (length(weights) != length(x))
# #'         stop("length of weights must equal length of x")
# #'     if (normwt) {
# #'         weights <- sum(present) * weights/sum(weights[present])
# #'         n <- sum(present)
# #'     }
# #'     else n <- round(sum(weights[present]), 2)
# #'     if (exclude.missing && n == 0)
# #'         return(structure(NULL, class = "describe"))
# #'     missing <- round(sum(weights[!present], na.rm = TRUE), 2)
# #'     atx <- attributes(x)
# #'     atx$names <- atx$dimnames <- atx$dim <- atx$special.miss <- NULL
# #'     atx$class <- atx$class[atx$class != "special.miss"]
# #'     isdot <- testDateTime(x, "either")
# #'     isdat <- testDateTime(x, "both")
# #'     x <- x[present, drop = FALSE]
# #'     x.unique <- sort(unique(x))
# #'     weights <- weights[present]
# #'     n.unique <- length(x.unique)
# #'     attributes(x) <- attributes(x.unique) <- atx
# #'     isnum <- (is.numeric(x) || isdat) && !is.category(x)
# #'     timeUsed <- isdat && testDateTime(x.unique, "timeVaries")
# #'     z <- list(descript = descript, units = un, format = fmt)
# #'     counts <- c(n, missing)
# #'     lab <- c("n", "missing")
# #'     if (length(special.codes)) {
# #'         tabsc <- table(special.codes)
# #'         counts <- c(counts, tabsc)
# #'         lab <- c(lab, names(tabsc))
# #'     }
# #'     if (length(atx$imputed)) {
# #'         counts <- c(counts, length(atx$imputed))
# #'         lab <- c(lab, "imputed")
# #'     }
# #'     if (length(pd <- atx$partial.date)) {
# #'         if ((nn <- length(pd$month)) > 0) {
# #'             counts <- c(counts, nn)
# #'             lab <- c(lab, "missing month")
# #'         }
# #'         if ((nn <- length(pd$day)) > 0) {
# #'             counts <- c(counts, nn)
# #'             lab <- c(lab, "missing day")
# #'         }
# #'         if ((nn <- length(pd$both)) > 0) {
# #'             counts <- c(counts, nn)
# #'             lab <- c(lab, "missing month,day")
# #'         }
# #'     }
# #'     if (length(atx$substi.source)) {
# #'         tabss <- table(atx$substi.source)
# #'         counts <- c(counts, tabss)
# #'         lab <- c(lab, names(tabss))
# #'     }
# #'     counts <- c(counts, n.unique)
# #'     lab <- c(lab, "unique")
# #'     x.binary <- n.unique == 2 && isnum && x.unique[1] == 0 &&
# #'         x.unique[2] == 1
# #'     if (x.binary) {
# #'         counts <- c(counts, sum(weights[x == 1]))
# #'         lab <- c(lab, "Sum")
# #'     }
# #'     if (isnum) {
# #'         xnum <- if (.SV4.)
# #'             as.numeric(x)
# #'         else oldUnclass(x)
# #'         if (isdot) {
# #'             dd <- sum(weights * xnum)/sum(weights)
# #'             fval <- formatDateTime(dd, atx, !timeUsed)
# #'             counts <- c(counts, fval)
# #'         }
# #'         else counts <- c(counts, format(sum(weights * x)/sum(weights),
# #'             ...))
# #'         lab <- c(lab, "Mean")
# #'     }
# #'     if (n.unique >= 10 & isnum) {
# #'         q <- if (any(weights != 1))
# #'             wtd.quantile(xnum, weights, normwt = FALSE, na.rm = FALSE,
# #'                 probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
# #'         else quantile(xnum, c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9,
# #'             0.95), na.rm = FALSE)
# #'         fval <- if (isdot)
# #'             formatDateTime(q, atx, !timeUsed)
# #'         else format(q, ...)
# #'         counts <- c(counts, fval)
# #'         lab <- c(lab, ".05", ".10", ".25", ".50", ".75", ".90",
# #'             ".95")
# #'     }
# #'     names(counts) <- lab
# #'     z$counts <- counts
# #'     counts <- NULL
# #'     if (n.unique >= 20) {
# #'         if (isnum) {
# #'             r <- range(xnum)
# #'             xg <- pmin(1 + floor((100 * (xnum - r[1]))/(r[2] -
# #'                 r[1])), 100)
# #'             z$intervalFreq <- list(range = as.single(r), count = as.integer(tabulate(xg)))
# #'         }
# #'         lo <- x.unique[1:5]
# #'         hi <- x.unique[(n.unique - 4):n.unique]
# #'         fval <- if (isdot)
# #'             formatDateTime(c(oldUnclass(lo), oldUnclass(hi)),
# #'                 atx, !timeUsed)
# #'         else format(c(format(lo), format(hi)), ...)
# #'         counts <- fval
# #'         names(counts) <- c("L1", "L2", "L3", "L4", "L5", "H5",
# #'             "H4", "H3", "H2", "H1")
# #'     }
# #'     if (n.unique > 1 && n.unique < 20 && !x.binary) {
# #'         tab <- wtd.table(if (isnum)
# #'             format(x)
# #'         else x, weights, normwt = FALSE, na.rm = FALSE, type = "table")
# #'         pct <- round(100 * tab/sum(tab))
# #'         counts <- t(as.matrix(tab))
# #'         counts <- rbind(counts, pct)
# #'         dimnames(counts)[[1]] <- c("Frequency", "%")
# #'     }
# #'     z$values <- counts
# #'     structure(z, class = "describe")
# #'   }
# #'
# describe.vector <-
# function (x, descript, exclude.missing = TRUE, digits = 4, weights = NULL,
#     normwt = FALSE, ...)
# {
#     # GM: modified by rounding n and missing
#     oldopt <- options(digits = digits)
#     on.exit(options(oldopt))
#     if (length(weights) == 0)
#         weights <- rep(1, length(x))
#     special.codes <- attr(x, "special.miss")$codes
#     labx <- attr(x, "label")
#     if (missing(descript))
#         descript <- as.character(sys.call())[2]
#     if (length(labx) && labx != descript)
#         descript <- paste(descript, ":", labx)
#     un <- attr(x, "units")
#     if (length(un) && un == "")
#         un <- NULL
#     fmt <- attr(x, "format")
#     if (length(fmt) && (is.function(fmt) || fmt == ""))
#         fmt <- NULL
#     if (length(fmt) > 1)
#         fmt <- paste(as.character(fmt[[1]]), as.character(fmt[[2]]))
#     present <- if (all(is.na(x)))
#         rep(FALSE, length(x))
#     else if (is.character(x))
#         (if (.R.)
#             x != "" & x != " " & !is.na(x)
#         else x != "" & x != " ")
#     else !is.na(x)
#     present <- present & !is.na(weights)
#     if (length(weights) != length(x))
#         stop("length of weights must equal length of x")
#     if (normwt) {
#         weights <- sum(present) * weights/sum(weights[present])
#         n <- sum(present)
#     }
#     else n <- round(sum(weights[present]),2)
#     if (exclude.missing && n == 0)
#         return(structure(NULL, class = "describe"))
#     missing <- round(sum(weights[!present], na.rm = TRUE),2)
#     atx <- attributes(x)
#     atx$names <- atx$dimnames <- atx$dim <- atx$special.miss <- NULL
#     atx$class <- atx$class[atx$class != "special.miss"]
#     isdot <- testDateTime(x, "either")
#     isdat <- testDateTime(x, "both")
#     x <- x[present, drop = FALSE]
#     x.unique <- sort(unique(x))
#     weights <- weights[present]
#     n.unique <- length(x.unique)
#     attributes(x) <- attributes(x.unique) <- atx
#     isnum <- (is.numeric(x) || isdat) && !is.category(x)
#     timeUsed <- isdat && testDateTime(x.unique, "timeVaries")
#     z <- list(descript = descript, units = un, format = fmt)
#     counts <- c(n, missing)
#     lab <- c("n", "missing")
#     if (length(special.codes)) {
#         tabsc <- table(special.codes)
#         counts <- c(counts, tabsc)
#         lab <- c(lab, names(tabsc))
#     }
#     if (length(atx$imputed)) {
#         counts <- c(counts, length(atx$imputed))
#         lab <- c(lab, "imputed")
#     }
#     if (length(pd <- atx$partial.date)) {
#         if ((nn <- length(pd$month)) > 0) {
#             counts <- c(counts, nn)
#             lab <- c(lab, "missing month")
#         }
#         if ((nn <- length(pd$day)) > 0) {
#             counts <- c(counts, nn)
#             lab <- c(lab, "missing day")
#         }
#         if ((nn <- length(pd$both)) > 0) {
#             counts <- c(counts, nn)
#             lab <- c(lab, "missing month,day")
#         }
#     }
#     if (length(atx$substi.source)) {
#         tabss <- table(atx$substi.source)
#         counts <- c(counts, tabss)
#         lab <- c(lab, names(tabss))
#     }
#     counts <- c(counts, n.unique)
#     lab <- c(lab, "unique")
#     x.binary <- n.unique == 2 && isnum && x.unique[1] == 0 &&
#         x.unique[2] == 1
#     if (x.binary) {
#         counts <- c(counts, sum(weights[x == 1]))
#         lab <- c(lab, "Sum")
#     }
#     if (isnum) {
#         xnum <- if (.SV4.)
#             as.numeric(x)
#         else oldUnclass(x)
#         if (isdot) {
#             dd <- sum(weights * xnum)/sum(weights)
#             fval <- formatDateTime(dd, atx, !timeUsed)
#             counts <- c(counts, fval)
#         }
#         else counts <- c(counts, format(sum(weights * x)/sum(weights),
#             ...))
#         lab <- c(lab, "Mean")
#     }
#     if (n.unique >= 10 & isnum) {
#         q <- if (any(weights != 1))
#             wtd.quantile(xnum, weights, normwt = FALSE, na.rm = FALSE,
#                 probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))
#         else quantile(xnum, c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9,
#             0.95), na.rm = FALSE)
#         fval <- if (isdot)
#             formatDateTime(q, atx, !timeUsed)
#         else format(q, ...)
#         counts <- c(counts, fval)
#         lab <- c(lab, ".05", ".10", ".25", ".50", ".75", ".90",
#             ".95")
#     }
#     names(counts) <- lab
#     z$counts <- counts
#     counts <- NULL
#     if (n.unique >= 20) {
#         if (isnum) {
#             r <- range(xnum)
#             xg <- pmin(1 + floor((100 * (xnum - r[1]))/(r[2] -
#                 r[1])), 100)
#             z$intervalFreq <- list(range = as.single(r), count = as.integer(tabulate(xg)))
#         }
#         lo <- x.unique[1:5]
#         hi <- x.unique[(n.unique - 4):n.unique]
#         fval <- if (isdot)
#             formatDateTime(c(oldUnclass(lo), oldUnclass(hi)),
#                 atx, !timeUsed)
#         else format(c(format(lo), format(hi)), ...)
#         counts <- fval
#         names(counts) <- c("L1", "L2", "L3", "L4", "L5", "H5",
#             "H4", "H3", "H2", "H1")
#     }
#     if (n.unique > 1 && n.unique < 20 && !x.binary) {
#         tab <- wtd.table(if (isnum)
#             format(x)
#         else x, weights, normwt = FALSE, na.rm = FALSE, type = "table")
#         pct <- round(100 * tab/sum(tab))
#         counts <- t(as.matrix(tab))
#         counts <- rbind(counts, pct)
#         dimnames(counts)[[1]] <- c("Frequency", "%")
#     }
#     z$values <- counts
#     structure(z, class = "describe")
# }
#
#
#
#
# #' Include NAs
# #'
# #'
# #'
# #'
# #'
# #' @param obj
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (obj)
# #' {
# #'     if (inherits(obj, "data.frame"))
# #'         for (i in seq(along = obj)) obj[[i]] <- na.include(obj[[i]])
# #'     else {
# #'         if (length(levels(obj)) && any(is.na(obj)))
# #'             obj <- factor(obj, exclude = NULL)
# #'     }
# #'     obj
# #'   }
# #'
# #' @export
# na.include  <- function (obj)
# {
#      # from library(Hmisc)
#     if (inherits(obj, "data.frame"))
#         for (i in seq(along = obj)) obj[[i]] <- na.include(obj[[i]])
#     else {
#         if (length(levels(obj)) && any(is.na(obj)))
#             obj <- factor(obj, exclude = NULL)
#     }
#     obj
# }
#
#
#
#
# #' Special version of summary
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param \dots
# #' @export
# summ <- function(x,...) UseMethod("summ")
# #' Summary for lmer objects
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param \dots
# #' @export
# summ.lmer <- function(x, ...) {
#                ret <- c(AIC = AIC(x@logLik), BIC= BIC(x@logLik), logLik=x@logLik)
#                ret
# }
#
#
#
# #' Alternative print -- generic
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, ...)
# #' {
# #'     UseMethod("pr")
# #'   }
# #'
# #' @export
# pr <- function(x,...) {
#    # print to cut and paste as input
#    UseMethod("pr")
# }
#
#
# #' Alternative print -- should probably use dput
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param pre
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, pre = "\t", ...)
# #' {
# #'     for (xx in x) cat(pre, "\"", xx, "\",\n", sep = "")
# #'     invisible(x)
# #'   }
# #'
# #' @export
# pr.default <- function(x,pre="\t",...) {
#           # cat('\nc(')
#           for ( xx in x) cat(pre,'"',xx,'",\n',sep='')
#           invisible(x)
# }
#
# ##  reorder.factor has been removed because the new method in
# ##  in the base package does the same thing
# ##reorder.factor <- function (x, v, FUN = mean, ...) {
# ##        warning("gm version produces ordinary factor -- not ordered factor")
# ##               factor(x, levels = levels(x)[order(tapply(v, x, FUN, ...))])
# ##}
# ## NOTE: reorder.factor in <base> is hidden. You must use is as reorder(x)
#
# ## reorder.default <- function(x,...) reorder( factor(x),...)  # so reorder works on non-factors
#
#
#
#
# #' Q matrix of QR decomposion with missing data
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param verbose
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, verbose = 0)
# #' {
# #'     miss <- apply(x, 1, function(xx) any(is.na(xx)))
# #'     xc <- x[!miss, ]
# #'     qf <- qr.Q(qqr <- qr(xc))
# #'     if (verbose > 0) {
# #'         cat("xc:", dim(xc), "\n")
# #'         cat("qf:", dim(qf), "\n")
# #'         print(qf)
# #'     }
# #'     if (ncol(xc) > ncol(qf))
# #'         xc <- xc[, 1:ncol(qf)]
# #'     ip <- sign(apply(qf * xc, 2, sum))
# #'     qf <- qf * rep(ip, rep(nrow(qf), length(ip)))
# #'     ret <- array(NA, dim = c(nrow(x), ncol(qf)))
# #'     rownames(ret) <- rownames(x)
# #'     colnames(ret) <- colnames(xc)
# #'     ret[!miss, ] <- qf
# #'     attr(ret, "rank") <- qqr$rank
# #'     attr(ret, "miss") <- miss
# #'     ret
# #'   }
# #'
# #' @export
# Q <- function(x, verbose = 0) {
#     # find the Q matrix of a qr decomposition with possible NAs
#     miss <- apply(x, 1, function(xx) any(is.na(xx)))
#     xc <- x[!miss,]
#     qf <- qr.Q(qqr<-qr(xc))
#
#     if(verbose > 0) {
#                cat("xc:", dim(xc),'\n')
#                cat("qf:", dim(qf), '\n')
#
#                print(qf)
#     }
#     if( ncol(xc) > ncol(qf)) xc <- xc[,1:ncol(qf)]
#     ip <- sign(apply(qf*xc,2,sum))   # sign of reln between Q and X
#     qf <- qf * rep(ip, rep( nrow(qf), length(ip)))   # change sign of each row
#     ret <- array(NA, dim = c(nrow(x),ncol(qf)))
#     rownames(ret) <- rownames(x)
#     colnames(ret) <- colnames(xc)
#     ret[!miss,] <- qf
#     attr(ret,'rank') <- qqr$rank
#     attr(ret,'miss') <- miss
#     ret
# }
#
#
#
#
#
# #' Generic function to extend 'contrasts' to 'lmer' objects.
# #'
# #' Generic function to extend 'contrasts' to 'lmer' objects.
# #'
# #'
# #'
# #' @param x
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x)
# #' UseMethod("Contrasts")
# #'
# #' @export
# Contrasts <- function(x) UseMethod("Contrasts")
#
# #' @export
# Contrasts.default <- function(x) contrasts(x)
#
# #' @export
# Contrasts.lmer <- function(x) {
#        dd <- x@frame[1,]
#        ret <- list()
#        for ( nn in names(dd)) {
#            vv <- dd[[nn]]
#            if (is.factor(vv)) ret[[nn]] <- contrasts(vv)
#        }
#        vv
# }
# #' Older version of comp
# #'
# #'
# #'
# #'
# #'
# #' @param fit
# #' @param form
# #' @param varname
# #' @param varpattern
# #' @param data
# #' @param \dots
# #' @export
# comp.old <- function(fit, form, varname, varpattern = vname, data = getElement(fit,'frame'), ...) {
#      ## Computing regression components
#      ## currently works form lmer only
#      ## idea is to return a data frame with value a variable and corresponding fitted values
#      ## for the 'component' of that variable or variables
#      ##
#      ## This can be used to fit a model to a prediction data.frame with only
#      ## the necessary variables defined!!
#      ## However BEWARE not to use 'unsafe' transformations (e.g. Q, poly)
#      ##
#      model.mat <- model.matrix(form, data,...)
#      #print(dim(model.mat))
#      ret <- data[rownames(model.mat),varname,drop = F]
#      fe <- fixef(fit)
#      pos.fe <- grep( varpattern, names(fe))
#      pos.mat <- grep( varpattern, colnames(model.mat))
#      ret$comp <- c( model.mat[,pos.mat] %*% fe[pos.fe] )
#      attr(ret,"predictors") <-  names(fe)[pos.fe]
#      ret
# }
#
#
#
#
#
#
#
# #' Backup of former version of comp
# #'
# #'
# #'
# #'
# #'
# #' @param fit
# #' @param form
# #' @param varname
# #' @param varpattern
# #' @param data
# #' @param \dots
# #' @export
# comp.bak <- function(fit, form, varname, varpattern = vname, data = getElement(fit,'frame'), ...) {
#      ## Computing regression components
#      ## currently works form lmer only
#      ## idea is to return a data frame with value a variable and corresponding fitted values
#      ## for the 'component' of that variable or variables
#      ##
#      ## This can be used to fit a model to a prediction data.frame with only
#      ## the necessary variables defined!!
#      ## However BEWARE not to use 'unsafe' transformations (e.g. Q, poly)
#      ##
#      ## Arguments:
#      ##
#      ## fit : model for which component to be estimated
#      ##
#      ## form : portion of model formula to generate model variables needed for fitting
#      ##
#      ## varname: variable names to be included in output data frame
#      ##
#      ## varpattern:  regular expression so select effects in component
#      ##
#      ## data:  the 'parent' data.frame -- often a prediction data frame.
#      ##        When using the original data frame, it is often necessary to
#      ##        'na.action = na.omit' for '...'
#      ##
#      model.mat <- model.matrix(form, data,...)
#      #print(dim(model.mat))
#      ret <- data[rownames(model.mat),varname,drop = F]
#      fe <- fixef(fit)
#      effect.names <- grep( varpattern, names(fe), value = T)
#      ret$comp <- c( model.mat[,effect.names] %*% fe[effect.names] )
#      attr(ret,"predictors") <-  effect.names
#      ret
# }
#
#
#
#
#
#
# #' Prediction for 'lmer' objects
# #'
# #'
# #'
# #'
# #'
# #' @param fit
# #' @param form
# #' @param varname
# #' @param varpattern
# #' @param data
# #' @param \dots
# #' @export
# comp <- function(fit, form = terms(fit), varname, varpattern = "", data = getData(fit), ...) {
# # this is the version that worked for RDC code
# # the version only returns the adjusted predicted value
# # The reason for returning a portion of the original data frame was to ensure
# # correspondence between components and data values.
# # With 'na.comp' the components are padded to match the original data.
# # and can do so in a form that conforms to the original data mat
# ## for lmer: data = getElement(fit,'frame')
#      ## Computing regression components
#      ## currently works for lmer only
#      ## idea is to return a data frame with value a variable and corresponding fitted values
#      ## for the 'component' of that variable or variables
#      ##
#      ## This can be used to fit a model to a prediction data.frame with only
#      ## the necessary variables defined!!
#      ## However BEWARE not to use 'unsafe' transformations (e.g. Q, poly)
#      ##
#      ## Arguments:
#      ##
#      ## fit : model for which component to be estimated
#      ##
#      ## form : portion of model formula to generate model variables needed for fitting
#      ##
#      ## varname: variable names to be included in output data frame
#      ##
#      ## varpattern:  regular expression to select effects in component
#      ##
#      ## data:  the 'parent' data.frame -- often a prediction data frame.
#      ##        When using the original data frame, it is often necessary to
#      ##        'na.action = na.omit' for '...'
#      ##
#
# #' @export
#      getData <- function(x,...) UseMethod("getData")
# #     getData.lme <- function(x) x$data
# #' @export
#      getData.lme <- function(x,...) nlme:::getData(x,...)
# #' @export
#      getData.lmer <- function(x) getElement(x,'frame')
#      disp(form)
#      disp( dim(data))
#      model.mat <- model.matrix(form, data,...)
#      disp(dim(model.mat))
#      ret <- data[rownames(model.mat),,drop = F]
#      fe <- fixef(fit)
#      effect.names <- grep( varpattern, names(fe), value = TRUE)
#      ret$comp <- c( model.mat[,effect.names] %*% fe[effect.names] )
#      attr(ret,"predictors") <-  effect.names
#      ret
# }
#
#
#
#
#
#
#
# #' Prediction with 'lmer' objects
# #'
# #'
# #'
# #'
# #'
# #' @param fit
# #' @param varpattern
# #' @param form
# #' @param data
# #' @param \dots
# #' @export
#     com <- function(fit, varpattern = "", form = terms(fit), data = getData(fit), ...) {
#     # this is a new version of 'comp' that
#     # only returns the adjusted predicted value
#     # The reason for returning a portion of the original data frame was to ensure
#     # correspondence between components and data values.
#     # With 'na.com' the components are padded to match the original data.
#     # and can do so in a form that conforms to the original data mat
#     ## for lmer: data = fit@frame
#          ## Computing regression components
#          ## currently works for lmer only
#          ## idea is to return a data frame with value a variable and corresponding fitted values
#          ## for the 'component' of that variable or variables
#          ##
#          ## This can be used to fit a model to a prediction data.frame with only
#          ## the necessary variables defined!!
#          ## However BEWARE not to use 'unsafe' transformations (e.g. Q, poly)
#          ##
#          ## Arguments:
#          ##
#          ## fit : model for which component to be estimated
#          ##
#          ## form : portion of model formula to generate model variables needed for fitting
#          ##
#          ## varname: variable names to be included in output data frame
#          ##
#          ## varpattern:  regular expression to select effects in component
#          ##
#          ## data:  the 'parent' data.frame -- often a prediction data frame.
#          ##        When using the original data frame, it is often necessary to
#          ##        'na.action = na.omit' for '...'
#          ##
#
#          getData <- function(x) UseMethod("getData")
#          getData.lme <- function(x) x$data
#          getData.lmer <- function(x) getElement(x,'frame')
#          #disp(form)
#          #disp( dim(data))
#          model.mat <- model.matrix(form, data,...)
#          #disp(dim(model.mat))
#          # ret <- data[rownames(model.mat),,drop = F]
#          fe <- fixef(fit)
#          effect.names <- grep( varpattern, names(fe), value = TRUE)
#          ret <- c( model.mat[,effect.names] %*% fe[effect.names] )
#          names(ret) <- rownames(model.mat)
#          # return only components corresponding to fit$resid
#          retnames <- rownames(fit$resid)   # this main only work for 'lme'
#          ret <- ret[retnames]
#          attr(ret,"predictors") <-  effect.names
#          ret
#     }
#
#
#
# #' Apply com with NAs
# #'
# #' @param fit
# #' @param \dots
# #' @export
# na.com <- function( fit, ...) {
#         ret <- com(fit,...)
#         prs <- attributes(ret)$predictors
#         ret <- na.pad( fit, ret)
#         attr(ret,"predictors") <- prs
#         ret
#     }
#
#
#
# #' Compute residuals from component with NAs
# #'
# #' @param fit
# #' @param varp
# #' @param level
# #' @param \dots
# #' @export
# na.comres <- function( fit, varp = "", level = 0, ...) na.com( fit, varp = varp, ...) + na.resid(fit, level = level)
# #' Response with NAs
# #'
# #' @param fit
# #' @export
# na.getResponse <- function( fit) na.pad( fit, getResponse(fit))
#
#     #length(names(na.comp(fit.moda.age,varp='Month')))
#     #dim(dd)
#
#
# #  glh( fit, Lall(fit, 'Wave'))
#
#
#
#
#
# #' Extended svd
# #'
# #' @param x
# #' @export
# cond <- function(x) {
#      # reporting on conditioning of matrix
#      Svd <- svd(x)
#      ret <- list(svd = Svd, dim = dim(x), rank = qr(x)$rank, condition = Svd$d[1]/Svd$d[length(Svd$d)])
#      #if(nrow(x) == ncol(x)) ret <- c(ret, cor= list(cor(x)))
#      ret
# }
# #' Version of Vcov used in RDC
# #'
# #' @param fit
# #' @param L
# #' @export
# Vcov.rdc <- function( fit, L  = Lmat(fit,"") ) {
#      # variance of etahat = L beta.hat
#      vc <- getFix(fit)$vcov
#      L %*% vc %*% t(L)
# }
# #' Version of Vcor used in RDC
# #'
# #' @param fit
# #' @param L
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (fit, L = Lmat(fit, ""))
# #' {
# #'     ret <- cov2cor(vc <- Vcov(fit, L))
# #'     sv <- svd(vc)$d
# #'     attr(ret, "condition") <- sv[1]/sv[length(sv)]
# #'     attr(ret, "class") <- "correl"
# #'     ret
# #'   }
# #'
# #' @export
# Vcor.rdc <- function(fit, L = Lmat(fit,"")) {
#      ret <- cov2cor(vc <- Vcov(fit, L))
#      sv <- svd(vc)$d
#
#      attr(ret,'condition') <- sv[1]/sv[length(sv)]
#      attr(ret,'class') <- "correl"
#      ret
# }
#
#
#
#
# #' print correlations
# #'
# #'
# #'
# #'
# #'
# #' @param x
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x)
# #' {
# #'     correl <- format(round(x, 3), nsmall = 3, digits = 4)
# #'     correl[!lower.tri(correl)] <- ""
# #'     if (!is.null(cond <- attr(correl, "condition"))) {
# #'         attr(correl, "condition") <- NULL
# #'     }
# #'     print(correl, quote = FALSE)
# #'     if (!is.null(cond))
# #'         cat("Condition:", cond, "\n")
# #'     invisible(x)
# #'   }
# #'
# #' @export
# print.correl <- function(x) {
#       correl <- format(round(x, 3), nsmall = 3,
#                   digits = 4)
#       correl[!lower.tri(correl)] <- ""
#       if (!is.null(cond <- attr(correl,"condition"))) {
#          attr(correl,"condition") <- NULL
#       }
#       print(correl, quote = FALSE)
#       if(!is.null(cond)) cat("Condition:",cond,"\n")
#       invisible(x)
# }
#
#
#
#
# #' Test version of glh in RDC
# #'
# #'
# #'
# #'
# #'
# #' @param fit
# #' @param Llist
# #' @param help
# #' @param clevel
# #' @param debug
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (fit, Llist, help = F, clevel = 0.95, debug = F)
# #' {
# #'     if (help) {
# #'         cat("help!!!")
# #'         return(0)
# #'     }
# #'     if (!is.list(Llist))
# #'         Llist <- list(Llist)
# #'     ret <- list()
# #'     fix <- getFix(fit)
# #'     beta <- fix$fixed
# #'     vc <- fix$vcov
# #'     dfs <- fix$df
# #'     for (ii in 1:length(Llist)) {
# #'         ret[[ii]] <- list()
# #'         L <- rbind(zz <- Llist[[ii]])
# #'         heading <- attr(zz, "heading")
# #'         nam <- names(Llist)[ii]
# #'         qqr <- qr(t(L))
# #'         L.rank <- qqr$rank
# #'         L.full <- t(qr.Q(qqr))[1:L.rank, , drop = F]
# #'         if (F) {
# #'             cat("\n+++++++++++++++++++\n")
# #'             print(nam)
# #'             print(L)
# #'             print(svd(L)$d)
# #'             print(L.full)
# #'             print(svd(L.full)$d)
# #'         }
# #'         if (debug) {
# #'             print(L)
# #'             print(dim(L.full))
# #'             print(dim(vc))
# #'         }
# #'         vv <- L.full %*% vc %*% t(L.full)
# #'         eta.hat <- L.full %*% beta
# #'         Fstat <- (t(eta.hat) %*% qr.solve(vv, eta.hat))/L.rank
# #'         Fstat2 <- (t(eta.hat) %*% solve(vv) %*% eta.hat)/L.rank
# #'         included.effects <- apply(L, 2, function(x) sum(abs(x))) !=
# #'             0
# #'         denDF <- min(dfs[included.effects])
# #'         numDF <- L.rank
# #'         ret.anova <- rbind(c(numDF, denDF, Fstat, Fstat2, pf(Fstat,
# #'             numDF, denDF, lower.tail = F)))
# #'         colnames(ret.anova) <- c("numDF", "denDF", "F value",
# #'             "F2", "Pr(>F)")
# #'         rownames(ret.anova) <- nam
# #'         ret[[ii]]$anova <- ret.anova
# #'         etahat <- L %*% beta
# #'         etavar <- L %*% vc %*% t(L)
# #'         etasd <- sqrt(diag(etavar))
# #'         denDF <- apply(L, 1, function(x, dfs) min(dfs[x != 0]),
# #'             dfs = dfs)
# #'         aod <- cbind(c(etahat), etasd, denDF, c(etahat/etasd),
# #'             2 * pt(-abs(etahat/etasd), denDF))
# #'         colnames(aod) <- c("Estimate", "Std.Error", "DF", "t value",
# #'             "Pr(>|t|)")
# #'         if (!is.null(clevel)) {
# #'             hw <- qt(1 - (1 - clevel)/2, denDF) * etasd
# #'             aod <- cbind(aod, LL = etahat - hw, UL = etahat +
# #'                 hw)
# #'             labs <- paste(c("Lower", "Upper"), format(clevel))
# #'             colnames(aod)[ncol(aod) + c(-1, 0)] <- labs
# #'         }
# #'         rownames(aod) <- rownames(L)
# #'         ret[[ii]]$estimate <- aod
# #'         ret[[ii]]$vcov <- Vcov(fit, L)
# #'         ret[[ii]]$vcor <- Vcor(fit, L)
# #'         ret[[ii]]$L <- L
# #'         ret[[ii]]$L.full <- L.full
# #'         if (!is.null(heading))
# #'             attr(ret[[ii]], "heading") <- heading
# #'     }
# #'     names(ret) <- names(Llist)
# #'     attr(ret, "class") <- "glh"
# #'     ret
# #'   }
# #'
# #' @export
# glh.rdc <- function(fit, Llist, help = FALSE, clevel = 0.95, debug = FALSE) {
#     if(help) {
#          cat("help!!!")
#          return(0)
#     }
#     if( !is.list(Llist) ) Llist <- list(Llist)
#     ret <- list()
#     fix <- getFix(fit)
#     beta <- fix$fixed
#     vc <- fix$vcov
#     dfs <- fix$df
#     for ( ii in 1:length(Llist)) {
#         ret[[ii]] <- list()
#
#         L <- rbind(zz <-  Llist[[ii]])
#         heading <- attr(zz, "heading")
#         nam <- names(Llist)[ii]
#
#         ## Anova
#
#         qqr <- qr(t(L))
#        #if(debug) print(qqr)
#         L.rank <- qqr$rank
#         L.full <- t(qr.Q(qqr)) [ 1:L.rank,,drop=F]
#         if(F) {
#                   cat("\n+++++++++++++++++++\n")
#                   print(nam)
#                   print(L)
#                   print(svd(L)$d)
#                   print(L.full)
#                   print(svd(L.full)$d)
#          }
#          if (debug) {
#             print(L)
#             print(dim(L.full))
#             print(dim(vc))
#          }
#         vv <- L.full %*% vc %*% t(L.full)
#         eta.hat <- L.full %*% beta
#         Fstat <- (t(eta.hat) %*% qr.solve(vv, eta.hat))/L.rank
#         Fstat2 <- (t(eta.hat) %*% solve(vv) %*% eta.hat)/L.rank
#         included.effects <- apply(L, 2, function(x) sum(abs(x))) != 0
#         denDF <- min(dfs[included.effects])
#         numDF <- L.rank
#
#         ret.anova <- rbind(c(numDF, denDF, Fstat,Fstat2, pf(Fstat, numDF, denDF, lower.tail = FALSE)))
#         colnames(ret.anova) <- c("numDF","denDF","F value","F2","Pr(>F)")
#         rownames(ret.anova) <-  nam
#         ret[[ii]]$anova <- ret.anova
#
#         ## Estimate
#
#         etahat <- L %*% beta
#         etavar <- L %*% vc %*% t(L)
#         etasd <- sqrt(diag(etavar))
#         denDF <- apply( L , 1, function(x,dfs) min(dfs[x!=0]), dfs = dfs)
#         aod <- cbind(
#                  c(etahat),
#                  etasd,
#                  denDF,
#                  c(etahat/etasd),
#                  2*pt(-abs(etahat/etasd), denDF))
#         colnames(aod) <- c("Estimate","Std.Error",'DF','t value','Pr(>|t|)')
#         if( !is.null(clevel) ) {
#             hw <- qt( 1-(1-clevel)/2, denDF) * etasd
#             aod <- cbind(aod, LL = etahat - hw, UL = etahat + hw)
#
#
# #' Add labels -- generic
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, ...)
# #' UseMethod("labs")
# #'
#             labs <- paste( c("Lower","Upper"), format(clevel))
#             colnames(aod)[ ncol(aod) + c(-1,0)] <- labs
#
#         }
#         #aod <- as.data.frame(aod)
#         rownames(aod) <- rownames(L)
#         ret[[ii]]$estimate <- aod
#
#         ## Vcov
#
#         ret[[ii]]$vcov <- Vcov( fit, L)
#         ret[[ii]]$vcor <- Vcor(fit,L)
#         ret[[ii]]$L <- L
#         ret[[ii]]$L.full <- L.full
#         if ( !is.null(heading)) attr(ret[[ii]],"heading") <- heading
#     }
#     names(ret) <- names(Llist)
#     attr(ret,"class") <- "glh"
#     ret
# }
#
#
#
#
#
#
#
# #' Format a coefficient table
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param digits
# #' @param pdigits
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, digits = 6, pdigits = digits - 1)
# #' {
# #'     pformat <- function(x, digits) {
# #'         x <- format(xx <- round(x, digits))
# #'         x[as.double(xx) == 0] <- paste(c("<.", rep("0", digits -
# #'             1), "1"), collapse = "")
# #'         x
# #'     }
# #'     xx <- array("", dim = dim(x), dimnames = dimnames(x))
# #'     for (i in 1:ncol(xx)) {
# #'         xx[, i] <- format(round(x[, i], digits), digits = digits)
# #'     }
# #'     if (length(isp <- grep("^Pr\(", colnames(x))) > 0) {
# #'         xx[, isp[1]] <- pformat(x[, isp[1]], digits = pdigits)
# #'     }
# #'     xx
# #'   }
# #'
# #' @export
# formatCoefmat <- function(x ,digits = 6, pdigits = digits-1 ) {
#      pformat <- function(x, digits) {
#              x <- format(xx <- round(x,digits))
#              x[ as.double(xx) == 0 ] <- paste(c("<.",
#                               rep('0',digits-1),"1"), collapse = "")
#              x
#      }
#      xx <- array("",dim=dim(x), dimnames = dimnames(x))
#      for ( i in 1:ncol(xx)) {
#          xx[,i] <- format(round(x[,i],digits),digits = digits)
#      }
#      if ( length( isp <- grep("^Pr\\(",colnames(x))) > 0) {
#          xx[,isp[1]] <- pformat( x[,isp[1]], digits = pdigits)
#      }
#      xx
# }
#
#
#
# #' Print a 'glh' object tested in RDC
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param round
# #' @param pround
# #' @param L
# #' @param cov
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, round = 6, pround = round - 1, L = T, cov = T, ...)
# #' {
# #'     rnd <- function(x, digits) {
# #'         if (is.numeric(x))
# #'             x <- round(x, digits = digits)
# #'         format(x)
# #'     }
# #'     for (ii in 1:length(x)) {
# #'         nn <- names(x)[[ii]]
# #'         tt <- x[[ii]]
# #'         ta <- tt$anova
# #'         tap <- array("", dim = dim(ta), dimnames = dimnames(ta))
# #'         cat("\n", nn, "\n", sep = "")
# #'         print(formatCoefmat(ta, digits = round, pdigits = pround),
# #'             quote = F, right = T)
# #'         cat("\n")
# #'         te <- tt$estimate
# #'         if (!is.null(zhead <- attr(tt, "heading")))
# #'             cat(zhead, "\n")
# #'         print(formatCoefmat(te, digits = round, pdigits = pround),
# #'             quote = F, right = T)
# #'         if (L == T) {
# #'             cat("\nL:\n")
# #'             print(tt$L)
# #'             if (dim(tt$L.full)[1] < dim(tt$L)[1]) {
# #'                 cat("\nL (full rank):\n")
# #'                 print(tt$L.full)
# #'             }
# #'         }
# #'         if (cov == T) {
# #'             cat("\nVar-Cov of estimates:\n")
# #'             print(tt$vcov)
# #'             cat("\nCorrelations:\n")
# #'             print(tt$vcor)
# #'         }
# #'     }
# #'     invisible(x)
# #'   }
# #'
# #' @export
# print.glh.rdc <- function(x, round = 6, pround = round - 1, L  = TRUE, cov = TRUE, ...) {
#      # should round by SD, i.e. keep 3 sig digits for sd and others rounded accordingly
#
#
#      rnd <- function(x, digits) {
#              if ( is.numeric( x)) x <- round(x, digits = digits)
#              format(x)
#      }
#      for ( ii in 1:length(x)) {
#          nn <- names(x)[[ii]]
#          tt <- x[[ii]]
#          ta <- tt$anova
#          tap <- array("", dim = dim(ta), dimnames = dimnames(ta))
#
#          # ta[["p-value"]] <- pformat( ta[["p-value"]], digits = pround)
#          ## print(as.data.frame(ta,row.names = nn))
#          cat("\n",nn,"\n",sep='')
#
#          print(formatCoefmat(ta,digits=round,pdigits=pround),quote=F,right=T)
#          cat("\n")
#          te <- tt$estimate
#          #tret <- te
#          #mode(tret) <- 'character'
#          #tret[,'p-value'] <- pformat( te[,'p-value'], digits = pround)
#          #if( !is.null(round)) {
#          #    for (i in 1:ncol(te) ) {
#          #        tret[,i] <- rnd(te[,i], digits = round)
#          #    }
#          #print(tret,quote=F)
#          if(!is.null(zhead <- attr(tt,'heading'))) cat(zhead,"\n")
#          print(formatCoefmat(te,digits=round,pdigits=pround),quote=F,right=T)
#          if (L == TRUE ) {
#             cat("\nL:\n")
#             print(tt$L)
#             if( dim(tt$L.full)[1] < dim(tt$L) [1]) {
#              cat("\nL (full rank):\n")
#              print(tt$L.full)
#             }
#          }
#          if ( cov == TRUE ) {
#             cat("\nVar-Cov of estimates:\n")
#             print(tt$vcov)
#             cat("\nCorrelations:\n")
#             print(tt$vcor)
#          }
#
#      }
#      invisible(x)
# }
#
#
# ##z <- glh(fit, list('language'=Lmat(fit, "language")))
# ##z
#
# ##z <- glh(fit, list('Language | Year=1998'=Ldiff(fit, "language",ref="Qc French Multi")))
# ##z
#
#
# ####################################################
# ####################################################  ABOVE is OKAY
#
#
# #' Extended anova
# #'
# #'
# #'
# #'
# #'
# #' @param object
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (object, ...)
# #' UseMethod("xanova")
# #'
# #' @export
# xanova <- function(object,...) UseMethod("xanova")
#
# #xanova.lmer <- function(fit, Llist ) {
# #     if ( is.list(Llist) )
# #}
#
#
#
#
# #' Modified Anova for lmer objects
# #'
# #'
# #'
# #'
# #'
# #' @param fit
# #' @param Llist
# #' @param df
# #' @param clevel
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (fit, Llist, df = NULL, clevel = 0.95)
# #' {
# #'     warning("xanova.lmer uses Chi-Square tests")
# #'     ret <- list()
# #'     for (ii in 1:length(Llist)) {
# #'         L <- rbind(Llist[[ii]])
# #'         QR <- qr(L)
# #'         R <- qr.R(QR)
# #'         dfH <- QR$rank
# #'         eta <- R %*% fixef(fit)
# #'         vv <- R %*% vcov(fit) %*% t(R)
# #'         chisq <- t(eta) %*% qr.solve(vv, eta)
# #'         test <- list(ChiSquare = chisq, DF = dfH, `p-value` = 1 -
# #'             pchisq(chisq, dfH))
# #'         ret[[ii]]$anova <- test
# #'         eta <- L %*% fixef(fit)
# #'         vv <- diag(L %*% vcov(fit) %*% t(L))
# #'         etasd <- sqrt(vv)
# #'         zval <- c(eta/etasd)
# #'         aod <- cbind(Estimate = c(eta), Std.Error = etasd, `z-value` = zval,
# #'             `p-value` = 2 * pnorm(-abs(zval)))
# #'         if (!is.null(clevel)) {
# #'             hw <- qnorm(1 - (1 - clevel)/2) * etasd
# #'             aod <- cbind(aod, LL = eta - hw, UL = eta + hw)
# #'             labs <- paste(c("Lower", "Upper", format(clevel)))
# #'             colnames(aod)[ncol(aod) + c(-1, 0)] <- labs
# #'         }
# #'         aod <- as.data.frame(aod)
# #'         class(aod) <- c("estimate.lme", "data.frame")
# #'         ret[[ii]]$estimate <- aod
# #'     }
# #'   }
# #'
# #' @export
# xanova.lmer <- function( fit, Llist , df = NULL, clevel = .95) {
#        # performs a Wald test on an object that has a fixef and a vcov methods
#        warning("xanova.lmer uses Chi-Square tests")
#        ret <- list()
#        for ( ii in 1:length(Llist) ) {
#            L <- rbind(Llist[[ii]])
#
#            # anova step
#
#            QR <- qr(L)
#            R <- qr.R(QR)
#            dfH <- QR$rank
#            eta <- R %*% fixef(fit)
#            vv <- R %*% vcov(fit) %*% t(R)
#            chisq <- t(eta) %*% qr.solve(vv, eta)
#            test <- list(ChiSquare = chisq, DF = dfH, "p-value" = 1-pchisq(chisq,dfH))
#            ret[[ii]]$anova <- test
#
#            # estimation
#
#            eta <- L %*% fixef(fit)
#            vv <- diag( L %*% vcov(fit) %*% t(L))
#            etasd <- sqrt(vv)
#            zval <- c(eta/etasd)
#            aod <- cbind(Estimate=c(eta), Std.Error = etasd,
#                "z-value" = zval, "p-value" = 2*pnorm(-abs(zval)))
#            if( !is.null(clevel) ) {
#                hw <- qnorm(1-(1-clevel)/2) * etasd
#                aod <- cbind( aod, LL = eta - hw, UL = eta + hw)
#                labs <- paste(c("Lower","Upper",format(clevel)))
#                colnames(aod) [ ncol(aod) + c(-1,0)] <- labs
#            }
#            aod <- as.data.frame(aod)
#            class(aod) <- c('estimate.lme','data.frame')
#
#            ret[[ii]]$estimate <- aod
#         }
# }
#
# ####################################   NEED TO TEST ABOVE
#
#
#
#
#
#
#
#
#
#
#
# #' Read coding tables in RDC
# #'
# #'
# #'
# #'
# #'
# #' @param x
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x)
# #' {
# #'     tonum <- function(x) as.numeric(gsub(",", "", as.character(x)))
# #'     ff <- function(x) format(x, big.mark = ",")
# #'     tran.table <- scan(what = list("character", "integer", "character",
# #'         "character"), flush = T)
# #'     sam <- tonum(tran.table[[3]])
# #'     pop <- tonum(tran.table[[4]])
# #'     samp.pct <- 100 * (sam/pop)/(sum(sam, na.rm = T)/sum(pop,
# #'         na.rm = T))
# #'     print(data.frame(Code = tran.table[[2]], Content = tran.table[[1]],
# #'         Sample = ff(sam), Popn = ff(pop), Sampling.Pct = round(samp.pct,
# #'             1)))
# #'     tran(tran.table[[2]], tran.table[[1]], x, tofactor = T)
# #'   }
# #'
# #' @export
# enc <- function(x) {
# 		## this function will use the coding table in Stats Can documentation
# 		## to create a factor with the right names
# 		#
# 		tonum <- function(x) as.numeric(gsub(",","",as.character(x)))
# 		ff <- function(x) format(x, big.mark=',')
# 		tran.table <- scan(what=list('character','integer','character','character'),
# 			flush = TRUE)
# 		# Brief report
# 		sam <- tonum(tran.table[[3]])
# 		pop <- tonum(tran.table[[4]])
# 		samp.pct <- 100 * (sam / pop) / ( sum(sam,na.rm=T)/sum(pop,na.rm=T))
# 		print( data.frame(
# 			Code = tran.table[[2]], Content = tran.table[[1]], Sample = ff(sam), Popn = ff(pop),
# 				Sampling.Pct = round(samp.pct,1)))
# 		tran( tran.table[[2]], tran.table[[1]], x , tofactor = TRUE)
# 	}
#
#
#
#
#
# #' Transform a frequency table to percentages
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param MARGIN
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, MARGIN = 1)
# #' {
# #'     if (length(dim(x)) == 1) {
# #'         ret <- cbind(N = x, pct = 100 * x/sum(x, na.rm = T))
# #'         ret <- rbind(ret, Total = apply(ret, 2, sum, na.rm = T))
# #'         print(round(ret, 1))
# #'         return(invisible(ret))
# #'     }
# #'     ret <- list(N = atotal(x), pct = 100 * acond(x, MARGIN))
# #'     cat("\nN:\n")
# #'     print(ret[[1]])
# #'     cat("\nPercentage:\n")
# #'     print(round(ret[[2]], 1))
# #'     invisible(ret)
# #'   }
# #'
# #' @export
# apct <- function(x,MARGIN=1) {
# 		if( length(dim(x)) == 1) {
# 			ret <- cbind( N = x, pct = 100*x/sum(x,na.rm=T))
# 			ret <- rbind( ret, Total = apply(ret,2,sum,na.rm=T))
# 			print( round(ret,1))
# 			return(invisible(ret))
# 		}
# 		# report a table
# 		ret <- list( N=atotal(x), pct = 100 * acond(x,MARGIN) )
# 		cat("\nN:\n")
# 		print( ret[[1]])
# 		cat("\nPercentage:\n")
# 		print( round( ret[[2]],1))
# 		invisible(ret)
# 	}
#
#
#
#
#
# #' Percentages of a column or row sum
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param MARGIN
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, MARGIN = 1)
# #' {
# #'     if (length(dim(x)) == 1) {
# #'         ret <- cbind(N = x, pct = 100 * x/sum(x, na.rm = T))
# #'         ret <- rbind(ret, Total = apply(ret, 2, sum, na.rm = T))
# #'         print(round(ret, 1))
# #'         return(invisible(ret))
# #'     }
# #'     ret <- list(N = atotal(x), pct = 100 * acond(x, MARGIN))
# #'     cat("\nN:\n")
# #'     print(ret[[1]])
# #'     cat("\nPercentage:\n")
# #'     print(round(ret[[2]], 1))
# #'     invisible(ret)
# #'   }
# #'
#   arep <- apct   # old names
#
# ## From library gm
#
#
#
# #' Modified version of write.sas
# #'
# #'
# #'
# #'
# #'
# #' @param df
# #' @param datafile
# #' @param codefile
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (df, datafile = "G:/SAS/R2SAS.txt", codefile = "G:/SAS/R2SAS.sas")
# #' {
# #'     debug <- F
# #'     pr <- function(x) {
# #'         cat(deparse(substitute(x)), "\n")
# #'         print(x)
# #'         cat("==========\n")
# #'         cat(x)
# #'         cat("\n=============\n")
# #'     }
# #'     lrecl <- 256
# #'     if (!debug) {
# #'         write.table(df, file = datafile, row = FALSE, col = FALSE,
# #'             sep = ";", na = ".")
# #'         lines <- scan(file = datafile, what = list("character"),
# #'             sep = "\n")
# #'         lrecl <- max(nchar(lines))
# #'     }
# #'     nms <- names(df)
# #'     nms.sas <- gsub("\.", "_", nms)
# #'     if (length(unique(toupper(nms.sas))) != length(nms.sas)) {
# #'         ind <- duplicated(toupper(nms.sas))
# #'         ind.rev <- duplicated(rev(toupper(nms.sas)))
# #'         cat("Warning:\n")
# #'         cat("The following R names may yield duplicate SAS names",
# #'             "\n", paste(nms[ind | rev(ind.rev)], collapse = " "),
# #'             "\n")
# #'         warning("Possible duplicate SAS names")
# #'     }
# #'     factors <- sapply(df, is.factor) | sapply(df, is.character)
# #'     classes <- sapply(df, class)
# #'     odd.classes <- setdiff(sapply(df, class), c("numeric", "factor",
# #'         "character"))
# #'     if (length(odd.classes) > 0) {
# #'         cat("Warning:\n")
# #'         cat("The following variables have classes that might not be handled properly by SAS\n")
# #'         print(classes[grep(odd.classes, classes)])
# #'         cat("\n")
# #'     }
# #'     factor.names <- nms[factors]
# #'     factor.names.sas <- nms.sas[factors]
# #'     dollarsign <- ifelse(factors, "$", "")
# #'     factor.lengths <- sapply(df[factor.names], function(x) {
# #'         if (is.factor(x))
# #'             max(nchar(levels(x)))
# #'         else max(nchar(x))
# #'     })
# #'     length.stmt <- paste(paste("   ", factor.names.sas, "$",
# #'         factor.lengths, "\n"), collapse = "")
# #'     length.stmt <- paste("LENGTH\n", length.stmt, ";\n")
# #'     if (debug)
# #'         pr(length.stmt)
# #'     input.stmt <- paste(paste("    ", nms.sas, dollarsign, "\n"),
# #'         collapse = "")
# #'     input.stmt <- paste("INPUT\n", input.stmt, ";\n")
# #'     if (debug)
# #'         pr(input.stmt)
# #'     code <- paste("filename r2sas '", datafile, "';\n", "libname to 'G:/SAS';\n",
# #'         "data to.r2sas;\n", "infile r2sas delimiter=';' dsd LRECL =",
# #'         lrecl + 100, ";\n", sep = "")
# #'     code <- paste(code, length.stmt, input.stmt, "\nrun;\n")
# #'     if (debug)
# #'         pr(code)
# #'     if (!debug)
# #'         cat(code, file = codefile)
# #'     invisible(0)
# #'   }
# #'
# #' @export
# write.sas <- function( df , datafile="G:/SAS/R2SAS.txt", codefile = "G:/SAS/R2SAS.sas"){
# 	debug <- F
# 	pr <- function(x) {
# 		cat(deparse(substitute(x)),"\n")
# 		print(x)
# 		cat("==========\n")
# 		cat(x)
# 		cat("\n=============\n")
# 	}
#
# 	lrecl <- 256
# 	if(!debug) {
# 		write.table(df, file = datafile, row = FALSE, col = FALSE,
#         		sep = ";", na = ".")
# 		# compute lrecl
# 		lines <- scan(file = datafile, what = list("character"), sep = "\n")
# 		lrecl <- max(nchar(lines))
# 	}
#     	nms <- names(df)
# 	nms.sas <- gsub("\\.","_",nms)
# 	## Check for duplicate SAS names
#
# 	if ( length( unique(toupper(nms.sas))) != length(nms.sas)) {
# 		ind <- duplicated( toupper(nms.sas))
# 		ind.rev <- duplicated( rev(toupper(nms.sas)))
# 		cat("Warning:\n")
# 		cat("The following R names may yield duplicate SAS names",
# 			"\n", paste(nms[ind | rev(ind.rev)],collapse=" "),"\n")
# 		warning("Possible duplicate SAS names")
# 	}
#
# 	factors <- sapply(df, is.factor) | sapply(df, is.character)
# 	## check for odd types
# 	classes <- sapply(df, class)
# 	odd.classes <- setdiff( sapply(df,class), c('numeric','factor','character') )
# 	if ( length(odd.classes) > 0 ) {
# 		cat("Warning:\n")
# 		cat("The following variables have classes that might not be handled properly by SAS\n")
# 		print( classes[ grep( odd.classes, classes)])
# 		cat("\n")
# 	}
#
# 	factor.names <- nms[factors]
# 	factor.names.sas <- nms.sas[factors]
# 	dollarsign <- ifelse( factors, "$","")
# 	factor.lengths <- sapply( df[factor.names], function(x) {
# 		if(is.factor(x)) max(nchar(levels(x) ) ) else
# 			max(nchar(x))
# 		})
# 	length.stmt <- paste(paste( "   ", factor.names.sas, "$", factor.lengths,"\n"),collapse = "")
# 	length.stmt <- paste( "LENGTH\n", length.stmt, ";\n")
# 	if (debug) pr(length.stmt)
# 	input.stmt <- paste(paste("    ", nms.sas, dollarsign,"\n"), collapse = "")
# 	input.stmt <- paste( "INPUT\n", input.stmt, ";\n")
# 	if (debug) pr(input.stmt)
#
# 	code <- paste("filename r2sas \'",datafile,"\';\n",   # might have to convert to backslashes
# 		"libname to \'G:/SAS\';\n",
# 		"data to.r2sas;\n",
# 		"infile r2sas delimiter=\';\' dsd LRECL =", lrecl+100, ";\n", sep = "")
# 	code <- paste(code,length.stmt, input.stmt, "\nrun;\n")
# 	if (debug) pr(code)
# 	if(!debug) cat(code, file = codefile)
# 	invisible(0)
# }
#
#
# ## date()
# ## write.sas(dd[dd$wave > 1,setdiff(sort(names(dd)),c('qday','bday')) ]) # approx 1 min.
# ## date()
#
# if(F) { # current version in /R/coursefun.R
# td <- function( basecol = NULL, col = c(3,5,4,6,7,8,2), lty = 1:7,
# 	lwd = 1, pch = 1:7, cex = 0.8, font = 1, len = 7, long = FALSE,
# 	new = FALSE, record = TRUE,  theme = col.whitebg,
# 	col.symbol = col, col.line = col,
# 	colsets = c( "plot.symbol","plot.line", "dot.symbol","dot.line","cloud.3d","box.dot"),...) {
# 	require(lattice)
# 	if (long) {
# 		col <- c(3,5,4,6,8,2)
# 		len <- 42
# 	}
# 	if (new) trellis.device( theme = theme, record = record, new = new, ...)
# 	len <- max( len, length(col.symbol), length(col.line), length( lty), length(lwd), length(pch),
# 		length(cex), length(font))
# 	spl <- trellis.par.get("superpose.line")
# 	spl$lty <- if ( is.null(lty)) spl$lty else lty
# 	spl$col <- if ( is.null(col.line)) spl$col else col.line
# 	spl$lwd <- if ( is.null(lwd)) spl$lwd else lwd
# 	trellis.par.set("superpose.line", Rows(spl, 1:len))
# 	sps <- trellis.par.get("superpose.symbol")
# 	sps$pch <- if ( is.null(pch)) sps$pch else pch
# 	sps$col <- if ( is.null(col.symbol)) sps$col else col.symbol
# 	sps$cex <- if ( is.null(cex)) sps$cex else cex
# 	sps$font <- if ( is.null(font)) sps$font else font
# 	trellis.par.set("superpose.symbol", Rows(sps, 1:len))
# 	if( !is.null(basecol) ) {
# 		for(ii in colsets) {
# 			tt <- trellis.par.get(ii)
# 			tt$col <- basecol
# 			trellis.par.set(ii,tt)
# 		}
# 	}
# 	invisible( attr( .Device, "trellis.settings"))
# }
# } # end of F
#
#
#
#
#
#
# #' Capitalize first character of each word
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param tofactor
# #' @param stop
# #' @param blanks
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, tofactor = is.factor(x), stop = c(" The", " Of",
# #'     " By", " To", " And"), blanks = c(" ", "(", "\"", "/", "+"))
# #' {
# #'     under2blank <- T
# #'     if (is.factor(x)) {
# #'         ret <- cap1(levels(x))
# #'         if (length(unique(ret)) != length(ret))
# #'             warning("factor levels have been shortened")
# #'         levels(x) <- ret
# #'         return(x)
# #'     }
# #'     ret <- as.character(x)
# #'     for (ii in 1:length(ret)) {
# #'         z <- ret[[ii]]
# #'         if (under2blank)
# #'             z <- gsub("_", " ", z)
# #'         n <- nchar(z)
# #'         z <- substring(z, 1:n, 1:n)
# #'         zu <- toupper(z)
# #'         zl <- tolower(z)
# #'         zb <- c(" ", zu[-n])
# #'         z <- paste(ifelse(zb %in% blanks, zu, zl), collapse = "")
# #'         for (ss in stop) z <- gsub(ss, tolower(ss), z)
# #'         ret[[ii]] <- z
# #'     }
# #'     ret
# #'   }
# #'
# #' @export
# cap1 <- function(x, tofactor = is.factor(x),
# 		stop=c(" The"," Of"," By"," To"," And"),
# 		blanks=c(" ","(","\"","/","+")) {
# 	# capitalizes first letters
# 	under2blank <- T
# 	if ( is.factor(x)) {
# 		ret <- cap1(levels(x))
# 		if ( length(unique(ret)) != length(ret)) warning("factor levels have been shortened")
# 		levels(x) <- ret
# 		return(x)
# 	}
# 	ret <- as.character(x)
# 	for ( ii in 1:length(ret)) {
# 		z <- ret[[ii]]
# 		if(under2blank) z <- gsub("_"," ",z)
# 		n <- nchar(z)
# 		z <- substring( z, 1:n, 1:n)
# 		zu <- toupper(z)
# 		zl <- tolower(z)
# 		zb <- c(" ",zu[-n])
# 		z <- 	paste( ifelse(zb %in% blanks, zu, zl), collapse ="")
# 		for ( ss in stop) z <- gsub(ss,tolower(ss), z)
# 		ret[[ii]] <- z
# 	}
# 	ret
# }
#
#
#

# #' @export
# abind.rdc <- function( arr1, arr2, d, facename = "") {
# 	# glue arr1 to arr2 along dimension d (i.e. face 'not d')
# 	# copied from library gm 05 05 03
# 	d1 <- dim(arr1)
# 	n1 <- length(d1)
# 	d2 <- dim(arr2)
# 	n2 <- length(d2)
# 	dn1 <- dimnames( arr1 )
# 	dn2 <- dimnames( arr2 )
#
# 	arenull <- is.null(dn1) & is.null(dn2)
# 	if ( is.null(dn1)) {
# 		dn1 <- lapply( as.list(d1), function(x) seq(1,x))
# 		dimnames(arr1) <- dn1
# 	}
# 	if ( n1 != n2 ) {
# 		d2 <- d1
# 		d2[d] <- 1
# 		dn2 <- dn1
# 		dn2[[d]] <- facename
# 		dimnames(arr2) <- NULL
# 		dim(arr2) <- d2
# 		dimnames(arr2) <- dn2
# 		n2 <- n1
# 	}
# 	if ( is.null(dn2)) {
# 		dn2 <- lapply( as.list(d2), function(x) seq(1,x))
# 		dimnames(arr2) <- dn2
# 	}
# 	perm <- 1:n1
# 	perm[c(d,n1)] <- c(n1,d)  # perm is an involution
#
# 	arr.p1 <- aperm(  arr1, perm )
#
# 	arr.p2 <- aperm(  arr2, perm )
# 	dret <- d1[perm]
# 	dret[n1] <- dret[n1] + d2[d]
# 	dnret <- dn1
# 	dnret[[d]] <- c(dnret[[d]],dn2[[d]])
# 	ret <- c(arr.p1, arr.p2)
# 	dim(ret) <- dret
#
# 	ret <- aperm( ret, perm )
# 	dimnames( ret ) <- dnret
# 	ret
# }
#
#
#
# # acond( with(dds, table( lang, prov)),2)
#
#
#
#
#
# #' Version of atotal in RDC
# #'
# #'
# #'
# #'
# #'
# #' @param arr
# #' @param FUN
# #' @param name
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (arr, FUN = sum, name = "Total", ...)
# #' {
# #'     d <- dim(arr)
# #'     if (length(d) == 1) {
# #'         arr <- c(arr)
# #'         d <- dim(arr)
# #'     }
# #'     if (is.character(FUN))
# #'         FUN <- get(FUN, mode = "function")
# #'     else if (mode(FUN) != "function") {
# #'         farg <- substitute(FUN)
# #'         if (mode(farg) == "name")
# #'             FUN <- get(farg, mode = "function")
# #'         else stop(paste("\"", farg, "\" is not a function", sep = ""))
# #'     }
# #'     if (is.null(d)) {
# #'         ret <- c(arr, FUN(arr))
# #'         names(ret)[length(ret)] = name
# #'         return(ret)
# #'     }
# #'     n <- length(d)
# #'     name <- rep(name, length = n)
# #'     ret <- arr
# #'     ind <- 1:n
# #'     for (i in n:1) {
# #'         new <- apply(ret, ind[-i], FUN, ...)
# #'         ret <- abind(ret, new, i, name[i])
# #'     }
# #'     ret
# #'   }
# #'
# #' @export
# atotal.rdc <- function( arr, FUN = sum, name = "Total",...) {
# 	# copied from library gm  05 05 03
# 	# 05 05 03: added option for name
# 	d <- dim(arr)
# 	if ( length(d) == 1) {
# 		arr <- c(arr)
# 		d <- dim(arr)
# 	}
# 	if ( is.character( FUN ) ) FUN <- get(FUN,mode='function')
# 	else if( mode(FUN) != "function") {
# 		farg <- substitute(FUN)
# 		if (mode(farg) == 'name' ) FUN <- get(farg,mode='function')
# 		else stop(paste("\"", farg, "\" is not a function", sep=""))
# 	}
# 	if ( is.null(d) ) {
# 		ret <- c(arr, FUN(arr))
# 		names(ret)[length(ret)] = name
# 		return(ret)
# 	}
# 	n <- length (d)
# 	name <- rep(name, length= n)
# 	ret <- arr
# 	ind <- 1:n
# 	for ( i in n:1 ) {
# 		new <- apply(ret, ind[ -i ], FUN,...)
# 		ret <- abind( ret, new, i, name[i])
# 	}
# 	ret
# }
#
# # aprop and acond moved to tab.R
#
# ###
# ###  Utility functions for regression
# ###
#
#     ##
#     ##  na.fitted na.resid   pads with NA to fit original data frame
#     ##
#     ##  extend to new classes by writing a 'na.pad' methoc
#
#
#
#
# #' Fitted values with NAs
# #'
# #'
# #'
# #'
# #'
# #' @param fit
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (fit, ...)
# #' na.pad(fit, fitted(fit, ...))
# #'
# #' @export
#     na.fitted <- function(fit,...) na.pad( fit, fitted(fit,...))
#
#
# #' Add NAs -- resid
# #'
# #'
# #'
# #'
# #'
# #' @param fit
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (fit, ...)
# #' na.pad(fit, resid(fit, ...))
# #'
# #' @export
#     na.resid <- function(fit,...) na.pad( fit, resid(fit,...))
#
#
# #' Residuals with NAs
# #'
# #'
# #'
# #'
# #'
# #' @param fit
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (fit, ...)
# #' na.pad(fit, residuals(fit, ...))
# #'
# #' @export
#     na.residuals <- function(fit,...) na.pad( fit, residuals(fit,...))
#     # na.predict <- function(fit,...) na.pad( fit, predict(x,...))
#
#
#
# #' Add NAs to match original data frame
# #'
# #'
# #'
# #'
# #'
# #' @param fit
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (fit, ...)
# #' UseMethod("na.pad")
# #'
# #' @export
#     na.pad <- function(fit,...) UseMethod('na.pad')
#
#
# #' Add NAs -- lme
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param obj
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, obj)
# #' {
# #'     ind <- rep(NA, nrow(x$data))
# #'     names(ind) <- rn <- rownames(x$data)
# #'     ind.part <- 1:nrow(x$resid)
# #'     names(ind.part) <- rownames(x$resid)
# #'     ind[names(ind.part)] <- ind.part
# #'     if (!is.null(dim(obj)))
# #'         ret <- obj[ind, ]
# #'     else ret <- obj[ind]
# #'     names(ret) <- rn
# #'     ret
# #'   }
# #'
# #' @export
#     na.pad.lme <- function( x, obj) {
#         ind <- rep(NA, nrow(x$data))
#         names(ind) <- rn <- rownames(x$data)
#
#         ind.part <- 1:nrow(x$resid)
#         names(ind.part) <- rownames(x$resid)
#         ind[ names(ind.part) ] <- ind.part
#         # disp(ind)
#         if( !is.null( dim (obj)) ) ret <- obj[ ind,]
#         else ret <- obj[ind]
#         names(ret) <- rn
#         ret
#     }
#
#
#
# #' Add NAs -- default
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param obj
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, obj)
# #' {
# #'     stop(paste("You need to write a method for na.pad for objects of class",
# #'         class(x), "\n", "see na.pad.lme for an example"))
# #'   }
# #'
# #' @export
#     na.pad.default <- function(x,obj) {
#         stop(paste("You need to write a method for na.pad for objects of class",class(x),"\n",
#             "see na.pad.lme for an example"))
#     }
#
#
#
#
#
# # Lagging
# #    - interesting if 'index' is continuous and/or if there
# #      are missing waves.
# #    - need some kind of intrapolation to compute a derivative
# #      at the time of observation
#
# # The first function here only works well with complete data
# # and an integer index ranging from 1 to n. (n can vary from
# # subject to subject)
#
#
#
#
# #' Lag within subject: older less efficient version of Lag
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param id
# #' @param idx
# #' @param lag
# #' @param at
# #' @param check
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, id, idx, lag = 1, at = NULL, check = T)
# #' {
# #'     if (check) {
# #'         comb <- paste(id, idx)
# #'         if (any(duplicated(comb)))
# #'             stop("id not unique in each level of idx")
# #'     }
# #'     if (any(is.na(idx)))
# #'         stop("NAs in idx")
# #'     if (any(round(idx) != idx))
# #'         stop("Non integers in idx")
# #'     ret <- x
# #'     id <- as.character(id)
# #'     names(x) <- id
# #'     for (i in max(idx):min(idx)) {
# #'         to.pos <- idx == i
# #'         if (is.null(at))
# #'             from <- x[idx == (i - lag)]
# #'         else from <- x[idx == at]
# #'         ids <- names(x[to.pos])
# #'         ret[to.pos] <- from[ids]
# #'     }
# #'     ret
# #'   }
# #'
# #' @export
# Lag.0 <- function(x,id,idx,lag = 1,at= NULL, check=T) {
# 	# computes Lagged values but without intrapolation
# 	# adds 'at'
#
# 	if (check) {
# 		comb <- paste(id,idx)
# 		if (any(duplicated(comb))) stop("id not unique in each level of idx")
# 	}
# 	if(any(is.na(idx))) stop("NAs in idx")
# 	if(any( round(idx) != idx)) stop("Non integers in idx")
# 	ret <- x
# 	id <- as.character(id)
# 	names(x) <- id
# 	for ( i in max(idx):min(idx)) {
# 		to.pos <- idx == i
# 		if( is.null(at) ) from <- x[idx == (i-lag)]
# 		else from <- x[idx == at ]
# 		ids <- names( x[to.pos])
# 		ret[to.pos] <- from[ids]
# 	}
# 	ret
# }
#
# # Previous version
# # cLag <- function(x, id = rep(1,length(x)), time = 1:length(x), lag=1, at = NULL, check = T, idx) {
# #    # renamed cLag to avoid conflict with Hmisc::Lag
# #    ## Takes approx. 30% of time of version Lag.0 on a small problem
# #   if (check) {
# # 		comb <- paste(id,time)
# # 		if (any(duplicated(comb))) stop("id not unique in each level of time")
# # 	}
# # 	if( !missing(idx)) time = idx    # for compatibility with previous argument names
# #   ret <- x
# #   names(ret) <- paste(id,time,sep='::')
# #   retnn <- if(is.null(at)) paste(id,time - lag,sep='::')  else paste(id,at,sep="::")
# #   ret [ retnn ]
# # }
# #
# #' @export
# cLag <-
#   function (x, id = rep(1, length(x)), time = capply(id,id,function(x) 1:length(x)), lag = 1,
#             at = NULL, check = T, idx)
#   {
#     if (check) {
#       comb <- paste(id, time)
#       if (any(duplicated(comb)))
#         stop("id not unique in each level of time")
#     }
#     if (!missing(idx))
#       time = idx
#     ret <- x
#     names(ret) <- paste(id, time, sep = "::")
#     retnn <- if (is.null(at))
#       paste(id, time - lag, sep = "::")
#     else paste(id, at, sep = "::")
#     ret[retnn]
#   }
#
#
#
# #' 'Contextual' Lag with respect to time within id
# #'
# #' Lag a vector with respect to time order
# #'
# #' The function can also be called as \code{Lag} which is now deprecated to
# #' avoid conflicts with \code{Hmisc::Lag}. Use \code{cLagI} for
# #' intra(extra)polative lagging.
# #'
# #' @aliases cLag Lag
# #' @param x vector of values to be lagged
# #' @param id identifies clusters within which lagging takes place
# #' @param time time variable for lagging, generally integer valued
# #' @param lag period for lagging: lag = 1 reaches one unit of \code{time} into
# #' the past. Negative numbers reach into the future.
# #' @param at uses the value of \code{x} at a particular value of \code{time}.
# #' e.g. \code{at = 1} returns the value of \code{x} corresponding to \code{time
# #' == 1}.
# #' @param check that \code{id/time} combinations are unique.
# #' @param idx for compatibility with a previous version.
# #' @return A vector of values of \code{x} lagged within levels of \code{id}.
# #' The value returned at \code{time == t} is the value of \code{x} at
# #' \code{time == t-1}.  If there is no observation for which \code{time == t-1}
# #' then the returned value is \code{NA}.  To lag to the previous value of
# #' \code{time}, one can use \code{rank}.  Consider, also, \code{\link{cLagI}}
# #' that intrapolates backwards one unit of \code{time}.
#
# #' @author Georges Monette
# #' @seealso \code{\link{cLagI}}, \code{\link{cDiffI}}, \code{\link{capply}},
# #' \code{\link{up}}
#
# #' @keywords ~kwd1 ~kwd2
# #' @export
# Lag <- cLag # historical name that conflicts with Hmisc
#
#
# # ## small test of Lag
# # zx <- c(2,3,2,4,2,4, 5,3,4,5,7,8,9)
# # zid<- c(1,1,2,2,3,3 ,3,4,4,4,4,4,4)
# # cbind( zid, zx, Lag(zx,zid,zx))
#
#
#
# #' @export
# cLagI <- function(x,id,time,lag=1,delta=.01,check=T) {
#     # renamed from LagI to avoid conflict with Hmisc
# 	# lags by intrapolating
# 	# with complete data at each value of index, this does the same thing
# 	# as Lag
# 	# If values of Lag are skipped then we get linear intrapolations.
# 	# Note that 'delta' should be small enough so that values of x are
# 	# at least delta apart. However, too small a value for delta introduces
# 	# numerical error
# 	#
# 	if (check) {
# 		comb <- paste(id,time)
# 		if (any(duplicated(comb))) stop("id not unique in each level of idx")
# 	}
# 	ret <- x
# 	id <- as.character(id)
# 	names(x) <- id
# 	names(time) <- id
# 	for (nn in unique(id)){
# 		pos <- id == nn
# 		xx <- x[pos]
# 		tt <- time[pos]
# 		topred <- tt-delta
# 		drop <- is.na(xx)|is.na(tt)
# 		xxc <- xx[!drop]
# 		ttc <- tt[!drop]
# 		nl <- length(xxc)
# 		if ( nl > 0) {
# 			if ( nl > 1 ) xx.p <- approx(ttc,xxc,topred)$y
# 			else xx.p <- NA
# 			xx.lag <- xx - lag*(xx - xx.p)/delta
# 			ret[pos] <- xx.lag
# 		}
# 	}
# 	ret
# }
#
#
#
# #' Lag with respect to time within id with interpolation for non-integer time
# #'
# #'
# #'
# #'
# #'
# #' @aliases LagI cLagI DiffI
# #' @param x the values to be lagged.
# #' @param id values are lagged within each level of id
# #' @param time measure of time used for lagging
# #' @param lag distance to lag: e.g. -1 takes the value of \code{x} 1 unit of
# #' \code{time} in the past.
# #' @param delta increment used for extrapolation
# #' @param check uniqueness of \code{id}/\code{time} combinations
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, id, time, lag = 1, delta = 0.01, check = T)
# #' {
# #'     if (check) {
# #'         comb <- paste(id, time)
# #'         if (any(duplicated(comb)))
# #'             stop("id not unique in each level of idx")
# #'     }
# #'     ret <- x
# #'     id <- as.character(id)
# #'     names(x) <- id
# #'     names(time) <- id
# #'     for (nn in unique(id)) {
# #'         pos <- id == nn
# #'         xx <- x[pos]
# #'         tt <- time[pos]
# #'         topred <- tt - delta
# #'         drop <- is.na(xx) | is.na(tt)
# #'         xxc <- xx[!drop]
# #'         ttc <- tt[!drop]
# #'         nl <- length(xxc)
# #'         if (nl > 0) {
# #'             if (nl > 1)
# #'                 xx.p <- approx(ttc, xxc, topred)$y
# #'             else xx.p <- NA
# #'             xx.lag <- xx - lag * (xx - xx.p)/delta
# #'             ret[pos] <- xx.lag
# #'         }
# #'     }
# #'     ret
# #'   }
# #'
# #' @export
# LagI <- cLagI
#
# #' @export
# cDiffI <- function(xx,...) {
# 	xx - LagI(xx,...)
# }
#
#
#
# #' Difference between current value and lagged value within subject.
# #'
# #'
# #'
# #'
# #'
# #' @param xx
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (xx, ...)
# #' {
# #'     xx - LagI(xx, ...)
# #'   }
# #'
# #' @export
# DiffI <- cDiffI
#
# # tt <- c(1,2,3,1,2,3,1,2,3)
# # xx <- c(1,2,3,1,2,3,1,2,3)
# # id <- rep(letters[1:3],rep(3,3))
# # xx-LagI(xx,id,tt)
# # DiffI(xx,id,tt,delta=.1)
#
# ###
# ###  Splines and polynomial: quadratic and cubic
# ###
#
#
# #' @export
# qs <- function(x, knots=quantile(x,pc), exclude = 0, pc = c(.25,.75)) {
# 	# simple quadratic spline
# 	ret <- cbind(x,x^2)
# 	nam <- deparse(substitute(x))
#     if ( missing(knots) ) warning("May be unsafe for prediction with newdata")
# 	for ( kk in knots ) {
# 		z <- x - kk
# 		z [z<0] <- 0
# 		ret <- cbind(ret, z^2)
# 	}
#
# 	if ( is.null(knots) ) dimnames(ret)[[2]] <- paste(nam,c('','^2'), sep='')
# 	else dimnames(ret)[[2]] <- paste(nam,c('','^2',paste('.',knots,sep='')),sep='')
# 	if (exclude > 0) ret <- ret[, -( 1:exclude)]
# 	ret
#
# 	ret
# }
#
#
#
# #' @export
# lsp <- function(x, knots=quantile(x,pc), exclude = 0, pc = c(.25,.75)) {
# 	# linear spline
# 	ret <- cbind(x)
# 	nam <- deparse(substitute(x))
#     if ( missing(knots) ) warning("May be unsafe for prediction with newdata")
# 	for ( kk in knots ) {
# 		z <- x - kk
# 		z [z<0] <- 0
# 		ret <- cbind(ret, z)
# 	}
# 	if ( is.null(knots) ) dimnames(ret)[[2]] <- paste(nam,c(''), sep='')
# 	else dimnames(ret)[[2]] <- paste(nam,c('',paste('.',knots,sep='')),sep='')
# 	if (exclude > 0) ret <- ret[, -( 1:exclude)]
# 	ret
#
# 	ret
# }
#
# #' @export
# cs <- function(x, knots=quantile(x,pc), exclude = 0, pc = c(.25,.75))  {
# 	# simple cubic spline
#     if ( missing(knots) ) warning("May be unsafe for prediction with newdata")
# 	ret <- cbind(x,(x)^2,(x)^3)
# 	nam <- deparse(substitute(x))
# 	for ( kk in knots ) {
# 		z <- (x) - kk
# 		z [z<0] <- 0
# 		ret <- cbind(ret, z^3)
# 	}
# 	if ( is.null(knots) ) dimnames(ret)[[2]] <- paste(nam,c('','^2','^3'), sep='')
# 	else dimnames(ret)[[2]] <- paste(nam,c('','^2','^3',paste('.',knots,sep='')),sep='')
# 	if (exclude > 0) ret <- ret[, -( 1:exclude)]
# 	ret
# }
#
#
#
# #' Matrix of powers
# #'
# #' @param x
# #' @param order
# #' @param exclude
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, order = 2, exclude = 0)
# #' {
# #'     ret <- cbind(rep(1, length(x)))
# #'     for (i in 1:order) ret <- cbind(ret, x^i)
# #'     nam <- deparse(substitute(x))
# #'     powers <- paste(nam, "^", 0:order, sep = "")
# #'     powers[1] <- "Intercept"
# #'     powers[2] <- nam
# #'     dimnames(ret)[[2]] <- powers
# #'     if (!is.null(exclude))
# #'         ret <- ret[, -(1:(exclude + 1))]
# #'     ret
# #'   }
# #'
# #' @export
# Poly <- function(x, order=2, exclude = 0)  {
# 	# polynomial
# 	# exclude: NULL to include intercept, otherwise exclude order up to
# 	# including exclude
# 	ret <- cbind(rep(1,length(x)))
# 	for ( i in 1:order) ret <- cbind(ret, x^i)
# 	nam <- deparse(substitute(x))
# 	powers <- paste(nam,"^",0:order,sep="")
# 	powers[1] <- "Intercept"
# 	powers[2] <- nam
# 	dimnames(ret)[[2]] <- powers
# 	if (!is.null(exclude)) ret <- ret[, -( 1:(exclude+1))]
# 	ret
# }
#
#
#
# if (F) {
# 	## shows that columns of cs span same space as columns of bs
#
# 	zz <- 0:20
# 	cs(zz, c(5,15))
#
# 	#search()
# 	#  detach(9)
# 	### comparison with ns
#
# 	zd <- data.frame( x = 0:20, y = (-10:10)^3 + rnorm(21))
#
# 	bss <- bs(zd$x,knots = c(5,15))
# 	css <- cs(zd$x,knots = c(5,15))
#
# 	fit <- lm(bss ~ css - 1)
# 	round(coef(fit),7)
# 	summary(fit)
#
# 	fit <- lm(css ~ bss - 1)
# 	summary(fit)
# 	round(coef(fit),7)
#
# 	fit <- lm(y ~ 1+cs(x,c(5,15)), zd)
# 	summary(fit)
# 	anova(fit)
#
# 	fit <- lm(y ~ x + I(x^2) + cs(x,c(5,15),2), zd)
# 	summary(fit)
# 	anova(fit)
#
# 	fit <- lm(y ~ x + I(x^2) + I(x^3) + cs(x,c(5,15),3), zd)
# 	summary(fit)
# 	anova(fit)
#
# 	fit <- lm(y ~ x + I(x^2) + cs(x,pc=c(.05,.95)), zd, singular.ok=T)
# 	summary(fit)
# 	anova(fit)
#
# 	qqr <- qr.Q(qr(cs(0:20,c(5,15))))
# 	fit <- lm( css ~ qqr-1)
# 	summary(fit)
# 	round(coef(fit),5) == round(qr.R(qr(cs(0:20, c(5,15)))),5)
#
# }
#
#
# ## much better approach to stack:
#
#
#
# #' Stack data frames
# #'
# #'
# #'
# #'
# #'
# #' @param \dots
# #' @param vname
# #' @param oname
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (..., vname = ".type", oname = ".order")
# #' {
# #'     z <- list(...)
# #'     for (ii in 1:length(z)) {
# #'         z[[ii]][, vname] <- rep(ii, nr <- nrow(z[[ii]]))
# #'         z[[ii]][, oname] <- 1:nr
# #'     }
# #'     ret <- z[[1]]
# #'     for (ii in 2:length(z)) ret <- merge(ret, z[[ii]], all = T)
# #'     ret[order(ret[[vname]], ret[[oname]]), ]
# #'   }
# #'
# #' @export
# mergec <- function( ... ,vname = '.type',oname = ".order") {
#
# 	## stacks data frames adding a new variable .type to identify each source data frame
# 	## good way to combine data and prediction data frames to show data and fitted
# 	## values in a common plot
#
# 	z <- list(...)
# 	for ( ii in 1:length(z)) {
# 		z[[ii]][,vname] <- rep(ii, nr <- nrow( z[[ii]]))
#         z[[ii]][,oname] <- 1:nr
# 	}
# 	ret <- z[[1]]
# 	for ( ii in 2:length(z)) ret <- merge(ret,z[[ii]], all = TRUE)
# 	ret [ order( ret[[vname]], ret[[oname]]),]
# }
#
# if(F) {
# 	## test mergec
# 	zd1 <- data.frame(x=1:4, y = letters[1:4], z = 1:4,d=LETTERS[1:4])
# 	zd2 <- data.frame(x=5:7, y = letters[5:7])
# 	zd3 <- data.frame(x=8:9,  w=1:2, d=1:2)
# 	sapply(zz <- mergec(zd1,zd2,zd3),class)
# 	zz
# 	levels(zz$d)
# }
#
# ###
# ###  extended merge for typical merging of data frames by 'ID' with possible
# ###  other common variables
# ###
# ###
#
#
#
# #' Extended merge with diagnostics
# #'
# #' Extended merge with diagnostics. This is a modification of \code{merge} that
# #' combines consistent variables even if not specified in 'by' to keep a common
# #' name.
# #'
# #'
# #'
# #' @param x,y data frames, or objects to be coerced to one
# #' @param by
# #' @param all
# #' @param dropdots
# #' @param verbose
# #' @param debug
# #' @param from
# #' @param \dots
# #' @author Georges Monette
# #' @seealso \code{\link[base]{merge}}
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, y, by, all = T, dropdots = F, verbose = F, debug = T,
# #'     from = F, ...)
# #' {
# #'     help <- "This is a modification of merge that combines consistent variables\neven if not specified in 'by' to keep a common name.\n-- Some errors fixed Apr 24, 2007"
# #'     xm <- function(a, b, tofac = is.factor(a) || is.factor(b)) {
# #'         if (tofac) {
# #'             levs <- union(levels(b), levels(a))
# #'             a <- as.character(a)
# #'             b <- as.character(b)
# #'         }
# #'         b[is.na(b)] <- a[is.na(b)]
# #'         if (tofac) {
# #'             levs <- union(levs, unique(b))
# #'             b <- factor(b, levels = levs)
# #'         }
# #'         b
# #'     }
# #'     na2f <- function(x) {
# #'         x[is.na(x)] <- F
# #'         x
# #'     }
# #'     consistent <- function(a, b) {
# #'         if (is.factor(a))
# #'             a <- as.character(a)
# #'         if (is.factor(b))
# #'             b <- as.character(b)
# #'         !na2f(a != b)
# #'     }
# #'     if (from) {
# #'         xname <- deparse(substitute(x))
# #'         yname <- deparse(substitute(y))
# #'         x[[".F"]] <- rep("x", nrow(x))
# #'         y[[".F"]] <- rep("y", nrow(y))
# #'     }
# #'     xby <- x[, by, drop = F]
# #'     yby <- y[, by, drop = F]
# #'     xby$.file <- rep("x", nrow(xby))
# #'     yby$.file <- rep("y", nrow(yby))
# #'     by2 <- rbind(xby, yby)
# #'     if (verbose)
# #'         cat("\nby in x and y:\n")
# #'     if (verbose)
# #'         print(atotal(do.call("tab", by2), sum, "Total"))
# #'     nams <- union(names(x), names(y))
# #'     if (verbose)
# #'         print(c(DimX = dim(x), DimY = dim(y)))
# #'     if (verbose)
# #'         cat("\nVariables in both:\n")
# #'     if (verbose)
# #'         print(intersect(names(x), names(y)))
# #'     if (verbose)
# #'         cat("\nVariables in X only:\n")
# #'     if (verbose)
# #'         print(setdiff(names(x), names(y)))
# #'     if (verbose)
# #'         cat("\nVariables in Y only:\n")
# #'     if (verbose)
# #'         print(setdiff(names(y), names(x)))
# #'     x$FromX <- 1:nrow(x)
# #'     y$FromY <- 1:nrow(y)
# #'     mm <- merge(x, y, by, all = T, ...)
# #'     newroots <- setdiff(intersect(names(x), names(y)), by)
# #'     if (verbose)
# #'         cat("\nDimension of merged data frames:\n")
# #'     if (verbose)
# #'         print(c(DimMerge = dim(mm)))
# #'     if (verbose)
# #'         cat("\nNames of variables in merged data frame:\n")
# #'     if (verbose)
# #'         print(names(mm))
# #'     if (F) {
# #'         dotx <- grep("\.x", names(mm), value = T)
# #'         if (verbose)
# #'             print(c(dotx = dotx))
# #'         doty <- grep("\.y", names(mm), value = T)
# #'         if (verbose)
# #'             print(c(doty = doty))
# #'         rootx <- substring(dotx, 1, nchar(dotx) - 2)
# #'         rooty <- substring(doty, 1, nchar(doty) - 2)
# #'         newroots <- intersect(rootx, rooty)
# #'     }
# #'     FromBoth <- !is.na(mm$FromX) & !is.na(mm$FromY)
# #'     Xonly <- !is.na(mm$FromX) & is.na(mm$FromY)
# #'     Yonly <- is.na(mm$FromX) & !is.na(mm$FromY)
# #'     if (verbose)
# #'         cat("\nRows in:\n")
# #'     if (verbose)
# #'         print(c(Both = sum(FromBoth), Xonly = sum(Xonly), Yonly = sum(Yonly)))
# #'     if (verbose)
# #'         cat("\nThe following variables occur in both data frames:\n")
# #'     if (verbose)
# #'         print(newroots)
# #'     drop.list <- character(0)
# #'     for (nn in newroots) {
# #'         nn.x <- paste(nn, ".x", sep = "")
# #'         nn.y <- paste(nn, ".y", sep = "")
# #'         mm[[nn]] <- xm(mm[[nn.x]], mm[[nn.y]])
# #'         if (all(same <- consistent(mm[[nn.x]], mm[[nn.y]]))) {
# #'             if (verbose)
# #'                 cat("Variable ", nn, " is consistent\n")
# #'             drop.list <- c(drop.list, nn)
# #'         }
# #'         else {
# #'             if (verbose)
# #'                 cat("Variable ", nn, " is inconsistent in the following rows:\n")
# #'             if (verbose)
# #'                 print(mm[same, c(by, nn.x, nn.y, nn)])
# #'         }
# #'     }
# #'     if (dropdots)
# #'         drop.list <- newroots
# #'     drop <- if (length(drop.list) > 0) {
# #'         c(paste(drop.list, "x", sep = "."), paste(drop.list,
# #'             "y", sep = "."))
# #'     }
# #'     else character(0)
# #'     if (verbose)
# #'         cat("\nDrop list:\n")
# #'     if (verbose)
# #'         print(drop)
# #'     if (length(drop) > 0) {
# #'         if (verbose)
# #'             print(c(drop = drop))
# #'         mm <- mm[, -match(drop, names(mm))]
# #'     }
# #'     onams <- 1:length(nams)
# #'     onams <- c(onams, onams + 0.1, onams + 0.2)
# #'     names(onams) <- c(nams, paste(nams, ".x", sep = ""), paste(nams,
# #'         ".y", sep = ""))
# #'     keep <- intersect(names(sort(onams)), names(mm))
# #'     mm[, keep]
# #'   }
# #'
# #' @export
# xmerge <- function(x, y, by , all = TRUE, dropdots = FALSE , verbose = FALSE, debug = TRUE, from = FALSE, ... ) {
#     help <-
# "This is a modification of merge that combines consistent variables
# even if not specified in 'by' to keep a common name.
# -- Some errors fixed Apr 24, 2007"
#
#     xm <- function( a, b, tofac = is.factor(a)||is.factor(b)) {
#           # replace values of a with non-missing values of b
#           if ( tofac ) {
#               levs <- union( levels(b), levels(a))
#               a <- as.character(a)
#               b <- as.character(b)
#           }
#           b [ is.na(b) ] <- a[is.na(b)]
#
#           if ( tofac ) {
# 				levs <- union( levs, unique(b))
# 				b <- factor(b,levels = levs)
# 			}
#           b
#     }
#     na2f <- function(x) {
#         x[is.na(x)] <- F
#         x
#     }
#     consistent <- function(a,b) {
#             # which values are consistent if neither is NA
#             if( is.factor(a)) a <- as.character(a)
#             if( is.factor(b)) b <- as.character(b)
#             !na2f(a != b)
#     }
#     if ( from ) {
#         xname <- deparse(substitute(x))
#         yname <- deparse(substitute(y))
#         # x[[paste("F.",xname,sep="")]] <- rep(T, nrow(x))
#         # y[[paste("F.",yname,sep="")]] <- rep(T, nrow(y))
#         x[[".F"]] <- rep("x", nrow(x))
#         y[[".F"]] <- rep("y", nrow(y))
#     }
# 	## ids in each file
# 	xby <- x[,by,drop=F]
# 	yby <- y[,by,drop=F]
# 	xby$.file <- rep('x',nrow(xby))
# 	yby$.file <- rep('y',nrow(yby))
# 	by2 <- rbind( xby, yby)
# 	if ( verbose ) cat("\nby in x and y:\n")
# 	if ( verbose ) print(atotal(do.call("tab",by2),sum,"Total"))
#     nams <- union( names(x), names(y))
#    	if ( verbose ) print( c( DimX = dim(x), DimY = dim(y)))
#     if ( verbose ) cat("\nVariables in both:\n")
#     if ( verbose ) print(intersect( names(x), names(y)))
#     if ( verbose ) cat("\nVariables in X only:\n")
#     if ( verbose ) print( setdiff( names(x), names(y)))
#     if ( verbose ) cat("\nVariables in Y only:\n")
#     if ( verbose ) print( setdiff( names(y), names(x)))
#     # compare two data frames, e.g. for merging
#     x$FromX <-  1:nrow(x)
#     y$FromY <- 1:nrow(y)
#     mm <- merge( x, y, by, all = TRUE, ...)
#     # names in both data frames
#     newroots <- setdiff( intersect(names(x),names(y)), by)
#
#     if ( verbose ) cat("\nDimension of merged data frames:\n")
#     if ( verbose ) print(c("DimMerge"=dim(mm)))
#     if ( verbose ) cat("\nNames of variables in merged data frame:\n")
#     if ( verbose ) print(names(mm))
#
#     # Are similarly named variables consistent
#
#     if(F){
#         dotx <- grep("\\.x", names(mm), value =T)
#         if ( verbose )print( c(dotx = dotx))
#         doty <- grep("\\.y", names(mm),value = TRUE)
#         if ( verbose )print( c(doty = doty))
#         rootx <- substring( dotx, 1, nchar(dotx) - 2 )
#         rooty <- substring( doty, 1, nchar(doty) - 2 )
#         newroots <- intersect( rootx, rooty )
#     }
#     FromBoth <- !is.na(mm$FromX) & !is.na(mm$FromY)
#     Xonly <-  !is.na(mm$FromX) & is.na(mm$FromY)
#     Yonly <-  is.na(mm$FromX) & !is.na(mm$FromY)
#     if ( verbose )cat("\nRows in:\n")
#     if ( verbose )print( c( Both = sum(FromBoth), Xonly = sum(Xonly), Yonly = sum(Yonly)))
#     if ( verbose )cat("\nThe following variables occur in both data frames:\n")
#     if ( verbose )print( newroots )
#     drop.list <- character(0)
#     for ( nn in newroots) {
#         nn.x <- paste(nn,".x", sep = '')
#         nn.y <- paste(nn,".y", sep = '')
#         mm[[nn]] <- xm( mm[[nn.x]], mm[[nn.y]])
#         # Note that consistent NAs should be okay here
#         if( all( same <- consistent(mm[[nn.x]],mm[[nn.y]] )) ) {
#             if ( verbose )cat ("Variable ",nn," is consistent\n")
#             drop.list <- c(drop.list, nn)
#         } else {
#             if ( verbose )cat("Variable ",nn," is inconsistent in the following rows:\n")
#             if ( verbose )print(mm[ same, c(by, nn.x, nn.y, nn)])
#         }
#     }
# 	if( dropdots ) drop.list <- newroots
#     drop <- if( length( drop.list)>0) {
#           c(paste(drop.list, "x", sep = "."), paste(drop.list, "y", sep = ".") )
#           } else character(0)
# 	if ( verbose )cat("\nDrop list:\n")
# 	if ( verbose )print( drop)
#     if( length(drop) > 0) {
#             if ( verbose )print( c(drop=drop))
#             mm <- mm[, - match( drop, names(mm))]
#     }
#     # nice ordering
#     onams <- 1:length(nams)
#     #print(onams)
#     onams <- c(onams, onams+.1, onams +.2)
#     #print(onams)
#     names(onams) <- c(nams, paste(nams,".x",sep=''), paste(nams,".y",sep =''))
#     keep <- intersect( names(sort(onams)), names(mm))
#     mm[, keep]
# }
#
#
# if ( F ) {  # test xmerge
#     d1 <- data.frame(id = 1:5, x = 1:5, d = 1:5, v1 = 1:5)
#     d2 <- data.frame(id = 3:7, x = 3:7, d = 1:5, v2 = 3:7)
#
#     # need to specify 'id'
#
#     xmerge(d1,d2,by='id',verbose=T) # right row but inconsistent variable 'd' is created with 2nd df replacing first -- no warning
#     xmerge(d1,d2,by='id',all=T)   # ditto
# }
#
#
#
#

# ###
# ###  Rbind to stack data and prediction frames
# ###
#
#
#
# #' Generic rbind that works on data frames
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, ...)
# #' UseMethod("Rbind")
# #'
# #' @export
# Rbind <- function(x, ...) UseMethod("Rbind")
#
#
#
# #' Rbind applied to data frames
# #'
# #'
# #'
# #'
# #'
# #' @param \dots
# #' @param vname
# #' @param oname
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (..., vname = ".which", oname = ".order")
# #' {
# #'     z <- list(...)
# #'     for (ii in 1:length(z)) {
# #'         z[[ii]] <- as.data.frame(z[[ii]], stringsAsFactors = FALSE)
# #'         z[[ii]][, vname] <- rep(ii, nr <- nrow(z[[ii]]))
# #'         z[[ii]][, oname] <- 1:nr
# #'     }
# #'     ret <- z[[1]]
# #'     for (ii in 2:length(z)) ret <- merge(ret, z[[ii]], all = TRUE,
# #'         sort = FALSE)
# #'     ret[order(ret[[vname]], ret[[oname]]), ]
# #'   }
# #'
# #' @export
# Rbind.data.frame <-
# function( ... ,vname = '.which',oname = ".order") {
#
#         ## stacks data frames adding a new variable .which to identify each source data frame
#         ## good way to combine data and prediction data frames to show data and fitted
#         ## values in a common plot
#         Unique <- function(x, ...) UseMethod("Unique")
#         Unique.factor <- function( x, ...) levels(x)
#         Unique.default <- function( x, ...) unique(x)
#         Maxp1 <- function( x ) {
#               mm <- max(c(0, as.numeric( as.character(x))), na.rm = TRUE)
#               ceiling( mm + 1)
#         }
#         z <- list(...)
#         cum.which <- numeric(0)
#         for ( ii in 1:length(z)) {
#                 ddi <- as.data.frame(z[[ii]],stringsAsFactors = FALSE )
#                 if ( vname %in% names(ddi)) cum.which <- c(cum.which, Unique(ddi[[vname]]))
#                 else {
#                      ddi[[vname]] <- rep(vv <- Maxp1(cum.which), nrow( ddi))
#                      cum.which <- c(cum.which, vv)
#                 }
#                 ddi[,oname] <- 1:nrow(ddi)
#                 z[[ii]] <- ddi
#         }
#         ret <- z[[1]]
#         if( length(z) > 1 )for ( ii in 2:length(z)) ret <- merge(ret,z[[ii]], all = TRUE, sort = FALSE)
#         ret [ order( ret[[vname]], ret[[oname]]),]
# }
#
#
#
# #' Rbind applied to lists
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param \dots
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, ...)
# #' {
# #'     ret <- Rbind.data.frame(x, ...)
# #'     as.list(ret)
# #'   }
# #'
# #' @export
# Rbind.list <- function(x,...) {
#     ret <- Rbind.data.frame(x, ...)
#     as.list(ret)
# }
#
#
#
#
# if (FALSE) {
#     df1 <- data.frame( a = c('z','y','b'), x = 1:3, y1 = 1:3)
#     df2 <- data.frame( a = c('a','d','b'), x = c("A","A","B"), y2 = 4:6)
#
#     mat1 <- cbind( x=1:3, y1 = 1:3)
#
#     Rbind( mat1, df2)
#
#     df1
#     df2
#     z <- Rbind( df1, df2)
#     z
#     levels(z$a)
#     levels(z$x)
#
#     list1 <- list( a = 1:4, b = letters[1:4])
#     list2 <- list( a = 5:6, b = letters[5:6], c = factor(c("A","B")))
#     Rbind(list1, list2)
#
# }
#
#
#
#
#
#
# ###
# ###   Reshaping data: function to call 'reshape'
# ###
#
# long <- function ( data , varying=NULL , sep.varying = "\\.", v.names = names(Varying),
#                  timevar = 'time', idvar = 'id', ids = 1:NROW(data),
#                  times = seq(length=length(Varying[[1]])),
#                  drop = NULL, new.row.names = NULL,
#                  split = list(regexp="\\.", include = FALSE),
#                  debug = FALSE){
#     help <- "
#     Example:   long(dwl, varying = c('wmrss', 'lm1tot','lm2tot'), sep='')
#     "
#            nn <- names(data)
#            if ( debug ) disp(nn)
#            Varying <- list()
#            for ( ii in 1:length(varying)){
#                if(length( prefix <- varying[[ii]]) == 1) {
#                           ns <- grep(paste("^",prefix,sep.varying,sep=""), nn,value = T)
#                           if (debug) disp(ns)
#                           Varying[[ii]] <- ns
#                           names(Varying)[ii] <- prefix
#                }  else {
#                     Varying[[ii]] <- varying[[ii]]
#                     names(Varying)[ii] <- names(varying)[ii]
#                }
#
#            }
#            #print( varying)
#            #print( v.names )
#            if ( debug ) disp(Varying)
#            if ( debug ) disp(times)
#     data <- data[,sort(names(data))] # for bug in reshape is time orders are different for different variables
#          ret <- stats::reshape( data, Varying, v.names , timevar, idvar, ids,
#                 times , drop, direction = 'long', new.row.names, sep.varying,
#                 split)
#          ret [ order( ret[[idvar]], ret[[timevar]]),]
#
# }
#
# # very simple reshape to long
# # tolong <- function(data, sep = "_",...){
# #  reshape(data, direction = 'long', sep = sep, varying = grep(sep, names(data)),...)
# #}
#
#
# # tolong <- function(data, sep = "_", expand = FALSE, safe_sep = "#>{{{{",...){
# #   #
# #   # This creates a long file from a wide file is each variable name for a 'varying variable'
# #   # has the form 'varname_time' and
# #   # 1.  the separator ('_' by default) does not occur elsewhere in any variable name
# #   # 2.  every possible 'varname' by 'time' combination occurs exactly once, i.e.
# #   #     every varying variable name exists for each possible time.
# #   # If expand == TRUE, then new variables are created to complete
# #   #     all possible combinations.
# #   # Note: there appears to be a bug in 'reshape' if the ordering of wide
# #   # variable names is not consistent wrt times. We attempt to address this
# #   # reordering the variables before calling reshape
# #   #
# #   if(expand) {
# #     namessafe <- sub(sep,safe_sep,names(data),fixed=TRUE)
# #     varnames <- grep( safe_sep, namessafe, value = TRUE, fixed = TRUE)
# #     names <- unique(sub(paste(safe_sep,".*$",sep=''),"",varnames))
# #     times <- unique(sub(paste("^.*",safe_sep,sep=''),"",varnames))
# #     allnames <- paste(rep(names,each=length(times)),sep,rep(times,length(names)),sep='')
# #     z <- data
# #     for ( nn in allnames ) {
# #       z[[nn]] <- if(is.null(data[[nn]])) NA else data[[nn]]
# #     }
# #     data <- z
# #   }
# #   data <- data[, sort(names(data))] # this seems to handle a bug in reshape
# #   # that appears to use the order of variables instead of the actual suffixes
# #   # to determine the time to which a value is ascribed.
# #   reshape(data, direction = 'long', sep = sep, varying = grep(sep, names(data),fixed=TRUE),...)
# # }
#
# #' @export
# tolong <- function(data, sep = "_", expand = FALSE, safe_sep = "#%@!",...){
#   #
#   # This creates a long file from a wide file is each variable name for a 'varying variable'
#   # has the form 'varname_time' and
#   # 1.  the separator ('_' by default) does not occur elsewhere in any variable name
#   # 2.  every possible 'varname' by 'time' combination occurs exactly once, i.e.
#   #     every varying variable name exists for each possible time.
#   #
#   # TO DO:
#   # - improve on safe_sep, e.g. generate a random string of safe characters
#   # then check whether it's in the names, if so expand length and repeat.
#   # remove 'safe_sep' from arguments.
#   # - Similarly do something smarter than ZZZZZ
#   # - warn if there is a variable named time
#   if( length( grep('^time$',names(data)))) warning("Variable 'time' in data is replaced by a variable to count occasions")
#   if(expand) {
#     namessafe <- sub(sep,safe_sep,names(data),fixed=TRUE)
#     varnames <- grep( safe_sep, namessafe, value = TRUE, fixed = TRUE)
#     names <- unique(sub(paste(safe_sep,".*$",sep=''),"",varnames))
#     times <- unique(sub(paste("^.*",safe_sep,sep=''),"",varnames))
#     allnames <- paste(rep(names,each=length(times)),sep,rep(times,length(names)),sep='')
#     z <- data
#     for ( nn in allnames ) {
#       z[[nn]] <- if(is.null(data[[nn]])) NA else data[[nn]]
#     }
#     data <- z
#   }
#   namessafe <- sub(sep,safe_sep,names(data),fixed=TRUE)
#   namestimes <- sub(paste("^.*",safe_sep,sep=''),"ZZZZZ",namessafe)
#   ord <- order(namestimes)
#   data <- data[, ord] # this seems to handle a bug in reshape
#   # that appears to use the order of variables instead of the actual suffixes
#   # to determine the time to which a value is ascribed.
#   reshape(data, direction = 'long', sep = sep, varying = grep(sep, names(data),fixed=TRUE),...)
# }
#
# #
# # zd <- data.frame( sub = c('a','b','c'), x.1 = 10+1:3, x.2 = 20+1:3, y.2 = c('a','b','c'), y.1 = factor(2:4),time=1:3, id = letters[1:3])
# # tolong( zd,sep='.')
# # tolong(zd, varying = list( y = c('y.1','y.2'), x = c("x.1","x.2")))
#
#
#
#

#
#
#
#
# ###
# ###  Trellis additions
# ###
#
#
# ###
# ###   Miscellaneous utility functions
# ###
#
#
#
#
# #' Transform selected values to NAs
# #'
# #'
# #'
# #'
# #'
# #' @param x
# #' @param val
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, val)
# #' {
# #'     x[match(x, val, 0) > 0] <- NA
# #'     x
# #'   }
# #'
# #' @export
# val2na <- function( x, val) {
# 	## val2na(1:10, c(3,5))
# 	## val2na(factor(c('a','b','b','c','a')), c('b','a'))
# 	x[match(x,val,0)>0] <- NA
# 	x
# }
#
# # zf <- factor(c('a','a','b',NA,'c','a',NA))
#
#
#
#

#
# # grepl <- function(pattern,x,...) match( 1:length(x) , grep(pattern,x,...), 0) >0
#
# #' @export
# p <- function(...) paste(..., sep ="")
#
# #' @export
# ch <- function(x) as.character(x)
#
# #' @export
# "%less%" <- function( a, b) setdiff(a,b)
# #' @export
# "%or%" <- function( a, b) union(a,b)
#
#
#
#
# #' Set operators
# #'
# #' Set operations written as binary operators
# #'
# #' @aliases %and% %less% %or%
# #' @param a,b vectors treated as sets
# #' @export
# "%and%" <- function( a, b) intersect( a, b)
#
#
#
# #' Library workaround for RDC
# #'
# #' @export
# lib <- function(x) {
#     xn <- deparse(substitute(x))
#     if ( !do.call("require", list(as.name(xn))) ) {
#         install.packages(xn)
#         do.call("library",list(as.name(xn)))
#     }
# }
#
# ###
# ###  Interfacing with SAS
# ###
#
#
#
#
# #' Read a SAS ODS file
# #'
# #' @param file input file
# #' @param tfile
# #' @export
# sasin <- function(file, tfile = tempfile() ) {
# # moved to "/R/fun.R"
#     help = "
#     sasin reads a .csv file created by SAS with
#        ODS CSV FILE = 'file';
#         < SAS procedure statements >
#        ODS CSV CLOSE;
#     The tables produced by SAS are elements in the list
#     returned by sasin.
#     "
#
#     todf <- function(ll) {
#         if ( length(ll) < 3) return (character(0))
#         if ( length(ll) == 3) return (ll[2])
#         cat(ll[2],"\n" , file = tfile)
#         for ( ii in 3:(length(ll)-1)) {
#             cat(ll[ii], "\n",file = tfile ,append = TRUE)
#         }
#         df <- read.csv( tfile , header = F)
#         if ( !any ( sapply( df, is.numeric ))) df <- read.csv(tfile)
#         df
#     }
#     readin <- scan(file,what="",sep="\n",blank.lines.skip=F)
#     blanks <- which(readin == "")
#     head.pos <- c(1,1+head(blanks,-1))
#     heads  <- gsub('\\"|,',"",readin[head.pos])
#     # disp(heads)
#     reps   <- diff( c(head.pos, 1+length(readin)))
#     # disp(reps)
#     heads  <- rep(heads, reps)
#     readin <- split( readin, heads)
#     readin <- lapply( readin , todf)
#     readin
# }
#
#
#
#
# #########
# #########   SCRAPS
# #########
#
# ###
# ### Stack
# ###
#
# ## NOTE THAT R HAS A FUNCTION CALLED stack
# if (FALSE) {
# stack <- function(...) {
# 	# simple version of stack, uses names of first df.
#       # 05-04-19: THis function might be obsolete d/t mergec below
# 	ll <- list(...)
# 	llclass <- sapply(ll,function(x) inherits(x,'data.frame'))
# 	if(!all(llclass)){
# 		dfr <- ll[[1]]
# 		keep <- ll[[2]]
# 		add <- ll[[3]]
# 		ll <- list(0)
# 		for ( i in 1:length(add) ) ll[[i]] <- dfr[,c(keep,add[i])]
# 	}
# 	nam <- names(ll[[1]])
# 	for ( ii in 2:length(ll)) names(ll[[ii]]) <- nam
# 	ret <- do.call('rbind',ll)
# 	nrows <- sapply(ll, function(x) dim(x)[1])
# 	ret$Type. <- rep(1:length(ll),nrows)
# 	ret
# }
# }
#
#
#
# #####
# #####  anova.lme
# #####
#
# #' @export
# xanova.lme <- function (object, ..., test = TRUE, type = c("sequential", "marginal"),
#     adjustSigma = TRUE, Terms, L, verbose = FALSE)
# {
#     warning("This is a modified version of anova.lme that uses min dfs for the denominator")
#     Lmiss <- missing(L)
#     dots <- list(...)
#     if ((rt <- (length(dots) + 1)) == 1) {
#         if (!inherits(object, "lme")) {
#             stop("Object must inherit from class \"lme\" ")
#         }
#         vFix <- attr(object$fixDF, "varFixFact")
#         if (object$method == "ML" && adjustSigma == TRUE) {
#             vFix <- sqrt(object$dims$N/(object$dims$N - ncol(vFix))) *
#                 vFix
#         }
#         c0 <- solve(t(vFix), fixef(object))
#         assign <- attr(object$fixDF, "assign")
#         nTerms <- length(assign)
#         if (missing(Terms) && Lmiss) {
#             type <- match.arg(type)
#             Fval <- Pval <- double(nTerms)
#             nDF <- integer(nTerms)
#             dDF <- object$fixDF$terms
#             for (i in 1:nTerms) {
#                 nDF[i] <- length(assign[[i]])
#                 if (type == "sequential") {
#                   c0i <- c0[assign[[i]]]
#                 }
#                 else {
#                   c0i <- c(qr.qty(qr(vFix[, assign[[i]], drop = FALSE]),
#                     c0))[1:nDF[i]]
#                 }
#                 Fval[i] <- sum(c0i^2)/nDF[i]
#                 Pval[i] <- 1 - pf(Fval[i], nDF[i], dDF[i])
#             }
#             aod <- data.frame(nDF, dDF, Fval, Pval)
#             dimnames(aod) <- list(names(assign), c("numDF", "denDF",
#                 "F-value", "p-value"))
#             attr(aod, "rt") <- rt
#         }
#         else {
#             nX <- length(unlist(assign))
#             if (Lmiss) {
#                 if (is.numeric(Terms) && all(Terms == as.integer(Terms))) {
#                   if (min(Terms) < 1 || max(Terms) > nTerms) {
#                     stop(paste("Terms must be between 1 and",
#                       nTerms))
#                   }
#                 }
#                 else {
#                   if (is.character(Terms)) {
#                     if (any(noMatch <- is.na(match(Terms, names(assign))))) {
#                       stop(paste("Term(s)", paste(Terms[noMatch],
#                         collapse = ", "), "not matched"))
#                     }
#                   }
#                   else {
#                     stop("Terms can only be integers or characters")
#                   }
#                 }
#                 dDF <- unique(object$fixDF$terms[Terms])
#                 if (length(dDF) > 1) {
#                   # stop("Terms must all have the same denominator DF")
#                   warning("Terms do not all have the same denominator DF -- using the minimum")
#                   dDF <- min(dDF)
#                 }
#                 lab <- paste("F-test for:", paste(names(assign[Terms]),
#                   collapse = ", "), "\n")
#                 L <- diag(nX)[unlist(assign[Terms]), , drop = FALSE]
#             }
#             else {
#                 L <- as.matrix(L)
#                 if (ncol(L) == 1)
#                   L <- t(L)
#                 nrowL <- nrow(L)
#                 ncolL <- ncol(L)
#                 if (ncol(L) > nX) {
#                   stop(paste("L must have at most", nX, "columns"))
#                 }
#                 dmsL1 <- rownames(L)
#                 L0 <- array(0, c(nrowL, nX), list(NULL, names(object$fixDF$X)))
#                 if (is.null(dmsL2 <- colnames(L))) {
#                   L0[, 1:ncolL] <- L
#                 }
#                 else {
#                   if (any(noMatch <- is.na(match(dmsL2, colnames(L0))))) {
#                     stop(paste("Effects", paste(dmsL2[noMatch],
#                       collapse = ", "), "not matched"))
#                   }
#                   L0[, dmsL2] <- L
#                 }
#                 L <- L0[noZeroRowL <- as.logical((L0 != 0) %*%
#                   rep(1, nX)), , drop = FALSE]
#                 nrowL <- nrow(L)
#                 if (is.null(dmsL1)) {
#                   dmsL1 <- 1:nrowL
#                 }
#                 else {
#                   dmsL1 <- dmsL1[noZeroRowL]
#                 }
#                 rownames(L) <- dmsL1
#                 dDF <- unique(object$fixDF$X[noZeroColL <- as.logical(c(rep(1,
#                   nrowL) %*% (L != 0)))])
#                 if (length(dDF) > 1) {
#                   ## stop("L may only involve fixed effects with the same denominator DF")
#                   warn <- paste( "L involves fixed effects with the different denominator DF:",
#                           paste(dDF, collapse=" "), collapse = " ")
#                   warning(warn)
#                   dDF <- min(dDF)
#                 }
#                 lab <- "F-test for linear combination(s)\n"
#             }
#             nDF <- sum(svd(L)$d > 0)
#             c0 <- c(qr.qty(qr(vFix %*% t(L)), c0))[1:nDF]
#             Fval <- sum(c0^2)/nDF
#             Pval <- 1 - pf(Fval, nDF, dDF)
#             aod <- data.frame(nDF, dDF, Fval, Pval)
#             names(aod) <- c("numDF", "denDF", "F-value", "p-value")
#             attr(aod, "rt") <- rt
#             attr(aod, "label") <- lab
#             if (!Lmiss) {
#                 if (nrow(L) > 1)
#                   attr(aod, "L") <- L[, noZeroColL, drop = FALSE]
#                 else attr(aod, "L") <- L[, noZeroColL]
#             }
#         }
#     }
#     else {
#         ancall <- sys.call()
#         ancall$verbose <- ancall$test <- NULL
#         object <- list(object, ...)
#         termsClass <- unlist(lapply(object, data.class))
#         if (!all(match(termsClass, c("gls", "gnls", "lm", "lmList",
#             "lme", "nlme", "nlsList", "nls"), 0))) {
#             stop(paste("Objects must inherit from classes \"gls\", \"gnls\"",
#                 "\"lm\",\"lmList\", \"lme\",\"nlme\",\"nlsList\", or \"nls\""))
#         }
#         resp <- unlist(lapply(object, function(el) deparse(getResponseFormula(el)[[2]])))
#         subs <- as.logical(match(resp, resp[1], FALSE))
#         if (!all(subs))
#             warning(paste("Some fitted objects deleted because",
#                 "response differs from the first model"))
#         if (sum(subs) == 1)
#             stop("First model has a different response from the rest")
#         object <- object[subs]
#         rt <- length(object)
#         termsModel <- lapply(object, function(el) formula(el)[-2])
#         estMeth <- unlist(lapply(object, function(el) {
#             val <- el[["method"]]
#             if (is.null(val))
#                 val <- NA
#             val
#         }))
#         if (length(uEst <- unique(estMeth[!is.na(estMeth)])) >
#             1) {
#             stop("All fitted objects must have the same estimation method.")
#         }
#         estMeth[is.na(estMeth)] <- uEst
#         REML <- uEst == "REML"
#         if (REML) {
#             aux <- unlist(lapply(termsModel, function(el) {
#                 aux <- terms(el)
#                 val <- paste(sort(attr(aux, "term.labels")),
#                   collapse = "&")
#                 if (attr(aux, "intercept") == 1) {
#                   val <- paste(val, "(Intercept)", sep = "&")
#                 }
#                 val
#             }))
#             if (length(unique(aux)) > 1) {
#                 warning(paste("Fitted objects with different fixed effects.",
#                   "REML comparisons are not meaningful."))
#             }
#         }
#         termsCall <- lapply(object, function(el) {
#             if (is.null(val <- el$call)) {
#                 if (is.null(val <- attr(el, "call"))) {
#                   stop("Objects must have a \"call\" component or attribute.")
#                 }
#             }
#             val
#         })
#         termsCall <- unlist(lapply(termsCall, function(el) paste(deparse(el),
#             collapse = "")))
#         aux <- lapply(object, logLik, REML)
#         if (length(unique(unlist(lapply(aux, function(el) attr(el,
#             "nall"))))) > 1) {
#             stop("All fitted objects must use the same number of observations")
#         }
#         dfModel <- unlist(lapply(aux, function(el) attr(el, "df")))
#         logLik <- unlist(lapply(aux, function(el) c(el)))
#         AIC <- unlist(lapply(aux, AIC))
#         BIC <- unlist(lapply(aux, BIC))
#         aod <- data.frame(call = termsCall, Model = (1:rt), df = dfModel,
#             AIC = AIC, BIC = BIC, logLik = logLik, check.names = FALSE)
#         if (test) {
#             ddf <- diff(dfModel)
#             if (sum(abs(ddf)) > 0) {
#                 effects <- rep("", rt)
#                 for (i in 2:rt) {
#                   if (ddf[i - 1] != 0) {
#                     effects[i] <- paste(i - 1, i, sep = " vs ")
#                   }
#                 }
#                 pval <- rep(NA, rt - 1)
#                 ldf <- as.logical(ddf)
#                 lratio <- 2 * abs(diff(logLik))
#                 lratio[!ldf] <- NA
#                 pval[ldf] <- 1 - pchisq(lratio[ldf], abs(ddf[ldf]))
#                 aod <- data.frame(aod, Test = effects, L.Ratio = c(NA,
#                   lratio), "p-value" = c(NA, pval), check.names = FALSE)
#             }
#         }
#         row.names(aod) <- unlist(lapply(as.list(ancall[-1]),
#             deparse))
#         attr(aod, "rt") <- rt
#         attr(aod, "verbose") <- verbose
#     }
#     class(aod) <- c("anova.lme", "data.frame")
#     aod
# }
#
#
#
#
#
#
#
#
#
#
# if( FALSE ) {
#
#
# #### NOT RUN: ####
#
#
# glmmPQL <-
# function (fixed, random, family, data, correlation, weights,
#     control, niter = 10, verbose = TRUE, ...)
# {
#     if (!require("nlme"))
#         stop("package 'nlme' is essential")
#     if (is.character(family))
#         family <- get(family)
#     if (is.function(family))
#         family <- family()
#     if (is.null(family$family)) {
#         print(family)
#         stop("'family' not recognized")
#     }
#     m <- mcall <- Call <- match.call()
#     nm <- names(m)[-1]
#     keep <- is.element(nm, c("weights", "data", "subset", "na.action"))
#     for (i in nm[!keep]) m[[i]] <- NULL
#     allvars <- if (is.list(random))
#         allvars <- c(all.vars(fixed), names(random), unlist(lapply(random,
#             function(x) all.vars(formula(x)))))
#     else c(all.vars(fixed), all.vars(random))
#     Terms <- if (missing(data))
#         terms(fixed)
#     else terms(fixed, data = data)
#     off <- attr(Terms, "offset")
#     if (length(off <- attr(Terms, "offset")))
#         allvars <- c(allvars, as.character(attr(Terms, "variables"))[off +
#             1])
#     m$formula <- as.formula(paste("~", paste(allvars, collapse = "+")))
#     environment(m$formula) <- environment(fixed)
#     m$drop.unused.levels <- TRUE
#     m[[1]] <- as.name("model.frame")
#     mf <- eval.parent(m)
#     off <- model.offset(mf)
#     if (is.null(off))
#         off <- 0
#     w <- model.weights(mf)
#     if (is.null(w))
#         w <- rep(1, nrow(mf))
#     mf$wts <- w
#     fit0 <- glm(formula = fixed, family = family, data = mf,
#         weights = wts, ...)
#     w <- fit0$prior.weights
#     eta <- fit0$linear.predictor
#     zz <- eta + fit0$residuals - off
#     wz <- fit0$weights
#     fam <- family
#     nm <- names(mcall)[-1]
#     keep <- is.element(nm, c("fixed", "random", "data", "subset",
#         "na.action", "control"))
#     for (i in nm[!keep]) mcall[[i]] <- NULL
#     fixed[[2]] <- quote(zz)
#     mcall[["fixed"]] <- fixed
#     mcall[[1]] <- as.name("lme")
#     mcall$random <- random
#     mcall$method <- "ML"
#     if (!missing(correlation))
#         mcall$correlation <- correlation
#     mcall$weights <- quote(varFixed(~invwt))
#     mf$zz <- zz
#     mf$invwt <- 1/wz
#     mcall$data <- mf
#     for (i in 1:niter) {
#         if (verbose) {
#             cat("iteration", i, "\n")
#             print(names(mcall))
#         }
#         fit <- eval(mcall)
#         etaold <- eta
#         eta <- fitted(fit) + off
#         if (sum((eta - etaold)^2) < 1e-06 * sum(eta^2))
#             break
#         mu <- fam$linkinv(eta)
#         mu.eta.val <- fam$mu.eta(eta)
#         mf$zz <- eta + (fit0$y - mu)/mu.eta.val - off
#         wz <- w * mu.eta.val^2/fam$variance(mu)
#         mf$invwt <- 1/wz
#         mcall$data <- mf
#     }
#     attributes(fit$logLik) <- NULL
#     fit$call <- Call
#     fit$family <- family
#     oldClass(fit) <- c("glmmPQL", oldClass(fit))
#     fit
# }
#
# #log <-function (x, base = exp(1)) {
# #    x <- pmax(x,.00000001)
# #    if (missing(base)) .Internal(log(x)) else .Internal(log(x, base))
# #}
#
#
#
#
# logLik.reStruct  <-
# function (object, conLin, ...)
# {
#     if (any(!is.finite(conLin$Xy)))
#         return(-1000000)
#     .C("mixed_loglik", as.double(conLin$Xy), as.integer(unlist(conLin$dims)),
#         as.double(pdFactor(object)), as.integer(attr(object,
#             "settings")), loglik = double(1), double(1), PACKAGE = "nlme")$loglik
# }
#
# }
#
#
# ##
# ##
# ##  General polynomial splines
# ##
# ##
# ##
#
#
# ##
# ##
# ##  fun-new.R   August 8, 2008
# ##
# ##
#
#
#
#
# oldtab = function(...) {
#     UseMethod("tab")
# }
#
# oldtab.default = function( ..., total.margins = TRUE ) {
#     # This might -- finally -- be equivalent to table(..., exclude = NULL)
#         aa <- list(...)
#         if ( length(aa) == 1 && is.list ( aa[[1]]) ) {
#              return( do.call("tab", aa[[1]]))
#     }
#         #disp( names(aa))
#         #disp( is.null( names(aa)))
#         if ( is.null(names(aa))) {
#             nns = names( match.call())
#             #disp( nns)
#             names(aa) = nns[ 2:(1+length(aa)) ]
#         }
#         for ( ii in 1:length(aa) ) aa[[ii]] <- factor(aa[[ii]], exclude = NULL)
#         #disp( "aa" )
#         #disp( aa)
#         ret <- do.call("table",aa)
#         if ( total.margins) ret = atotal(ret)
#         ret
# }
#
# #' @export
# oldtab.data.frame = function( dd, fmla , total.margins = TRUE) {
#      if ( missing( fmla )) return ( do.call('tab', as.list( dd , total.margins = total.margins)))
#      xx = model.frame( fmla , dd , na.action = na.include)
#      xx = c(xx, total.margins = total.margins)
#      do.call( 'tab', xx)
# }
#
# #' @export
# oldtab.formula <- function( fmla, dd, total.margins = TRUE, ...){
#                tab( dd, fmla, total.margins = total.margins, ...)
# }
#
# if (FALSE) {
#
#     tt = function ( ... ) {
#           disp( list(...))
#           disp( nargs())
#           disp( names(match.call()))
#           disp( deparse( match.call()) )
#     }
#
#     tt( 1:3, b = 3:4)
#     tt( a = 1:3, b = 3:4)
#
#
#     zdd = data.frame( x = c(1,2,3,2,1,2,NA,1,2,NaN,Inf), A = c(NA,rep( c('a','b'), 5)), B = c(rep( c('x','z'),each = 5),NA))
#     tab( zdd, ~ x + A + B)
#     with( rima, tab( Category, Year))
#     tab( rima[c('Category','Year')])
# }
#
#
# ###
# ###
# ###   ith derivative of an expression
# ###
# ###
# ###
# ###
#
# #' @export
# d <- function( ex, xv,i ) if ( i == 0) ex else D( d(ex, xv, i-1), xv)
#
# # Example
# if (FALSE) {
#    d( quote( exp( x^2 /2 )), 'x', 3)      # admittedly not efficient code
#    d( quote( 3*x^3 + a*x^2), 'x', 2)
# }
#
#
# # reorder factor removed because functional version now in base package
# #
# # reorder.factor <- function (x, v, FUN = mean, ...) {
# #
# #       # Hmisc returns an ordered factor,
# #       # This returns the same as its input
# #       if ( inherits(x,'ordered')) ordered(x, levels(x)[order(tapply(v, x, FUN, ...))])
# #       else factor(x, levels(x)[order(tapply(v, x, FUN, ...))])
# # }
# #
# #
#
# # misscode added May 20, 2012
# # recodes NAs as a value below the range of observed
# # values to produce 3-D marginal value plots
#
# #' @export
# misscode <- function(x,...) UseMethod('misscode')
#
#
# #' Turns NAs into a value below the range of non-missing data for plotting %%
# #' ~~function to do ... ~~
# #'
# #' \code{misscode} turns NAs of numerical variables into a value slightly less
# #' that non-missing values. Missing values for factors are made into a
# #' non-missing level. When applied to a data frame, a new variable
# #' \code{.nmiss} is added to the data frame indicating the number of variables
# #' with missing data in each row of the data frame.
# #'
# #'
# #'
# #' @aliases misscode.default misscode.data.frame misscode.factor misscode
# #' @param x
# #' @param \dots
# #' @param offset
#
# #' @author Georges Monette
#
#
# #' @keywords ~kwd1 ~kwd2
# #' @examples
# #'
# #' ##---- Should be DIRECTLY executable !! ----
# #' ##-- ==>  Define data, use random,
# #' ##--	or do  help(data=index)  for the standard data sets.
# #'
# #' ## The function is currently defined as
# #' function (x, ..., offset = 0.1)
# #' {
# #'     rr <- range(x, na.rm = TRUE)
# #'     vmiss <- min(x, na.rm = TRUE) - offset * diff(rr)
# #'     nas <- is.na(x)
# #'     x[nas] <- vmiss
# #'     attr(x, "nas") <- nas
# #'     x
# #'   }
# #'
# #' @export
# misscode.default <- function(x,...,offset = .1) {
#   rr <- range(x, na.rm = TRUE)
#   vmiss <- min(x,na.rm = TRUE) - offset * diff(rr)
#   nas <- is.na(x)
#   x[nas] <- vmiss
#   attr(x,'nas') <- nas
#   x
# }
# #' @export
# misscode.factor <- function(x, ...) {
#   nas <- is.na(x)
#   x <- addNA(x, ifany = TRUE)
#   attr(x,'nas') <- nas
#   x
# }
# #' @export
# misscode.data.frame <- function(x,...) {
#   x[] <- lapply(x[],misscode,...)
#   isna <- lapply( x, function(x) attr(x,'nas'))
#   #   disp(isna)
#   isna <- do.call(cbind,isna)
#   isna <- apply(isna, 1, sum)
#   #   disp(isna)
#   x$.nmiss <- isna
#   x
# }
#
#
#
#
#
#
