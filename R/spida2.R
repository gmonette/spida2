#' spida2: Functions used in the Statistical Consulting Service at York University
#' 
#' The spida2 package is a collection of functions and datasets intended primarily
#' for statistical consulting, with a particular emphasis on
#' longitudinal and hierarchical data analysis. Some documents related to this 
#' package can be found at \url{http://blackwell.math.yorku.ca/R/spida2/doc}
#'
#' @section General Bugs:
#' \itemize{
#' \item \code{\link[latticeExtra]{+.trellis}} in xyplot + xyplot requires same data set to avoid misregistration of panel if some levels of panel factors are absent.
#'         Use \code{\link{Rbind}} to create common data frame and different names for plotted variables, e.g. to plot both lines and points. 
#' }
#' 
#' @section Code Trick:
#' \itemize{
#' \item To create a unique name: 
#' \code{
#' id <- '.id'
#' while( any(id %in% nams)) id <- paste0(id,".")
#' }
#' }
#' 
#' @section New functions:
#' \itemize{
#' \item \code{\link[spida2]{Rbind}} works on data.frame, missing variables get NA: use to avoid shifted panels in latticeExtra
#' \item \code{\link[spida2]{Apply}} returns list with same structure, e.g. list array
#' \item \code{\link[spida2]{tabf}} returns result of applying function as a list array
#' \item \code{\link[spida2]{waldx}} temporary name for rank-deficient aware version of wald
#' \item \code{\link[spida2]{waldf}} rank-deficient aware version of wald that returns a data frame and L matrix
#' \item \code{\link[spida2]{subrow}} subtracts selected rows of L matrix from ranges of rows, e.g. to compare with comparators
#' \item \code{\link[spida2]{lchol}} returns lower-triangular L so that G = L'L.
#' \item \code{\link[spida2]{getR}}, \code{\link[spida2]{getG}} and \code{\link[spida2]{getV}} return R, G and V matrices for \code{\link[nlme]{lme}} objects.
#' \item \code{\link[spida2]{pdInd}} constructs a pdClass for a G matrix with patterns of zero covariances. See
#' a vignette at \href{http://blackwell.math.yorku.ca/R/spida2/doc/pdInd.html}{pdInd: G matrix with pattern of zeros}.
#' }
#' 
#' @section Wald tests and linear hypothesis matrices:
#' \itemize{
#' \item \code{\link[spida2]{wald}} Wald tests with L matrices optionally created with regular expressions. Uses SVD to handle linear dependencies in rows of L
#' \item \code{\link[spida2]{walddf}} version of wald that returns a data frame
#' \item \code{\link[spida2]{as.data.frame.wald}} return a data frame from a wald object
#' \item \code{\link[spida2]{coef.wald}} method to extract estimated coefficients
#' \item \code{\link[spida2]{print.wald}} printing method
#' \item \code{\link[spida2]{Lfx}}  creates hypothesis matrices for derivatives and differences. Add example for factor differences
#' \item \code{\link[spida2]{M}} constructor for M objects to generate portions of design and hypothesis matrices. Used with  \code{\link[spida2]{Lfx}}
#' \item \code{\link[spida2]{rpfmt}} format estimated values and p-values from a wald test
#' \item \code{\link[spida2]{Lall}} for lmer objects
#' \item \code{\link[spida2]{Lc}} for lmer objects
#' \item \code{\link[spida2]{Lmu}} for lmer objects
#' }
#' @section Utilities for fitted objects:
#' \itemize{
#' \item \code{\link[spida2]{getD}} get data frame from a fitted object
#' \item \code{\link[spida2]{getData}} older version using methods. Might work if previous fails
#' \item \code{\link[spida2]{getFix}} get fixed effects from a fitted object
#' \item \code{\link[spida2]{getX}} get X matrix from fitted object
#' \item \code{\link[spida2]{getV}} get V matrix from a mixed model
#' \item \code{\link[spida2]{getG}} get G matrix from a mixed model
#' \item \code{\link[spida2]{getR}} get R matrix from a mixed model
#' \item \code{\link[spida2]{Vcov}} get estimated variance covariance of fixed effects from a fitted object
#' }
#' @section Multilevel data frames:
#' \itemize{
#' \item \code{\link[spida2]{capply}}: \code{capply(x,id,FUN)} applies the
#'   function \code{FUN} to chunks of 'x' formed by levels of 'id'.
#'   The result has the same form as 'x' with replication within
#'   chunks, if needed.
#' \item \code{\link{up}}, \code{\link[spida2]{agg}} and \code{\link[spida2]{up_apply}} create summary
#'   data sets consisting, by default, of within-id-invariant variables. 
#'   Summaries of id-varying variables can also be included. \code{\link[spida2]{agg}} can
#'   create mean incidence matrices for lower-level factors.
#' \item \code{\link[spida2]{cvar}} and \code{\link[spida2]{dvar}} are designed to be used
#'   in linear model formulas to generate 'centered-within-group' and 'within-group
#'   deviation' variables. WIth factors, they generate mean incidence matrices.
#' \item \code{link[spida2]{varLevel}} and \item \code{link[spida2]{gicc}}: the level of a variable with respect
#'   to a clustering formula and the 'generalized' intra-class correlation coefficient.
#' \item \code{\link[spida2]{tolong}} and \code{\link[spida2]{towide}} are 
#'   interfaces to \code{\link[stats]{reshape}} to facilitate
#'   the typical uses of reshape for longitudinal data.
#' }
#' @section Splines -- parametric and non-parametric:
#' \itemize{
#' \item \code{\link[spida2]{gsp}} creates a function for a generalized 
#'       spline that can
#'       be included in a linear model formula. 
#' \item \code{\link[spida2]{sc}} creates a spline contrast matrix for 
#'       general spline hypotheses. The matrix can be included
#'       in hypothesis matrices for the \code{\link[spida2]{wald}} function.
#' \item \code{\link[spida2]{smsp}} creates a matrix for a smoothing spline.
#' }
#' @section Datasets:
#' \itemize{
#' \item \code{\link[spida2]{hsfull}} Classical data set on high school math achievement and ses. See Bryk and Raudenbush and many other sources.
#' \item \code{\link[spida2]{iq}} Recovery after traumatic brain injury
#' \item \code{\link[spida2]{Drugs}} Longitudinal data on drugs and schizophrenia symptoms. Illustrates role of control variables with non-random assignment. 
#' \item \code{\link[spida2]{Indonesia}} Xerophthalmia.
#' \item \code{\link[spida2]{migraines}} Longitudinal data on migraine treatment and weather
#' \item \code{\link[spida2]{coffee}} Artificial data on coffee, heart damage and stress. Illustrates Simpson's Paradox with continuous predictors.
#' \item \code{\link[spida2]{hw}} Artificial data on height, weight and health. Illustrates suppression.
#' \item \code{\link[spida2]{Unemp}} U.S. monthly unemployment from January 1995 to February 2019.
#' }
#' @section Graphics:
#' \itemize{
#' \item \code{\link[spida2]{gd}} and \code{\link[spida2]{td}} easy interface to set
#'       graphical parameters for lattice and graphics. \code{gd} sets
#'       parameters to make graphs look like ggplot2 graphics. 
#' \item \code{\link[spida2]{panel.fit}} add fitted values and error bands with \code{\link[latticeExtra]{layer}} or \code{\link[latticeExtra]{glayer}}
#' \item \code{\link[spida2]{panel.dell}} add data ellipse \code{\link[latticeExtra]{layer}} or \code{\link[latticeExtra]{glayer}}
#' }
#' @section Miscellaneous utility functions:
#' \itemize{
#' \item \code{\link[spida2]{sortdf}} sort rows of a data frame -- useful in a magrittr pipeline
#' \item \code{\link[spida2]{assn}} assign -- useful in a magrittr pipeline
#' \item \code{\link[spida2]{disp}} utility to display value of a variable -- useful for debugging
#' \item \code{\link[spida2]{getFactorNames}} get names of variables that are factors in a data frame
#' \item \code{\link[spida2]{\%less\%}} synonym for \code{\link{setdiff}} as well as
#'       \code{\link[spida2]{\%and\%}} and \code{\link[spida2]{\%or\%}}
#' \item \code{\link[spida2]{labs}} assign, extract and print labels for various objects
#' \item \code{\link[spida2]{pch}} generate plotting character mnemonically
#' \item \code{\link[spida2]{pfmt}} format p-values
#' \item \code{\link[spida2]{print.cat}}
#' \item \code{\link[spida2]{rnd}} round a vector to keep significant digits in variation in values
#' \item \code{\link[spida2]{run}} evaluate a string as a command with try
#' \item \code{\link[spida2]{grepv}} grep(..., value = TRUE)
#' }
#' @section String manipulation functions:
#' Function designed to work smoothly with pipes in \pkg{magrittr}:
#' \itemize{
#' \item \code{\link[spida2]{sub_}} and \code{\link[spida2]{gsub_}} handle substitution in a pipeline and return a factor if the input is a factor. 
#' \item \code{\link[spida2]{name}} changes the names of an object and returns the renamed object 
#' }
#' @docType package
#' @name spida2
NULL