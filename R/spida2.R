#' spida2: Functions used in the Statistical Consulting Service at York University
#' 
#' The spida2 package is a collection of functions and datasets intended primarily
#' for statistical consulting, with a particular emphasis on
#' longitudinal and hierarchical data analysis.
#'
#' @section Wald tests and linear hypothesis matrices:
#' \itemize{
#' \item \code{\link{wald}} quick Wald tests -- optionally building hypothesis matrices from regular expressions
#' \item \code{\link{walddf}} version of wald that returns a data frame
#' \item \code{\link{as.data.frame.wald}} return a data frame from a wald object
#' \item \code{\link{coef.wald}} method to extract estimated coefficients
#' \item \code{\link{print.wald}} printing method
#' \item \code{\link{Lfx}}  creates hypothesis matrices for derivatives and differences. Add example for factor differences
#' \item \code{\link{M}} constructor for M objects to generate portions of design and hypothesis matrices. Used with  \code{\link{Lfx}}
#' \item \code{\link{rpfmt}} format estimated values and p-values from a wald test
#' \item \code{\link{Lall}} for lmer objects
#' \item \code{\link{Lc}} for lmer objects
#' \item \code{\link{Lmu}} for lmer objects
#' }
#' @section Utilities for fitted objects:
#' \itemize{
#' \item \code{\link{getData}} get data frame from a fitted object
#' \item \code{\link{getFix}} get fixed effects from a fitted object
#' \item \code{\link{getX}} get X matrix from fitted object
#' \item \code{\link{gpanel.fit}} plot fitted values and confidence band
#' \item \code{\link{Vcov}} get estimated variance covariance of fixed effects from a fitted object
#' }
#' @section Multilevel data frames:
#' \itemize{
#' \item \code{\link{tolong}} and  \code{\link{towide}} are 
#'   interfaces to \code{\link{stats::reshape}} to facilitate
#'   the typical uses of reshape for longitudinal data.
#' \item \code{\link{capply}}: \code{capply(x,id,FUN)} applies the
#'   function \code{FUN} to chunks of 'x' formed by levels of 'id'.
#'   The result has the same form as 'x' with replication within
#'   chunks, if needed.
#' \item \code{\link{up}} and \code{\link{up_apply}} create summary
#'   data sets consisting, by default, of within-id-invariant variables. 
#'   Summaries of id-varying variables can also be included.
#' }
#' @section Splines -- parametric and non-parametric:
#' \itemize{
#' \item \code{\link{gsp}} creates a function for a generalized 
#'       spline that can
#'       be included in a linear model formula. 
#' \item \code{\link{sc}} creates a spline contrast matrix for 
#'       general spline hypotheses. The matrix can be included
#'       in hypothesis matrices for the \code{\link{wald}} function.
#' \item \code{\link{smsp}} creates a matrix for a smoothing spline.
#' }
#' @section Datasets:
#' \itemize{
#' \item \code{\link{hsfull}} Classical data set on high school math achievement and ses. See Bryk and Raudenbush and many other sources.
#' \item \code{\link{iq}} Recovery after traumatic brain injury
#' \item \code{\link{Drugs}} Longitudinal data on drugs and schizophrenia symptoms. Illustrates role of control variables with non-random assignment. 
#' \item \code{\link{Indonesia}} Xerophthalmia.
#' \item \code{\link{migraines}} Longitudinal data on migraine treatment and weather
#' \item \code{\link{coffee}} Artifical data on coffee, heart damage and stress. Illustrates Simpson's Paradox with continuous predictors.
#' \item \code{\link{hw}} Artifical data on height, weight and health. Illustrates suppression.
#' }
#' @section Miscellaneous utility functions:
#' \itemize{
#' \item \code{\link{gd}} and \code{\link{td}} easy interface to set
#'       graphical parameters for lattice graphsics.
#' \item \code{\link{sortdf}} sort rows of a data frame -- useful in a magrittr pipeline
#' \item \code{\link{assn}} assign -- useful in a magrittr pipeline
#' \item \code{\link{disp}} utility to display value of a variable -- useful for debugging
#' \item \code{\link{getFactorNames}} get names of variables that are factors in a data frame
#' \item \code{\link{\%less\%}} synonym for \code{\link{setdiff}} as well as
#'       \code{\link{\%and\%}} and \code{\link{\%or\%}}
#' \item \code{\link{labs}} assign, extract and print labels for various objects
#' \item \code{\link{pch}} generate plotting character mnemonically
#' \item \code{\link{pfmt}} format p-values
#' \item \code{\link{print.cat}}
#' \item \code{\link{rnd}} round a vector to keep significant digits in variation in values
#' \item \code{\link{run}} evaluate a string as a command with try
#' \item \code{\link{grepv}} grep(..., value = TRUE)
#' }
#' @docType package
#' @name spida2
#' @alias spida
#' @alias yscs
NULL
