% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{anova.list}
\alias{anova.list}
\alias{AIC.list}
\alias{BIC.list}
\title{anova, AIC and BIC on model lists}
\usage{
\method{anova}{list}(x)

\method{AIC}{list}(x)

\method{BIC}{list}(x)
}
\arguments{
\item{object}{an optionally named list of models of a type
for which \code{anova(mod1, mod1, ...)} would produce
meaningful output.}
}
\description{
These methods work on lists of models in a manner similar
to that in which \code{do.call(anova, modelList)} should work -- but does not.
}
\section{Functions}{
\itemize{
\item \code{AIC(list)}: AIC method

\item \code{BIC(list)}: BIC method

}}
