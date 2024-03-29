% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pdInd.R
\name{pdInd}
\alias{pdInd}
\title{Construct pdInd object}
\usage{
pdInd(
  value = numeric(0),
  form = NULL,
  nam = NULL,
  data = sys.parent(),
  cov = NULL
)
}
\arguments{
\item{value}{an optional initialization value}

\item{form}{an optional one-sided linear formula specifying the row/column
names for the matrix represented by \code{object}.}

\item{nam}{and optional vector of character strings specifying the
row/column names for the matrix represented by \code{object}.}

\item{data}{and optional data frame i which to evaluate the variables names
in \code{value} and \code{form}. ...}

\item{cov}{optional position in lower triangle of covariances that are
estimated and, thus, possibly non-zero. The default is that the covariances
in the first column are estimated and possibly non-zero.}

\item{object}{an object inheriting from the class \code{pdInd}, representing
a positive definite matrix with zero covariances except in the first row and
column.}
}
\description{
This function is a constructor for the \code{pdInd} class used to
represent a positive-definite random effects variance matrix
with some specified patterns of zero covariances.
}
\details{
Mixed models in which many predictors have random slopes often fail to converge
in part because of the large number of parameters in the full covariance (G)
matrix for random effects. One way of fitting a more parsimonious model that
includes random slopes is to use \code{\link[nlme]{pdDiag}} with zeros off the
diagonal. However, this also forces zero covariances between random slopes and
and the random intercept, resulting in a model that is not equivariant
with respect to location transformations of the predictors with random
slopes. The alternative remedy of omitting random slopes for some
predictors can lead to biased estimates and incorrect standard errors of
regression coefficients.

The default covariance pattern for \code{pdInd} produces a G matrix with
zero covariances except in the first row and column. If the first random
effect is the intercept, the resulting model assumes independence between random
slopes without imposing minimality of variance over the possibly
arbitrary origin. This imposition is
the reason that having all covariances equal to zero results in a
model that fails to be equivariant under location transformations.

The optional \code{cov} parameter can be used to allow selected non-zero
covariance between random slopes.
and biased eIt is often desirable to fit a parsimonious model with more than one
variable with a random slope.
}
