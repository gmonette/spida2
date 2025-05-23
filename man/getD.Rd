% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wald-lrt.R, R/wald.R
\name{getD}
\alias{getD}
\title{Get model data frame with optional additional rows}
\usage{
getD(model, add = NULL)

getD(model, add = NULL)
}
\arguments{
\item{model}{whose model frame is used to obtain
variable names and, particularly, correct
levels for factors used as predictors in the model.}

\item{add}{a data frame with values of predictor variables}
}
\value{
a data frame concatenating the rows of the
model frame with those provided in 'add' with
an additional variable '.source' with value
'model' or 'add' to indicated the source of the row.

a data frame concatenating the rows of the
model frame with those provided in 'add' with
an additional variable '.source' with value
'model' or 'add' to indicated the source of the row.
}
\description{
Facilitates obtaining the data frame for a model
along with optional additional creating rows of
predictor values to graph fitted values.

Facilitates obtaining the data frame for a model
along with optional additional creating rows of
predictor values to graph fitted values.
}
\details{
Might supersede \code{\link{getData}} and other
processes to create a predictor data frame. Relies on
\code{\link{getModelData}}, which, in strong contrast
with \code{\link{getData}}, does not use methods
for each modelling methods. If some modelling methods
don't work, then 'getD' will be come a generic function
with the current definition as the default method.

Might supersede \code{\link{getData}} and other
processes to create a predictor data frame. Relies on
\code{\link{getModelData}}, which, in strong contrast
with \code{\link{getData}}, does not use methods
for each modelling methods. If some modelling methods
don't work, then 'getD' will be come a generic function
with the current definition as the default method.
}
