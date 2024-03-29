% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tab.R
\name{dropLast}
\alias{dropLast}
\title{Drop last facets of array}
\usage{
dropLast(arr, drop = FALSE, keep = NULL)
}
\arguments{
\item{arr}{array}

\item{drop}{drop parameter in subsetting, default FALSE}

\item{keep}{names of facets to be kept, default NULL}
}
\description{
Primarily used to strips totals from a table bordered by totals
by dropping the last facet.
}
\examples{
(arr <- array(1:24,2:4))
arr \%>\% dropLast
arr \%>\% dropLast(drop = TRUE)
tab(iris, ~ Species)
tab(iris, ~ Species) \%>\% dropLast 
tab(iris, ~ Species + I(Sepal.Length > 5)) \%>\% dropLast 
# row percentages:
tab(iris, ~ Species + I(Sepal.Length > 5), pct= 1)
tab(iris, ~ Species + I(Sepal.Length > 5), pct= 1) \%>\% dropLast
tab(iris, ~ Species + I(Sepal.Length > 5), pct= 1) \%>\% dropLast(keep = "All")
}
