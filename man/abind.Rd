% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tab.R
\name{abind}
\alias{abind}
\title{Bind comformable arrays}
\usage{
abind(arr1, arr2, d, facename = "")
}
\arguments{
\item{arr1}{first array}

\item{arr2}{second array}

\item{d}{dimension along which arr1 and arr2 are joined, or, if the
\code{length(dim(arr1)) == length(dim(arr2)) + 1} the dimension in arr1 that
is extended by arr2.}

\item{facename}{Name for the new array dimension}
}
\value{
The returned value, ret, is an array with dimension dim(arr1) except
for the dth dimension where dim(ret)\link{d} == dim(arr1)\link{d} + dim(arr2)\link{d}.
}
\description{
\code{abind} binds two conformable arrays along a dimension.
}
\details{
Warning: abind::abind should be preferred to spida2::abind.

dim( arr1 ) and dim( arr2 ) must be equal except in the dth dimension. If
the length of dim( arr2 ) is 1 less than that of dim( arr1 ), then 'arr2' is
treated as if it had dropped the dth dimension with size 1.

If length(dim(arr2)) == length(dim(arr1)) - 1, then arr2 is treated as if it
dropped a dimension of size 1 in the dth dimension. 'facename' is used as
the name of the dimnames list for this dimension.
}
\seealso{
\code{\link[base]{aperm}}, to permute arrays
}
\author{
Georges Monette
}
