% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/capply.R
\name{napply}
\alias{napply}
\title{napply: lapply with access to the names of the list argument}
\usage{
napply(X, FUN, ..., DO.CLEAN = TRUE)
}
\arguments{
\item{X}{a list to be passed to \code{\link{lapply}}}

\item{FUN}{a function to be passed to  \code{\link{lapply}}}

\item{...}{optional arguments to FUN}

\item{DO.CLEAN}{remove ,name attribute in result}
}
\description{
lapply allows a function to operate on each element of a list but does not
provide a way to access the name of the elements in the list. napply
provides a mechanism for doing this by setting a ".name" attribute
for each element. This attribute can be extracted with the .name function.
}
\details{
\strong{WARNING:} This sets and uses a '.name' attribute and is likely
to conflict with functions that use such an attribute or that
rely on 'is.vector'.
}
\examples{
\dontrun{
# A common application occurs in manipulating lists of data frames
# (such as might have been read in as a list of similarly
# structured Excel sheets, each for a different year) 
# where the name of the data frame needs to be added as a variable
# to each data frame before merging or 'rbind'ing the data frames.

dlist <- list(
    "2019" = data.frame(x=1:3, y = 1:3 ),
    "2020" = data.frame(x=11:13, y = 11:13 )
)
dlist <- napply( dlist,
   function(d) {
      d$year <- as.numeric(.name(d))
      d
   })
dmerged <- do.call(rbind, dlist)    
dmerged
}
}
