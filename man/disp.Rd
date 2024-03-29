% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/disp.R
\name{disp}
\alias{disp}
\title{Print the name and value of a variable in a function}
\usage{
disp(x, head = deparse(substitute(x)))
}
\arguments{
\item{x}{value to print}

\item{head}{(default: deparse(substitute(x))) heading preceding printed value}
}
\description{
Useful for debugging as a substitute for printing the value of a variable in a function. It also prints the name of the variable.
}
\details{
Prints object if \code{options(verbose=TRUE)}
}
