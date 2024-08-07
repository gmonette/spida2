% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tab.R
\name{tab}
\alias{tab}
\alias{tab.formula}
\alias{tab.data.frame}
\alias{tab.default}
\alias{Tab}
\alias{tab.table}
\alias{tab.matrix}
\alias{tab.array}
\title{Table of frequencies or relative frequencies bordered with totals and
including NAs}
\usage{
tab(x, ...)

\method{tab}{table}(x, ...)

\method{tab}{matrix}(x, ...)

\method{tab}{array}(x, ...)

\method{tab}{formula}(fmla, data = sys.frame(sys.parent()), ...)

\method{tab}{data.frame}(
  x,
  fmla,
  total.margins = TRUE,
  useNA = "ifany",
  pct = NULL,
  pr = NULL,
  test = FALSE,
  weights = NULL,
  na.rm = NULL,
  all.label = "All",
  simulate = FALSE,
  B = 2000,
  ...
)

\method{tab}{default}(
  ...,
  total.margins = TRUE,
  pct = NULL,
  pr = NULL,
  useNA = "ifany",
  test = simulate,
  weights = NULL,
  na.rm = NULL,
  all.label = "All",
  simulate = FALSE,
  B = 2000
)
}
\arguments{
\item{\dots}{as with \code{table}, one or more objects which can be interpreted as factors (including character strings), or a list (or data frame) whose components can be so interpreted.}

\item{fmla}{a formula whose right-hand side names the variables to be used for tabulation. The optional left-hand side specifies a variable to be used for weights.}

\item{data}{a data frame in which formula are interpreted}

\item{total.margins}{if FALSE, generate table without margins}

\item{useNA}{whether to include NA levels. The default is "ifany". Can also
be set to "no" or "always".}

\item{pct}{margins to be scaled to sum to 100. This is the vector of margin indices on which percentages are conditioned. For example, with the call \code{tab(~ A + B + C, data, pct = 2:3)}, the table will contain conditional percentages for variable
\code{A} conditional of combinations of \code{B} and \code{C}.}

\item{pr}{margins to be scaled to sum to 1. This is the vector of margin indices on which percentages are conditioned.}

\item{test}{(default FALSE) use \code{\link{chisq.test}}}

\item{weights}{(not working temporarily) instead of generating a frequency table, generate a table
with the sum of the weights}

\item{simulate}{(default FALSE) simulate p-value with  \code{\link{chisq.test}}}

\item{B}{(default 2000) number of replications for simulation}

\item{keep}{names of margins to keep with 'Tab', default = "All". To drop
all margins, use default = "".}
}
\value{
An object of class 'table' of dimension equal to the number of
variables, with optional margins showing totals. Elements of the matrix can
be frequencies, relative frequencies, percentages or sums of weights.
}
\description{
Generates a table of frequencies or relative frequencies or relative
percentages
}
\section{Methods (by class)}{
\itemize{
\item \code{tab(table)}: method class table

\item \code{tab(matrix)}: method for matrices

\item \code{tab(array)}: method for arrays

\item \code{tab(formula)}: method for formulas

\item \code{tab(data.frame)}: method for data frames

\item \code{tab(default)}: default method

}}
\examples{
titanic <- as.data.frame(Titanic)
head(titanic)
tab(titanic, Freq ~ Sex + Survived + Age , test = T)
tab(titanic, Freq ~ Sex + Age)
tab(titanic, Freq ~ Sex + Survived + Age)
round(tab(titanic, Freq ~ Sex + Survived + Age,
    pct = c(1,3)),2)
round(Tab(titanic, Freq ~ Sex + Survived + Age,
    pct = c(1,3)),2)
round(Tab(titanic, Freq ~ Sex + Survived + Age,
    pct = c(1,3), keep = ""),2)
}
\seealso{
\code{\link{tab_}} to drop "Total" margins and
\code{\link{tab__}} to drop "Total" and "All" margins.
}
\author{
Georges Monette
}
