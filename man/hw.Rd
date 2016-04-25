\name{hw}
\Rdversion{1.1}
\alias{hw}
\alias{hwoutliers}
\docType{data}
\title{
Articial data with classical outliers.
}
\description{
Articial data on health as predicted by height and weight with outliers to show the effect of
different types of outliers.
}
\usage{
data(hw)

data(hwoutliers)
}

\format{
  A data frame with 4 datasets, indexed by \code{Type}, each with 15 'good' observations. The last 3 datasets eadh
  contain, in addition to the 'good' data, an outlier.
  \describe{
    \item{\code{Weight}}{a measure of body weight}
    \item{\code{Height}}{a measure of height}
    \item{\code{Health}}{an index of health: higher is better}
    \item{\code{Type}}{numerical type, 0, 1, 2 or 3}
    \item{\code{Outlier}}{a factor with levels \code{none} \code{Type 1} \code{Type 2} \code{Type 3}}

  }
}
\details{
\code{hw}  has four datasets identified by \code{Outlier}
indicating the 'type' of outlier depicted in each dataset: none, Type 1: large residual for y with
typical value of x, Type 2: small residual for y with atypical value of x, and Type 3: atypical x
and poorly fitting y.  \code{hwoutliers} contains the rows of 'acceptable' data and one row for each of
the three types of outlier. This format lends itself better to 3-D plots.  Note that this nomenclature
for outlier is not widely adopted!
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(hw)
## maybe str(hw) ; plot(hw) ...
}
\keyword{datasets}
