\name{coffee}
\Rdversion{1.1}
\alias{coffee}
\concept{artifical data}
\concept{causal inference}
\concept{spida 2009}
\concept{measurement error}
\docType{data}
\title{
Artificial data set on Coffee Consumption, Heart Damage Index and Stress
}
\description{
Artifical data on 20 subjects to illustrate issues of causality with observational data
}
\usage{data(coffee)}
\format{
  A data frame with 20 observations on the following 5 variables.
  \describe{
    \item{\code{Occupation}}{a factor with levels \code{Grad_Student} \code{Professor} \code{Student}}
    \item{\code{Coffee}}{Index of coffee consumption}
    \item{\code{Stress}}{Index of stress}
    \item{\code{Heart}}{Heart Index: higher indicates a less healthy heart}
    \item{\code{Stress2}}{Index of stress that is less reliable than 'Stress'}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(coffee)
## maybe str(coffee) ; plot(coffee) ...
}
\keyword{datasets}
