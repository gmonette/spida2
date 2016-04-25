\name{Drugs}
\Rdversion{1.1}
\alias{Drugs}
\concept{spida 2009}
\concept{longitudinal}
\concept{observational}
\docType{data}
\title{
Longitudinal study of schizophrenia symptoms with non-randomized treatment
}
\description{
This is a subset of a data set analyzed by Y. Korman on schizophrenia symptoms
among patients who were assessed at 6 annual checkups. At each checkup, the
type of drug currently prescribed by the treating physician was recorded.
}
\usage{data(Drugs)}
\format{
  A data frame with 318 observations on 53 subjects each observed at 6 annual
  checkups with the following 11 variables.
  \describe{
    \item{\code{age1}}{age in years at first checkup}
    \item{\code{yrsill1}}{years since initial diagnosis of schizophrenia at time of first checkup}
    \item{\code{pos}}{number of positive symptoms (symptoms that are not typically present in a normal subject}
    \item{\code{neg}}{number of negative symptoms (characteristics that are absent but normally present in normal subject}
    \item{\code{gen}}{other symptoms}
    \item{\code{total}}{sum of previous three}
    \item{\code{year}}{time of measurement (1 to 6)}
    \item{\code{drug}}{a factor with levels \code{Atypical} \code{Clozapine} \code{Typical}}
    \item{\code{status}}{a factor with levels \code{Married} \code{Separated/Divorced} \code{Single}}
    \item{\code{Sex}}{a factor with levels \code{F} \code{M}}
    \item{\code{Subj}}{a factor identifying each subject}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(Drugs)
## maybe str(Drugs) ; plot(Drugs) ...
}
\keyword{datasets}