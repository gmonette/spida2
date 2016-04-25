\name{Indonesia}
\Rdversion{1.1}
\alias{Indonesia}
\concept{longitudinal}
\concept{observational}
\concept{SPIDA 2009}
\concept{dichtomous response}
\docType{data}
\title{
Respiratory infections among Indonesian children. Longitudinal observational study with no treatment 
}
\description{
275 preschool children in Indonesia were
examined for up to 6 consecutive quarters
for the present of respiratory infection.
Also recorded presence of xerophthalmia,
an eye condition indicating a vitamin A
deficiency.
}
\usage{data(Indonesia)}
\format{
  A data frame with 1194 observations on 275 children each measured on 1 to 6 occasions
  on the following 10 variables.
  \describe{
    \item{\code{id}}{subject id}
    \item{\code{resp}}{indicator for a respiratory infection}
    \item{\code{age}}{age in months at time of measurement}
    \item{\code{baseline_age}}{age at start of program}
    \item{\code{xerophthalmia}}{indicator for the presence of xerophthalmia}
    \item{\code{sex}}{indicator for sex}
    \item{\code{height_for_age}}{height relative to age}
    \item{\code{stunted}}{indicator}
    \item{\code{time}}{time of measurement in quarters (3 month intervals)}
    \item{\code{season}}{1=Spring, 2=Summer, 3=Autumn, 4=Winter}
  }
}
\details{
  Desc:  275 preschool children in Indonesia were
            examined for up to 6 consecutive quarters
            for the present of respiratory infection.
            Also recorded presence of xerophthalmia,
            an eye condition indicating a vitamin A
            deficiency.


}
\source{
 Subset of cohort from Summer et al. (1984)
            reported in Diggle et al. (2002)
            'Analysis of Longitudinal Data', 2nd ed.
            See: Zeger and Karim, JASA (1991)
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(Indonesia)
## maybe str(Indonesia) ; plot(Indonesia) ...
}
\keyword{datasets}
