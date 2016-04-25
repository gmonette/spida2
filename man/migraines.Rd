\name{migraines}
\Rdversion{1.1}
\alias{migraines}
\concept{longitudinal dichotomous data}
\concept{one group pre post design}
\concept{binary response}
\concept{dichotomous response}
\concept{seasonal periodic effects}
\concept{SPIDA 2009}
\docType{data}
\title{
Effectiveness of migraine treatment
}
\description{
Subset of data on migraine treatments collected by Tammy Kostecki-Dillon. The data consists of headache logs kept by 133 patients in
a treatment program in which bio-feedback was used to attempt to reduce migraine frequency and severity. Patients entered the program at different times
over a period of about 3 years. Patients were encouraged to begin their logs four weeks before the onset of treatment and to continue for one month afterwards,
but only 55 patients have data preceding the onset of treatment.
}
\usage{data(migraines)}
\format{
  A data frame with 4152 observations on the following 9 variables.
  \describe{
    \item{\code{id}}{patient id}
    \item{\code{ha}}{headache indictor: 1 if a headache has occurred on the day recorded, 0 otherwise}
    \item{\code{time}}{time in days relative to the onset of treatment which occurs at time == 0}
    \item{\code{dos}}{day of study in days from January 1 in the first year of the study}
    \item{\code{hatype}}{type of migraine headache experienced by a subject: \code{Aura} \code{Mixed} \code{No Aura}}
    \item{\code{female}}{an indicator variable}
    \item{\code{age}}{age at onset of treatment}
    \item{\code{airq}}{a measure of air quality}
    \item{\code{medication}}{medication status: \code{Reduced} \code{Continuing} \code{None}}
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
data(migraines)
## maybe str(migraines) ; plot(migraines) ...
}
\keyword{datasets}
