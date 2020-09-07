#' Math achievement and ses in a sample of 160 U.S. schools from the 1982 study ``High School and Beyond''
#'
#' This is a classical data set from the field of education used to illustrate
#' multilevel data and models. It is used in the first edition of Bryk and Raudenbush.
#' \code{hsfull} is the complete data set with 160 high schools,
#' \code{hs} is a random subset of 40 high schools, \code{hs1} is a random subset of
#' 80 schools and \code{h2} contains the complement of \code{hs1}.  These two subsets
#' can be used to illustrate split sample validation: develop a model on one half of the data
#' and assess its performance on the other.
#'
#' Each row consists of the data for one student.   \code{hsfull} is the complete data set. \code{hs1} and \code{hs2}
#' are complementary split halves of the schools in the data. \code{hs} is a selection of 40 schools
#' which seems to be a good number of clusters for presentations in class.
#'
#' @format  A data frame with 7185 observations on the following 9 variables.
#' \describe{
#'   \item{\code{school}}{school id}
#'  \item{\code{mathach}}{measure of math achievment}
#'    \item{\code{ses}}{socio-economic status of family}
#'    \item{\code{Sex}}{a factor with levels \code{Female} \code{Male}}
#'    \item{\code{Minority}}{a factor with levels \code{No} \code{Yes}}
#'    \item{\code{Size}}{the size of the school}
#'    \item{\code{Sector}}{a factor with levels \code{Catholic} \code{Public}}
#'    \item{\code{PRACAD}}{a measure of the priority given by the school to academic subjects}
#'    \item{\code{DISCLIM}}{a measure of the disciplinary climate in the school}
#'  }
#' @source Bryk and Raudenbush (COMPLETE)
#' @references Raudenbush, Stephen and Bryk, Anthony (2002),
#'             Hierarchical Linear Models: Applications and Data Analysis Methods, Sage (chapter 4).
#' @concept hierarchical data
#' @concept observational data
#' @keywords datasets
#' @examples
#' \dontrun{
#'   xqplot(hsfull)
#'   xqplot( up( hsfull, ~ school) )
#' }
"hsfull"
#' Longitudinal study of IQ after traumatic brain injuries
#'
#' A subset of data gathered on
#' verbal and performance IQ of patients recovering from coma after traumatic brain injuries
#' \code{iq} is a subset of the original data, \code{iqsim} is simulated set of 10 hypothetical
#' observations on one subject.
#'
#' @format A longitudinal with 331 observations 200 subjects measured on a varying number
#'        of occasions (ranging from 1 to 5) on the following 7 variables.
#'  \describe{
#'    \item{\code{DAYSPC}}{time of measurement in days post recovery from coma}
#'    \item{\code{DCOMA}}{duration of coma rounded to nearest day}
#'    \item{\code{SEX}}{a factor with levels \code{Female} \code{Male}}
#'    \item{\code{AgeIn}}{age at time of injury}
#'    \item{\code{ID}}{identifying subjects}
#'    \item{\code{PIQ}}{performance (or mathematical) IQ}
#'    \item{\code{VIQ}}{verbal IQ}
#'  }
#' @references Wong, Monette, Wiener COMPLETE
#' @keywords datasets
"iq"
#' Monthly unemployment data in the U.S. from 1995-01-01 to 2019-02-01
#' 
#' @format A data set with 290 rows an two variables:
#'   \describe{
#'     \item{\code{date}}{from 1995-01-01 to 2019-02-01}
#'     \item{\code{unemployment}}{rate in the U.S.}
#'   }
#' @keywords datasets
"Unemp"
#' Copy Monthly unemployment data in the U.S. from 1995-01-01 to 2019-02-01
#' 
#' @format A data set with 290 rows an two variables:
#'   \describe{
#'     \item{\code{date}}{from 1995-01-01 to 2019-02-01}
#'     \item{\code{unemployment}}{rate in the U.S.}
#'   }
#' @keywords datasets
"Unemp2"
#' Florida state death penalty data
#' 
#' Dataset detailing death penalty 674 homicide trials in the
#' state of Florida from 1976-1987 with respect to verdict, and
#' victim and defendant race. The data were previously used
#' (Agresti 2012) to demonstrate Simpson's Paradox. This
#' dataset was obtained from the \code{\link{asbio}} package.
#' 
#' A reversal of associations or comparisons may occur as a
#' result of lurking variables or aggregating groups. This is
#' called Simpson's Paradox.
#' 
#' @format A data frame with 8 observations on the following 4
#'   variables.
#'   \describe{
#'     \item{\code{count}}{Counts from cross classification}
#'     \item{\code{verdict}}{Death penalty verdict: No Yes}
#'     \item{\code{d.race}}{Defendant's race: Black, White}
#'     \item{\code{v.race}}{Victim's race: Black, White}
#'   }
#' @references
#' Agresti, A. (2012) Categorical Data Analysis, 3rd edition. New York. Wiley.
#' Radelet, M. L., and G. L. Pierce (1991) Choosing those who will die: race and the death penalty in Florida. Florida Law Review 43(1):1-34.
#' Simpson, E. H. (1951) The Interpretation of interaction in contingency tables. Journal of the Royal Statistical Society Ser. B 13: 238-241.
#' @keywords datasets
"death.penalty"
#' Florida state death penalty data
#' 
#' Dataset detailing death penalty 674 homicide trials in the
#' state of Florida from 1976-1987 with respect to verdict, and
#' victim and defendant race. The data were previously used
#' (Agresti 2012) to demonstrate Simpson's Paradox. This
#' dataset was obtained from the \code{\link{asbio}} package.
#' 
#' A reversal of associations or comparisons may occur as a
#' result of lurking variables or aggregating groups. This is
#' called Simpson's Paradox.
#' 
#' @format A data frame with 8 observations on the following 4
#'   variables.
#'   \describe{
#'     \item{\code{Freq}}{Counts from cross classification}
#'     \item{\code{sentence}}{Death penalty verdict: Internment Death}
#'     \item{\code{defendant}}{Defendant's race: Black, White}
#'     \item{\code{victim}}{Victim's race: Black, White}
#'   }
#' @references
#' Agresti, A. (2012) Categorical Data Analysis, 3rd edition. New York. Wiley.
#' Radelet, M. L., and G. L. Pierce (1991) Choosing those who will die: race and the death penalty in Florida. Florida Law Review 43(1):1-34.
#' Simpson, E. H. (1951) The Interpretation of interaction in contingency tables. Journal of the Royal Statistical Society Ser. B 13: 238-241.
#' @keywords datasets
"Florida_sentencing"
