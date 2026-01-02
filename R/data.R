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
#' Florida state death penalty data 1976-1987
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
"Florida_sentences"
#' Berkeley Admissions Data by Gender and Department
#'
#' This is a famous data set used to illustrate Simpson's Paradox
#' using admissions data at the University of California at Berkeley
#' in 1973. The admissions figures showed that men applying
#' to graduate school were more likely
#' to be admitted than women.  Studies showed that qualifications
#' were similar for men and women, so the data were further
#' analyzed to discover
#' which departments should be held responsible for the gender
#' bias in admissions.
#' 
#' Surprisingly, in some departments women had a higher 
#' probability of admission. There were a few departments in which men
#' were more likely to be admitted but the difference was small.
#' 
#' The data and analysis are often held as an example in which a more careful
#' analysis showed that there was no gender bias favouring men.
#' Careful reflection might lead to a different conclusion.  
#'
#' @format A data set with four variables:
#'  \describe{
#'    \item{\code{Dept}}{Department applied to. For anonymity the departments are identified by letters.}
#'    \item{\code{Gender}}{of the applicant.}
#'    \item{\code{Status}}{whether \code{Admitted} or \code{Denied}}
#'    \item{\code{count}}{number of applicants}
#'  }
#' @references https://en.wikipedia.org/wiki/Simpson%27s_paradox
#' @keywords datasets
"Berkeley"
#' Smoking data by country and gender in approximately 2004
#' 
#' @format A data set with 24 variables with self-descriptive names:
#'  \describe{
#' \item{\code{country}}{Country}
#' \item{\code{iso3}}{iso 3-letter country abbreviation}
#' \item{\code{region}}{Region: 3-letter region designation}
#' \item{\code{HealthExpPC.Govt.exch}}{Health Expenditures per capita by government in US dollars using current currency exchange rates}
#' \item{\code{HealthExpPC.Tot.ppp}}{Total Health Expenditures per capita in US dollars using purchasing power parity exchange rates}
#' \item{\code{HealthExpPC.Govt.ppp}}{Health Expenditures per capita by government in US dollars using purchasing power parity exchange rates}
#' \item{\code{HealthExpPC.Tot.exch}}{HealthExpPC.Tot.exch}
#' \item{\code{total}}{same as \code{HealthExpPC.Tot.ppp}}
#' \item{\code{govt}}{same as \code{HealthExpPC.Govt.ppp}}
#' \item{\code{private}}{equal to \code{total - govt}}
#' \item{\code{sex}}{Sex}
#' \item{\code{lifeexp.Birth}}{lifeexp.Birth}
#' \item{\code{lifeexp.At60}}{lifeexp.At60}
#' \item{\code{smoking.tobacco.current}}{smoking.tobacco.current}
#' \item{\code{smoking.tobacco.daily}}{smoking.tobacco.daily}
#' \item{\code{smoking.cig.current}}{smoking.cig.current}
#' \item{\code{smoking.cig.daily}}{smoking.cig.daily}
#' \item{\code{Pop.Total}}{Total population}
#' \item{\code{Pop.MedAge}}{Median age}
#' \item{\code{Pop.pCntUnder15}}{percent of population under 15 years of age}
#' \item{\code{Pop.pCntOver60}}{Pop.pCntOver60}
#' \item{\code{Pop.pCntAnnGrowth}}{Pop.pCntAnnGrowth}
#' \item{\code{consumption.cigPC}}{consumption.cigPC}
#' \item{\code{hiv_prev15_49}}{hiv_prev15_49}
#'  }
#' @examples
#' if(!require("p3d")) remotes::install_github('gmonette/p3d')
#' library(p3d)
#' library(spida2) 
#' data(Smoking3)
#' dd <- Smoking3
#' head(dd)
#' dd$cig <- dd$consumption.cigPC
#' dd$Cigarettes <- dd$consumption.cigPC
#' dd$pop <- dd$Pop.Total
#' dd$Life <- dd$lifeexp.Birth
#' dd$region <- factor(dd$region)
#' dd$log.total <- log(dd$total)
#' dd$Health <- dd$total
#' 
#' ds <- subset( dd, sex == "BTSX")  # Both sexes combined
#' rownames(ds) <- ds$country
#' Init3d(cex = 2)
#' Plot3d(Life ~ Cigarettes + Health | region, ds, phi = 0, theta = 0, fov = 0)
#' Plot3d(lifeexp.At60 ~ govt + private | region, ds, phi = 0, theta = 0, fov = 0)
#' Id3d()
#' Plot3d_par(col = c('#550000','#005500','#000055'))
#' Plot3d(lifeexp.At60 ~ govt + private | region, ds, phi = 0, theta = 0, fov = 0)
#' 
#' @keywords datasets
"Smoking3"
#' Smoking data by country and gender in approximately 2004
#' 
#' @format A data set with 24 variables with self-descriptive names:
#'  \describe{
#' \item{\code{country}}{Country}
#' \item{\code{iso3}}{iso 3-letter country abbreviation}
#' \item{\code{region}}{Region: 3-letter region designation}
#' \item{\code{HealthExpPC.Govt.exch}}{Health Expenditures per capita by government in US dollars using current currency exchange rates}
#' \item{\code{HealthExpPC.Tot.ppp}}{Total Health Expenditures per capita in US dollars using purchasing power parity exchange rates}
#' \item{\code{HealthExpPC.Govt.ppp}}{Health Expenditures per capita by government in US dollars using purchasing power parity exchange rates}
#' \item{\code{HealthExpPC.Tot.exch}}{HealthExpPC.Tot.exch}
#' \item{\code{total}}{same as \code{HealthExpPC.Tot.ppp}}
#' \item{\code{govt}}{same as \code{HealthExpPC.Govt.ppp}}
#' \item{\code{private}}{equal to \code{total - govt}}
#' \item{\code{sex}}{Sex}
#' \item{\code{lifeexp.Birth}}{lifeexp.Birth}
#' \item{\code{lifeexp.At60}}{lifeexp.At60}
#' \item{\code{smoking.tobacco.current}}{smoking.tobacco.current}
#' \item{\code{smoking.tobacco.daily}}{smoking.tobacco.daily}
#' \item{\code{smoking.cig.current}}{smoking.cig.current}
#' \item{\code{smoking.cig.daily}}{smoking.cig.daily}
#' \item{\code{Pop.Total}}{Total population}
#' \item{\code{Pop.MedAge}}{Median age}
#' \item{\code{Pop.pCntUnder15}}{percent of population under 15 years of age}
#' \item{\code{Pop.pCntOver60}}{Pop.pCntOver60}
#' \item{\code{Pop.pCntAnnGrowth}}{Pop.pCntAnnGrowth}
#' \item{\code{consumption.cigPC}}{consumption.cigPC}
#' \item{\code{hiv_prev15_49}}{hiv_prev15_49}
#'  }
#' @examples
#' if(!require("p3d")) remotes::install_github('gmonette/p3d')
#' library(p3d)
#' library(spida2) 
#' data(Smoking3)
#' dd <- Smoking3
#' head(dd)
#' dd$cig <- dd$consumption.cigPC
#' dd$Cigarettes <- dd$consumption.cigPC
#' dd$pop <- dd$Pop.Total
#' dd$Life <- dd$lifeexp.Birth
#' dd$region <- factor(dd$region)
#' dd$log.total <- log(dd$total)
#' dd$Health <- dd$total
#' 
#' ds <- subset( dd, sex == "BTSX")  # Both sexes combined
#' rownames(ds) <- ds$country
#' Init3d(cex = 2)
#' Plot3d(Life ~ Cigarettes + Health | region, ds, phi = 0, theta = 0, fov = 0)
#' Plot3d(lifeexp.At60 ~ govt + private | region, ds, phi = 0, theta = 0, fov = 0)
#' Id3d()
#' Plot3d_par(col = c('#550000','#005500','#000055'))
#' Plot3d(lifeexp.At60 ~ govt + private | region, ds, phi = 0, theta = 0, fov = 0)
#' 
#' @keywords datasets
"Smoking3"
#' Smoking data by country and gender in approximately 2004
#' combining data from various sources.
#' 
#' @format A data set with 34 variables with mainly self-descriptive names:
#'  \describe{
#' \item{\code{Country}}{Country}
#' \item{\code{Continent}}{Continent}
#' \item{\code{LE}}{Overall life expectancy at birth}
#' \item{\code{CigCon}}{Cigarette consumption per capita per annum}
#' \item{\code{LE.q}}{Life expectancy quartile}
#' \item{\code{Cont}}{Continent (6 abbreviated)}
#' \item{\code{Cont2}}{Continent (6 less abbreviated)}
#' \item{\code{HealthExpPC}}{Health Expenditures per capita per annum}
#' \item{\code{Year}}{year}
#' \item{\code{HE}}{Health Expenditure quartile}
#' \item{\code{country}}{Country}
#' \item{\code{iso3}}{iso 3-letter country abbreviation}
#' \item{\code{region}}{Region: 3-letter region designation}
#' \item{\code{HealthExpPC.Govt.exch}}{Health Expenditures per capita by government in US dollars using current currency exchange rates}
#' \item{\code{HealthExpPC.Tot.ppp}}{Total Health Expenditures per capita in US dollars using purchasing power parity exchange rates}
#' \item{\code{HealthExpPC.Govt.ppp}}{Health Expenditures per capita by government in US dollars using purchasing power parity exchange rates}
#' \item{\code{HealthExpPC.Tot.exch}}{HealthExpPC.Tot.exch}
#' \item{\code{total}}{same as \code{HealthExpPC.Tot.ppp}}
#' \item{\code{govt}}{same as \code{HealthExpPC.Govt.ppp}}
#' \item{\code{private}}{equal to \code{total - govt}}
#' \item{\code{sex}}{Sex}
#' \item{\code{lifeexp.Birth}}{lifeexp.Birth}
#' \item{\code{lifeexp.At60}}{lifeexp.At60}
#' \item{\code{smoking.tobacco.current}}{smoking.tobacco.current}
#' \item{\code{smoking.tobacco.daily}}{smoking.tobacco.daily}
#' \item{\code{smoking.cig.current}}{smoking.cig.current}
#' \item{\code{smoking.cig.daily}}{smoking.cig.daily}
#' \item{\code{Pop.Total}}{Total population}
#' \item{\code{Pop.MedAge}}{Median age}
#' \item{\code{Pop.pCntUnder15}}{percent of population under 15 years of age}
#' \item{\code{Pop.pCntOver60}}{Pop.pCntOver60}
#' \item{\code{Pop.pCntAnnGrowth}}{Pop.pCntAnnGrowth}
#' \item{\code{consumption.cigPC}}{consumption.cigPC}
#' \item{\code{hiv_prev15_49}}{hiv_prev15_49}
#'  }
#' @examples
#' if(!require("p3d")) remotes::install_github('gmonette/p3d')
#' library(p3d)
#' library(spida2) 
#' data(Smoking4)
#' dd <- Smoking4
#' head(dd)
#' dd$cig <- dd$consumption.cigPC
#' dd$Cigarettes <- dd$consumption.cigPC
#' dd$pop <- dd$Pop.Total
#' dd$Life <- dd$lifeexp.Birth
#' dd$region <- factor(dd$region)
#' dd$log.total <- log(dd$total)
#' dd$Health <- dd$total
#' ds <- subset( dd, sex == "BTSX")  # Both sexes combined
#' rownames(ds) <- ds$country
#' Init3d(cex = 2)
#' Plot3d(Life ~ Cigarettes + Health | region, ds, phi = 0, theta = 0, fov = 0)
#' Plot3d(lifeexp.At60 ~ govt + private | region, ds, phi = 0, theta = 0, fov = 0)
#' Id3d()
#' Plot3d(lifeexp.At60 ~ govt + private | region, ds, phi = 0, theta = 0, fov = 0)
#' Plot3d(LE ~ govt + private | region, ds, phi = 0, theta = 0, fov = 0)
#' fit <- lm(LE ~ govt + private, ds)
#' Fit3d(fit, resid = TRUE)
#' summary(fit)
#' @keywords datasets
"Smoking4"
