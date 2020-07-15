#' ---
#' title: "Regression in R"
#' author: "Georges Monette"
#' date: "July 2020"
#' fontsize: 12pt
#' header-includes:
#' - \usepackage{amsmath}
#' - \usepackage{geometry}
#' - \usepackage{longtable}
#' - \geometry{papersize={6in,4in},left=.2in,right=.2in,top=.1in,bottom=.1in}
#' - \newcommand{\var}{\mathrm{Var}}
#' - \raggedright
#' output:
#'   pdf_document:
#'     toc: true
#'     number_sections: true
#'     toc_depth: 3
#'     highlight: tango
#' ---


#' 
#' Print array with kableExtra
#' 
library(knitr)
library(spida2)
# remotes::install_github('gmonette/spida2@wald-lrt')


# install.packages('kableExtra')
# install.packages('car')

# library(kableExtra)

#' 
#' 
#+ results='asis'
kable(Titanic)
Titanic

Titanic[,,1,1] %>% kable

