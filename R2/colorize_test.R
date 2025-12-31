#' ---
#' title: "default.R"
#' # This is ~/.config/rstudio/templates/default.R
#' # pdf output  ~/.config/rstudio/templates/default.pdf
#' author: ""
#' date: "2025-03-24"
#' # minimalist xournalpp screen: 13.18 in x 7.375 in
#' geometry: top=.4in,left=.4in,right=.4in,bottom=.45in,paperwidth=6.59in,paperheight=3.682in
#' output: bookdown::pdf_document2
#' header-includes:
#'   - \usepackage{pdfpages}
#'   - \pagenumbering{arabic}
#' ---
#' \pagenumbering{gobble}
#' \pagenumbering{arabic}
#' \raggedright
#' \tableofcontents
#'
#' # Reading a Stata file
#'
#' rio doesn't do a good job here.
#'
library(haven)

dd <- read_dta('~/Courses_taken/2025_Discrimination_Research/session1.dta')


# Conversion of Stat (.dta) objects as read by 'haven::read_sta()'
conv <- function(object, ...) UseMethod("conv")
conv.haven_labelled <- function(object,...) haven::as_factor(object)
conv.default <- function(object,...) object
conv.data.frame <- function(object,...){
  object[] <- lapply(object, conv)
  as.data.frame(object)
}
dd <- conv(dd)
head(dd)
lapply(dd, class)
#'
#' # Raster pdf for fast plotting of very large vector graphics
#'
# install.packages('rasterpdf')
library(spida2)
library(rasterpdf)
# help(p=rasterpdf)
ls(2)
#'
system.time({
  file <- paste0("/tmp/tmp",paste0(sample(letters,16,T), collapse = ''), '.pdf')
  raster_pdf(file, width = 12, height = 6, res = 600)
  print(xqplot(dd, mfrow = c(2,3)))
  dev.off()
}
)
#' \includepdf[pages=-,pagecommand={}]{`r file`}
#'
#' # BUGS
#'
#' - 2025-07-18: R version 4.5.1 (2025-06-13),
#'   - `as.data.frame.matrix` has methods in different packages
#'     some of which produce weird results. Can use:
#'     `base:::as.data.frame.matrix` explicitly.
#'
#' # lme classes
#'
#' random:
#' list( id = pdBlocked(list(pdIdent(~control-1),pdSymm(~patient+years_post-1))
#' list( id = pdBlocked(list(pdIdent(~ 0 + control),pdSymm(~ 0 + patient + years_post))
#' 
#' 
#' weights:
#'  weights=varIdent(form = ~1 | gender)
#'  
#' correlation:
#'    correlation = corCompSymm(form = ~ 1 | id) 
#'    
#' SOMETHING ABOUT HAVING INTEGERS FROM 1 to ...
#'  
#' 
#' # Colorize
#' 
library(colorize)
#'
#' This should be `r colorize('red','red')` and this 
#' should be `r colorize('cyan','cyan')`.
#' 
#' \Huge works in pdf but not html
#' 
#' Is this huge?
#' 
#' This should be `r colorize('red','red')` and this 
#' should be `r colorize('cyan','cyan')`.
#' 
#' \bf works in pdf but not html
#' 
#' Is this bold?
#' 
#' This should be `r colorize('red','red')` and this 
#' should be `r colorize('cyan','cyan')`.
#' 
#' \normalfont
#' 
#' Is this normal?
#' 
#' This should be `r colorize('red','red')` and this 
#' should be `r colorize('cyan','cyan')`.


#'

