#' 
#' 
#' make a grid from a fit based on getD
#' 
library(spida2)
fit <- lm(mathach ~ ses * Sex, hs)
getD(fit) %>% head
# get makegrid from UTFA

#' kable_array(Titanic, row_spec = list(0, angle = -30), 
#'   add_header_above =list(c(' '= 1,'Gender'=2,Freq = 1)),
#'   column_spec(1, bold = T, border_right = T)) %>%
#'   save_kable(file = 'test.html', self_contained = T)
#'   
#'   
#'   #' ---
#' title: "Inital Analyses of Salary Data"
#' subtitle: "Breakdown of annual increases by year and employee group"
#' author: "Georges Monette"
#' date: "July 2020"
#' fontsize: 12pt
#' header-includes:
#' - \usepackage{amsmath}
#' - \usepackage{geometry}
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
#+ setup, include=FALSE
knitr::opts_chunk$set(echo=F,comment='   ')
knitr::opts_chunk$set(comment=' ')
library(ggplot2)  # to prevent it from hiding function in lattice[Extra]
library(readxl)
library(spida2)
library(lattice)
library(latticeExtra)
# set directory 
if(rstudioapi::isAvailable()) setwd(here())
getwd()
options(max.print = 10000)
# 
# functions
#
fmt <- function(x,dig=0,keep.nas = FALSE,...) {
  xx <- format(round(x,dig), big.mark = ',',...)
  if(!keep.nas) {
    xx[] <- sub('(^ *)NaN( *)$','\\1   \\2', xx)
    xx[] <- sub('(^ *)NA( *)$','\\1  \\2', xx)
  }
  invisible(print(xx, quote = FALSE, right = T))
}

Rbind <-
  function (..., verbose = FALSE) {
    nam <- ".unique."
    x <- list(...)
    x <- lapply(x, as.data.frame)
    nams <- unique(unlist(sapply(x, names)))
    # print(nams)
    while (nam %in% nams) nam <- paste0(nam, ".")
    if (verbose) 
      print(nam)
    x <- lapply(seq_along(x), function(ii) {
      x[[ii]][[nam]] <- rep(ii, length.out = nrow(x[[ii]]))
      x[[ii]]
    })
    fun <- function(x, y) merge(x, y, all = TRUE, sort = FALSE)
    ret <- Reduce(fun, x)
    ret[[nam]] <- NULL
    ret
  }
'add<-' <- function(x,value) Rbind(x,value)


valid <- function(x) !is.na(x)
rd <- function(...) as.data.frame(read_excel(..., col_types = 'text'))
tonum <- function(x) as.numeric(as.character(x))
num <- function(x) {
  if(is.numeric(x)) return(x)
  # to numeric from formats in data
  library(spida2)
  as.character(x) %>%
    gsub('[$,]','',.) %>%
    as.numeric
}
df2num <- function(df, nams) {
  for(n in nams) df[[n]] <- num(df[[n]])
  df
}
if(FALSE) { # test
  num('$1,345.00')
  num('000000323')
}
#
# functions
#
exgrid <- function(...) expand.grid(..., stringsAsFactors = FALSE)
Levels <- function(x) {
  switch(class(x), 
         'factor' = levels(x), 
         numeric = sort(unique(x)), 
         logical = sort(unique(x)),
         levels(factor(x)))
}


#' Cartesian product of predi
makegrid <- function(fit, ...) {
  dots <- list(...)
  dd <- getD(fit)[,-1] # remove response
  dd <- as.list(dd)
  dd <- lapply(dd, unique)
  if(length(dots)) {
    for(nn in names(dots)){
      dd[nn] <- dots[nn] 
    }
  }
  dd <- c(dd, stringsAsFactors = FALSE)
  dd <- do.call(expand.grid, dd)
  dd <- getD(fit, add = dd )
  dd <- subset(dd, .source == 'add')
  dd$.source <- NULL
  dd
}
getD
pr <- function(x,...) UseMethod('pr')
pr.data.frame <- function(x, ...) {
  x %>% head %>% print
  x %>% lapply(Levels) %>% print
  x %>% sapply(class) %>% print
  invisible(x)
}
pr.default <- function(x,...) print(x,...)
:q
