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
#'   html_document:
#'     toc: true
#'     number_sections: false
#'     toc_depth: 3
#'     highlight: tango
#' ---

options(knitr.table.format = "html") 
#' 
#' Print array with kableExtra
#' 
library(knitr)
library(spida2)
# remotes::install_github('gmonette/spida2@wald-lrt')


# install.packages('kableExtra')
# install.packages('car')

library(kableExtra)

# kpr <- function(x,...) UseMethod('kpr', x, ...)
# kpr.

z <- Titanic
names(dimnames(z)[1]) <-c('blodk\\ndjfkdj') 
dimnames(z)[[1]] <- 1:4
dimnames(z)[[1]] <- 1:4
zz <- dimnames(z)
names(zz)[1] <- 'djf   kldj'
dimnames(z) <- zz
z[,,2,1]
dimnames(z)
apply(z, 1:2, as.matrix)

kpr <- function(a,caption = '', ...) {
  # ...: row_spec
  dots <- list(...)
  print(length(args))
  dd <- dim(a)
  ddn <- dimnames(a)
  ret <-
    if(length(dd) > 2) {
      ap <- aperm(a, c(1,3:length(dd),2))
      ncols <- dd[2]
      nelts <- prod(dd)
      nrows <- nelts / ncols
      dim(ap) <- c(nrows,ncols)
      dnap <- list(rep(ddn[[1]], nrows/dd[1]), ddn[[2]])
      names(dnap) <- names(ddn)[1:2]
      dimnames(ap) <- dnap
      ap
    } else a
  args <- c(ddn[-(1:2)], stringsAsFactors = F)
  gnams <- do.call(expand.grid, args)
  gnams[] <- lapply(seq_along(gnams), function(ii) {
    paste0(names(gnams)[ii],': ', gnams[[ii]])
  })
  gnams <- c(rev(gnams), sep = ' / ')
  gnams <- do.call(paste, gnams)
  gnams
  ret <- kable(ret)  
  ret <- kable_styling(ret, full_width=F)
#  aarg <- list(ret, 0 , angle = +45)
#  ret <- row_spec(ret, 0, angle = -45)
#  ret <- do.call('row_spec', aarg)
   for(nn in names(dots)) {
  #   print(nn)
     aargs <- list(ret)
     aargs <- c(aargs, dots[[nn]])
     ret <- do.call(nn, aargs)
   }
  for( i in seq_along(gnams)) {
    ret <- pack_rows(ret, gnams[i], (i - 1) * dd[1] + 1, i * dd[1])
  }
  ret <- kable_styling(ret, fixed_thead = T)
  ret
}

#' results='asis'
kpr(tab(z), row_spec = list(0, angle = +45), add_header_above =list(c(' '= 1,Gender'=2,Freq = 1) ))  %>% 
  column_spec(1, bold = T, border_right = T) 
#  save_kable(file = 'test.html', self_contained = T)
kpr(tab(z))  %>% 
  column_spec(1, bold = T, border_right = T) 
#  save_kable(file = 'test.html', self_contained = T)
#' 
#' 
#+ results='asis'
kable(Titanic)
Titanic %>%
  kable() %>%
  kable_styling(full_width=F) %>%
  add_header_above(c(' '=1,"Demographics" = 2, Status = 2))
#  save_kable(file = 'test.html')


x <- knitr::kable(head(mtcars), "html")
# Put Row 2 to Row 5 into a Group and label it as "Group A"
pack_rows(x, "Group A", 2, 5) %>%
  kable_styling(full_width=F) %>%
  save_kable(file = 'test.html')
group_rows(x, "Group A", 2, 5) 


Titanic[,,1,1] %>% kable

#+ results='asis', echo = F
for( i in 1:10) {
  cat('\n\n# Header',i, '\n\n')
  cat('text', i, 'continuation','\n\n')
}


