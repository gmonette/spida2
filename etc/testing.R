#' Testing
#'
library(yscs)
library(xtable)
library(spida)
rm(list=ls())
options(xtable.type = 'latex')
#debug(yscs::rpfmt.wald)
fit <- lm(Heart ~ Coffee+ Stress, coffee)
w <- wald(fit)
#+ results='asis'
xtable(fit)
#+ results='asis'
w <- rpfmt(w)
xtable(t(w))

