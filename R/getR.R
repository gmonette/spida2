# 
# New functions

#' Extract R, G or V matrix in a mixed or GLS model
#' 
#' These functions call \code{\link[nlme]{getVarCov}} in the 
#' \pkg{nlme} package. They are intended to have names and functions
#' that are easy to remember.
#' 
#' @param fit model created with \code{\link[nlme]{lme}} or 
#' \code{\link[nlme]{gls}}
#' @param individuals for \code{getR} or \code{getV} with
#'    \code{\link[nlme]{lme}} objects, select the clusters for which
#'    variances are returned. If not specified, the variance of the
#'    first cluster is returned.
#' @return For \code{\link[nlme]{lme}} objects, \code{getG} returns
#'    the between-cluster variance of random effects, \code{getV},
#'    and \code{getR} returns a list with the within-cluster marginal
#'    variance and the within-cluster conditional variance respectively
#'    for the the clusters listed in \code{individuals}. If 
#'    \code{individuals} is missing, the variance of the first
#'    cluster is returned. ISSUE: For 
#'    \code{\link[nlme]{gls}} objects all functions return the same thing but
#'    uninformatively if correlation is clustered and if weights
#'    produce differenct variances in the corresponding positions in 
#'    different clusters. 
#' @examples
#' library(spida2)
#' library(nlme)
#' library(gnew)
#' data <- expand.grid( Xdev = c(-3,-2,-1,0,1,2,3), id = 1:5 )
#'  
#' set.seed(12345)
#' data <- within(data, {
#'   Xmean <- 2*id
#'   X <- Xdev + Xmean
#'   Y <- (-1 + .1*rnorm(max(id)))[id] * Xdev + 
#'     2 * Xmean + .3 * id * rnorm(length(id))
#' })
#' 
#' library(lattice)
#' gd()
#' xyplot(Y ~ X, data, groups = id)
#' fit0 <- lme(Y ~ X, data,
#'             random = ~ 1+ X |id)
#' fit <- lme(Y ~ X, data, 
#'            random = ~ 1 + X | id,
#'            weights = varConstPower(form = ~ fitted(.)),
#'            correlation = corAR1(form = ~ 1 | id),
#'            control = list(returnObject = TRUE))
#' fitgls <- gls(Y ~ X, data, weights = varConstPower(form = ~ fitted(.)),
#'            correlation = corAR1(form = ~ 1|id),
#'            control = list(returnObject = TRUE, maxIter= 1000, 
#'            verbose = TRUE, msMaxIter = 1000,
#'                msVerbose=TRUE))
#' summary(fit)
#' getVarCov(fit)
#' getVarCov(fit, individuals = '2')
#' getVarCov(fit, individuals = '2', type = 'conditional') %>% 
#'   .[[1]] %>% 
#'   diag
#' getVarCov(fit,  type = 'conditional')%>% 
#'   .[[1]] %>% 
#'   diag
#' 
#' getG(fit)
#' getR(fit)[[1]]
#' getV(fit)[[1]]
#' 
#' 
#' (Z <- cbind(1, 2+seq(-3,3)))
#' Z
#' (getG(fit))
#' 
#' Z %*% getG(fit) %*% t(Z)
#' 
#' getV(fit)[[1]]
#' getR(fit)[[1]]
#' sigma(fit)
#' getVarCov(fit, type = 'random.effects')
#' getVarCov(fit)
#' Z %*% getG(fit) %*% t(Z)
#' getV(fit)[[1]] - Z %*% getG(fit) %*% t(Z) - getR(fit)[[1]]
#' 
#' getG(fit0)
#' Z %*% getG(fit0) %*% t(Z)
#' Z %*% getG(fit0) %*% t(Z) %>% svd %>% .$d
#' getR(fit0)
#' sigma(fit0)
#' getV(fit0)
#' Z %*% getG(fit0) %*% t(Z) + getR(fit0)[[1]]
#' @export
getG <- function(fit,...) UseMethod('getG')
#' @describeIn getG default method
#' @export
getG.default <- function(fit,...) "Unknown class"
#' @describeIn getG lme method
#' @export
getG.lme <- function(fit,...) getVarCov(fit,...)
#' @describeIn getG gls method
#' @export
getG.gls <- function(fit,...) getVarCov(fit,...)
#' @rdname getG
#' @export
getR <- function(fit,...) UseMethod('getR')
#' @describeIn getG default method
#' @export
getR.default <- function(fit,...) "Unknown class"
#' @describeIn getG lme method
#' @export
getR.lme <- function(fit, ...) {
  getVarCov(fit, type = 'conditional', ...)
}
#' @describeIn getG gls method
#' @export
getR.gls <- function(fit, ...) {
  getVarCov(fit, type = 'conditional', ...)
}
#' @rdname getG
#' @export
getV <- function(fit,...) UseMethod('getV')
#' @describeIn getG gls method
#' @export
getV.default <- function(fit, ...) "Unknown class"
#' @describeIn getG lme method
#' @export
getV.lme <- function(fit, ...) {
  getVarCov(fit, type = 'marginal', ...)
}
#' @describeIn getG gls method
#' @export
getV.gls <- function(fit, ...) {
  getVarCov(fit, type = 'marginal', ...)
}

if(FALSE) {
library(spida2)
library(nlme)
  library(gnew)
data <- expand.grid( Xdev = c(-3,-2,-1,0,1,2,3), id = 1:5 )

set.seed(12345)
data <- within(data, {
  Xmean <- 2*id
  X <- Xdev + Xmean
  Y <- (-1 + .1*rnorm(max(id)))[id] * Xdev + 
    2 * Xmean + .3 * id * rnorm(length(id))
})

library(lattice)
gd()
xyplot(Y ~ X, data, groups = id)
fit0 <- lme(Y ~ X, data,
            random = ~ 1+ X |id)
fit <- lme(Y ~ X, data, 
           random = ~ 1 + X | id,
           weights = varConstPower(form = ~ fitted(.)),
           correlation = corAR1(form = ~ 1 | id),
           control = list(returnObject = TRUE))
fitgls <- gls(Y ~ X, data, 
           weights = varConstPower(form = ~ fitted(.)),
           correlation = corAR1(form = ~ 1 | id),
           control = list(returnObject = TRUE, maxIter = 1000))

summary(fit)
getVarCov(fit)
getG(fit)
getVarCov(fit, individuals = '2')
getG(fit, ind = '2')
getR(fit) %>% 
  .[[1]] %>% 
  diag
getR(fit)%>% 
  .[[1]] %>% 
  diag

getG(fit)
getR(fit)[[1]]
getV(fit)[[1]]


(Z <- cbind(1, 2+seq(-3,3)))
Z
(getG(fit))

Z %*% getG(fit) %*% t(Z)

getV(fit)[[1]]
getR(fit)[[1]]
sigma(fit)
getG(fit)
getV(fit)[[1]] - Z %*% getG(fit) %*% t(Z) - getR(fit)[[1]]

getG(fit0)
Z %*% getG(fit0) %*% t(Z)
Z %*% getG(fit0) %*% t(Z) %>% svd %>% .$d
getR(fit0)
sigma(fit0)
getV(fit0)

getG(fitgls)

getVarCov(fitgls) 
getG(fitgls)
getR(fitgls)
getV(fitgls, individuals = '5')

# Non-linear
fitnl <- nlme( PIQ ~ b0 + b1*exp(-a*DAYSPC),
      data = iq,
      fixed = list( b0 ~ 1+sqrt( DCOMA ),
                    b1 ~ 1,
                    a ~ 1),
      random = list( ID = b0 ~ 1),
      start = list( fixed = c(
        100, 0,-20.,.05)),
      control = list( maxIter = 100,
                      returnObject = TRUE),
      verbose = T)
summary(fitnl)
# Error message 'getVarCov' not implemented for nlme objects
library(MASS)
fitpql <- fit <- glmmPQL(Y ~ X, data, 
                     random = ~ 1 + X | id,
                     family = gaussian,
                     # weights = varConstPower(form = ~ fitted(.)),
                     # NOTE: weights expects a variable in 'data'
                     # probably because the lme weight mechanism
                     # is used to implement PQL
                     correlation = corAR1(form = ~ 1 | id),
                     control = list(returnObject = TRUE))
getG(fitpql)
getR(fitpql)
#
# Does pdIdent work?



}

