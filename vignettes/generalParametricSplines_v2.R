## ----include = FALSE----------------------------------------------------------
library(spida2)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include=FALSE-----------------------------------------------------
# library(carEx)
library(spida2)
library(magrittr)
library(latticeExtra)

sp <- gspline(0)
Dmat <- environment(sp)$Dmat
Emat <- environment(sp)$Emat
basis <- environment(sp)$basis
Xmat <- environment(sp)$Xmat
Cmat <- environment(sp)$Cmat
Pcon <- environment(sp)$Pcon
Xf <- environment(sp)$Xf

## -----------------------------------------------------------------------------
Xf(0:9, knots = c(3,7), degree = 3)

## -----------------------------------------------------------------------------
Cmat(knots = c(3, 7), degree = c(2, 3, 2), smooth = c(1, 2))

## -----------------------------------------------------------------------------
Emat(knots = c(3, 7), degree = c(2, 3, 2), smooth = c(1, 2))

## -----------------------------------------------------------------------------
sp <- gspline(knots = c(3, 7), degree = c(2, 3, 2), smoothness = c(1, 2))
sp(0:9)

## -----------------------------------------------------------------------------
df <- data.frame(x = 0:10)
set.seed(123)
df <- within(df, y <- -2* (x-5) + .1 * (x-5)^3 + rnorm(x))
df <- rbind(df, data.frame(x = seq(0,10,.1), y = NA))
df <- sortdf(df, ~ x)
plot(y~x, df, pch = 16)
fit <- lm(y ~ sp(x), data = df)
summary(fit)
lines(df$x , predict(fit, df))

## -----------------------------------------------------------------------------
sp <- gspline(knots = c(3,7), degree = c(2,3,2), smoothness = c(1,2))
sp(0:9)

## -----------------------------------------------------------------------------
sp(c(2, 3, 7), D = 1)

## -----------------------------------------------------------------------------
sp(c(3, 3, 3), D = 2, limit = c(-1,0,1))

## -----------------------------------------------------------------------------

xpred <- seq(0,10, .05)

A.1 <- cbind(0, sp(xpred, D = 1))
ww.1 <- as.data.frame(wald(fit, A.1))

A.2 <- cbind(0, sp(xpred, D = 2))
ww.2 <- as.data.frame(wald(fit, A.2))

plot(xpred, ww.1$coef, type = 'l')
plot(xpred, ww.2$coef, type = 'l')
library(latticeExtra)
ww.1$x <- xpred
xyplot(coef ~ x, ww.1, type = 'l',
 	   lower = ww.1$L2, upper = ww.1$U2,
 	   ylab = 'first derivative',
 	   subscripts = TRUE) +
 	layer(gpanel.fit(...))
head(ww.1)

## -----------------------------------------------------------------------------
unemp <- read.csv('http://blackwell.math.yorku.ca/data/USUnemployment.csv')
unemp$date <- as.Date(unemp$date)
head(unemp)
library(lattice)
library(latticeExtra)
xyplot(unemployment ~ date, unemp) + layer(panel.abline(v = as.Date('2008-12-15', col = 'blue')))

toyear <- function(x) {
  # number of years from January 1, 2000
  (as.numeric(x) - as.numeric(as.Date('2000-01-01')))/365.25
}
unemp <- within(
  unemp,
  {
    year <- toyear(date)
    month <- as.numeric(format(date, '%m'))
  })
summary(unemp)

## -----------------------------------------------------------------------------
quintiles <- quantile(unemp$year, 1:4/5)
sp2 <- gspline(quintiles, 2, 1) # quadratic spline
sp3 <- gspline(quintiles, 3, 2) # cubic spline    

## -----------------------------------------------------------------------------
fit2 <- lm(unemployment ~ sp2(year), unemp)
summary(fit2)
unemp$fit2 <- predict(fit2)
fit3 <- lm(unemployment ~ sp3(year), unemp)
summary(fit3)
unemp$fit3 <- predict(fit3)

fit2d <- lm(unemployment ~ sp2d(year), unemp)
summary(fit2d)
unemp$fit2d <- predict(fit2d)
fit3d <- lm(unemployment ~ sp3d(year), unemp)
summary(fit3d)
unemp$fit3d <- predict(fit3d)

pp <- xyplot(unemployment ~ date, unemp, type = 'b',
             col = 'gray',
             key = list(
               corner = c(1,1),
               text = list(lab = c('quadratic','cubic','quadratic discontinuous','cubic discontinuous')),
               lines = list(col= c('red','blue','red','blue'),
                            lty = c(3,3,2,1))
               )) +
  layer(panel.lines(x, unemp$fit3, col = 'blue', lty = 3)) +
  layer(panel.lines(x, unemp$fit2, col = 'red', lty = 3)) + 
  layer(panel.lines(x, unemp$fit3d, col = 'blue', lty = 1)) +
  layer(panel.lines(x, unemp$fit2d, col = 'red', lty = 2)) 
pp


## -----------------------------------------------------------------------------
AIC(fit2, fit3, fit2d, fit3d)

## -----------------------------------------------------------------------------
per3 <- gspline(12 * 1:5/5, 3, 2, periodic = TRUE)
fitper3 <-  lm(unemployment ~ sp3d(year) + per3(month), unemp)
summary(fitper3)
unemp$fitper3 <- predict(fitper3)
pp <- xyplot(unemployment ~ date, unemp, type = 'b',
             col = 'gray') +
  layer(panel.lines(x, unemp$fitper3, col = 'blue'))  
pp

## -----------------------------------------------------------------------------
pred <- data.frame(month = seq(0,24,.01), year = 0)
pred$fitper3 <- predict(fitper3, newdata = pred)
xyplot(fitper3 ~ month, pred, type = 'l',
       xlim = c(0,24),
       scales = list( x = list( at =3 * 0:8)),
       ylab = 'seasonal unemployment') +
  layer(panel.abline(v = 12))

## -----------------------------------------------------------------------------
derivs <- expand.grid(month = seq(0, 24, .01), D = 1:3)
Ld <- with(derivs, per3(month, D = D, limit = -1))
Ld <- cbind(0*Ld[,rep(1,12)], Ld)
derivs <- cbind(derivs, walddf(fitper3, Ld))
derivs $ order <- 
  factor(c('first','second','third')[derivs$D])
derivs $ order <- with(derivs, reorder(order, D))
inds <- which(derivs$month %% (12/5) < .0001 & derivs$D == 3)
derivs$coef[inds] <- NA
xyplot(coef ~ month | order , derivs, type = 'l', layout = c(1,3),
       xlim = c(0,24),
       ylab = 'derivative',
       strip.left = T, strip = F) 


Ldi <- per3(seq(0,24,12/10), D = 3)


## -----------------------------------------------------------------------------
summary(fitper3)
Lfx(fitper3)
derivs$year <- 0
ww <- walddf(
  fitper3, 
  Lfx(fitper3,
      list( 0,
0 * M(sp3d(year)),
1 * M(per3(month, D = D))), 
data = derivs))
head(ww)
xyplot(coef ~ month|order, ww) ####################  WORKED !!!!

## -----------------------------------------------------------------------------
sp1d <- gspline(quintiles_with_crash, 1, 0)  
fit_int <- lm(
  unemployment ~ sp3d(year) + per3(month) + year:per3(month), 
  unemp)
summary(fit_int)
car::Anova(fit_int)
wald(fit_int, ':')


## -----------------------------------------------------------------------------
Lfx(fit_int)
quintiles_with_crash  
range(unemp$year)
pred <- expand.grid(
  month = seq(0,24,.1), 
  date = as.Date(c('1995-01-01','2002-01-01','2009-01-01','2016-01-01')),
  D = 1)
pred$year <- toyear(pred$date)
ww <- walddf(
  fit_int,
  Lfx(fit_int,
      list( 0,
            0 * M(sp3d(year)),
            1 * M(per3(month)),
            1 * M(per3(month)) * year 
      ), pred)
)

xyplot(coef ~ month | sub("-01-01","",date), ww,
       type = 'l',
       xlim = c(0,24),
       scales = list(
         x = list(
           at = c(0,4,8,12,16,20,24),
           labels = c('Jan','May','Sep','Jan','May','Sep','Jan'))),
       main = 'Seasonal component of U.S. unemployment')

unemp$fit_int <- predict(fit_int, unemp)
###################################################################### REDO
pp <- xyplot(unemployment ~ date, unemp, type = 'b',
             col = 'gray',
             key = list(
               corner = c(1,1),
               text = list(lab = c('quadratic','cubic','quadratic disc.','cubic disc.')),
               lines = list(col= c('red','blue','red','blue'),
                            lty = c(1,1,3,3))
               )) +
  layer(panel.lines(x, unemp$fit3d, col = 'blue', lty = 2)) +
  layer(panel.lines(x, unemp$fit_int, col = 'red', lty = 2)) 
pp

# diagnostics
plot(fit_int)
unemp[169,]  # as expected
library(nlme)
fit_int_ar <- 
  gls(
    unemployment ~ sp3d(year) + per3(month) + year:per3(month), 
    unemp, correlation = corAR1(form = ~ 1))
summary(fit_int_ar)
anova(fit_int_ar, fit_int)

## -----------------------------------------------------------------------------
car::Anova(fit_int_ar)
wald(fit_int, ':')
fit_int_ar2 <- update(fit_int_ar, correlation = corARMA(form = ~ 1, p = 2, q = 0))
anova( fit_int_ar2, fit_int_ar, fit_int)



## -----------------------------------------------------------------------------
Sin <- function(x) cbind(sin=sin(2*pi*x),cos=cos(2*pi*x))
fit_fourier_int <- gls(
  unemployment ~ sp3d(year) + Sin(month/12) + Sin(2*month/12) + Sin(3*month/12)+
    Sin(4*month/12) +
    year:(Sin(month/12) + Sin(2*month/12)), 
  unemp, correlation = corAR1(form = ~ 1))
car::Anova(fit_fourier_int)
AIC(fit_fourier_int, fit_int_ar)


fit_factor_int <- lm(
  unemployment ~ sp3d(year) + factor(month) + year:factor(month),
  unemp)

car::Anova(fit_factor_int)
anova(fit_factor_int, fit_int)


## -----------------------------------------------------------------------------


