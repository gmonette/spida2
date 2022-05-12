#' ---
#' title: 'wald with rank deficient models'
#' bibliography: bib.bib
#' link-citations: true
#' ---
#' 
library(spida2)

{
  
dd <- expand.grid(a=letters[1:2], b = letters[5:7] , x = 1:10)

set.seed(123)
within(dd, {
  y <- as.numeric(a)*as.numeric(b)*x + .1*rnorm(x)
}) %>% 
  subset(!(a=='b'&b=='g')) -> dd
}


fit <- lm(y ~ a*b*x, dd)
summary(fit)
dd %>% tab(~a+b)
waldx(fit)
ww <- as.data.frame(waldx(fit))
waldx(fit, ww$L)
ww <- as.data.frame(waldx(fit))
ww$L %>% dim
waldx(fit, ww$L)

pred <- up(dd, ~a + b)
pred
wald(fit, pred = pred)
pred <- merge(pred, list(x=5))
ww <- as.data.frame(waldx(fit, pred = pred))
ww$L
waldx(fit, ww$L)

# prediction with cartesian frame
# 
pred <- with(dd, pred.grid(a, b, x =5))



# Workflow for singular models
# 


dd
# differences in each cell from a reference level
# 
pred.grid
pred <- with(dd, pred.grid(a,b, x = 10:12))
pred
waldx(fit, pred = pred)
library(kableExtra)

dd %>% tab(~a+b)
# comparisons with a = a and b = e
# 
pred <- with(dd, pred.grid(a, b, x = 10:12))
ww <- as.data.frame(waldx(fit, pred = pred))


# 
# A small experiment with a simple singular model missing the reference level
# - result: fit just ignore missing level as if it had been dropped!
# 
ds <- expand.grid(a = letters[1:3], x = 1:10) 
ds %>% lapply(class)
set.seed(123)
ds$y <- with(ds, as.numeric(a) * x + rnorm(x))
fit <- lm(y ~ a * x, ds)
fit %>% summary
dss <- subset(ds, a != 'a')             
fits <- lm(y ~ a * x, dss)
car::Anova(fits)
fits %>% summary
