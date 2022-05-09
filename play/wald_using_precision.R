library(spida2)
zd <- data.frame(x = 1:4)
zd <- within(
  zd,
  {
    x <- 1:4
    z <- c(1,2,1,2)
    y <- x + rnorm(x)
    a <- c('a','b','a','b')
    w <- c(1,1,2,1)
    yc <- c(1,1,1,0)
  }
)


fit <- lm(y ~ x + z + a, zd)
summary(fit)
fitw <- lm(y ~ x + z + a, zd, weights = w)
summary(fitw)
getFix(fitw)

fitg <- glm(yc ~ x + z + a, zd, weights = w, family = binomial )
summary(fitg)
vcov(fitg)
getX(fitg)
weights(fitg)
xpx <- crossprod( getX(fitg) * sqrt(weights(fitg)))
xpx
svd(xpx)
summary.glm
str(fitg)
fitg$weights

g <- getAnywhere


g(vcov.glm)
g(vcov.summary.glm)
g(summary.glm)
g(glm)
chol2inv
c2i <- chol2inv(diag(4)+.2)
#'
#' Working with factors
#'
dd <-
  data.frame(
  fruit = c('orange', 'apple','apple','lemon','orange'),
  weights = c(2, 3, 4, 1, 3),
  cost = c(2, 1, 2, 1, 3)
)
#' The default order 
tab(dd, ~ fruit)
#' To make 'lemon' the reference level
dd <-
  within(
    dd,
    {
      fruit2 <- relevel(factor(fruit), 'lemon')
    }
  )
tab(dd, ~ fruit2)
#'
#' To sort the levels according to the mean values of another variable
#' 
dd <-
  within(
    dd,
    {
      fruit3 <- reorder(factor(fruit), weights)
    }
  )
tab(dd, ~ fruit3)
#'


summary(lm(cost ~ fruit3, dd))