#' Meaningful contrasts for factors
#'
#' Generate coding values for a factor
#' so factor contrasts have specified interpretations
#' 
#' @param E a \eqn{k-1 \times k} matrix of 
#'        linearly independent factor
#'        level contrasts to be estimated
#'        in regression output. The row names
#'        of E are the names of the estimates
#'        in the regression output. The rows of `E` must
#'        also be linearly independent of the vector of 
#'        coefficients for the `intercept` argument.
#'        Usually, the rows of `E` are contrasts 
#'        summing to 0, and the intercept consists
#'        of a row of constants each equal to \eqn{1/k}.
#' @param intercept coefficients for a linear combination of
#'        of factor levels that defines the intercept. Default
#'        the average of factor levels.
#' @examples
#' days <- c('Sat','Sun','Mon','Tue','Wed','Thu','Fri')
#' data <- expand.grid(id = 1:100, day = factor(days)) 
#' data <- within(
#'   data,
#'   {
#'      y <- 10*rnorm(100)[id] + rnorm(7)[day] + rnorm(days)
#'   })
#' require(nlme)
#' fit <- lme(y ~ day, data, random = ~ 1 |id)
#' summary(fit)
#' #
#' # In the output above, the intercept is the expected value 
#' # for the the 'reference level', which is the level omitted in
#' # the list, i.e. 'Fri' because it is first alphabetically.
#' # The other values are comparisons of each day with
#' # Friday.
#' #
#' data$day <- factor(data$day, 
#'             levels = c('Sat','Sun','Mon',
#'                      'Tue','Wed','Thu','Fri'))
#' E <- rbind(
#'     ":Weekend-Weekday" = c(.5,.5 , -.2, -.2, -.2, -.2, -.2),
#'     ":Sun-Sat"   = c(-1, 1 , 0, 0, 0, 0, 0),
#'     ":Tue-Mon"   = c(0, 0, -1, 1, 0, 0, 0),
#'     ":Wed-Tue"   = c(0, 0, 0, -1, 1, 0, 0),
#'     ":Thu-Wed"   = c(0, 0, 0,  0,  -1, 1, 0),
#'     ":Fri-Thu"   = c(0, 0,  0,  0,  0, -1, 1)
#' )
#' colnames(E) <- levels(data$day)
#' require(spida2)
#' contrasts(data$day) <- coding(E)
#' fit <- lme(y ~ day, data, random = ~ 1 | id)
#' summary(fit)
#' pred <- with(data, pred.grid(day = day))
#' pred
#' contrasts(pred$day) # using pred.grid preserves contrasts
#' pred$fit <- predict(fit, newdata = pred, level = 0)
#' contrasts(pred$day)
#' require(lattice)
#' require(latticeExtra)
#' require(latex2exp)
#' xyplot(fit ~ day, pred, type = 'b')
#' 
#' ww <- waldf(fit, pred = pred)
#' contrasts(ww$day)
#' xyplot(coef ~ day, ww, type = 'b')
#' xyplot(coef ~ day, ww, type = 'b',
#'   fit = ww$coef,
#'   ylab = TeX('Estimated y $\\pm$ SE'),
#'   upper = ww$coef + ww$se,
#'   lower = ww$coef - ww$se,
#'   subscripts = TRUE ) +
#'   layer(panel.band(...))
#'   
#' # A pooled analysis would produce different results
#' # for contrasts.
#' # Because days are balanced within subjects (`id`)
#' # The estimates are the same but the standard
#' # errors and p-value are different.
#' 
#' fit.pooled <- lm(y ~ day, data)   
#' summary(fit.pooled)
#' 
#' # The default treatment contrasts correpond to the
#' # following estimates
#' 
#' E <- rbind(
#'   'Sun-Sat' = c(-1,1,0,0,0,0,0),
#'   'Mon-Sat' = c(-1,0,1,0,0,0,0),
#'   'Tue-Sat' = c(-1,0,0,1,0,0,0),
#'   'Wed-Sat' = c(-1,0,0,0,1,0,0),
#'   'Thu-Sat' = c(-1,0,0,0,0,1,0),
#'   'Fri-Sat' = c(-1,0,0,0,0,0,1)
#' )
#' Cmat <- coding(E, intercept = c(1,0,0,0,0,0,0))
#' Cmat
#' # which is the familiar 'treatment' coding matrix
#' 
#' # A quick way to generate the contrast matrix
#' # from the E matrix without using the "spida2" 
#' # package:
#' Cmat <- solve(rbind(1/ncol(E), E))[, -1]
#' contrasts(data$day) <- Cmat
#' @export   
coding <- function(E, intercept = 1/ncol(E)) {
  # factor coding to get between factor level
  # contrasts defined by E, intercept defined
  # by 'intercept'
  solve(rbind(intercept,E))[,-1]
}
