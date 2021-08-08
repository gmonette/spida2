#' Overdispersion with glmer
#' 
#' Test for overdispersion and compute overdispersion factor for generalized linear mixed models.
#' 
#' Two functions `overdisp_fun` and `quasi_table` are adapted from the
#' from the 
#' \href{https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#testing-for-overdispersioncomputing-overdispersion-factor}{GLMM FAQ}
#' by Ben Bolker.
#' 
#' @param model a model fitted using a generalized linear mixed model family without
#'              an independent scale parameter, i.e. `poisson` or `binomial` with counts
#'              (not with a 0/1 response for logistic regression)
#'
#' @aliases overdispersion
#' @examples
#' library(lme4)
#' library(spida2)
#' set.seed(123)
#' zd <- expand.grid(plate = 1:100) #treat = c('A','B'), batch = 1:10, plate = 1:5)
#' zd <- within(zd,
#'      {  
#'          treat <- plate %% 2 + 1
#'          batch <- plate %% 10 + 1
#'          eta <- rnorm(10)[batch]+ 0 * rnorm(100) + 5 * (treat == 2)
#'          eta_od <- rnorm(10)[batch]+ rnorm(100) + 5 * (treat == 2)
#'          count <- rpois(nrow(zd), exp(eta))
#'          count_od <- rpois(nrow(zd), exp(eta_od))
#'      }
#' )
#' fit <- glmer(count ~ treat + (1|batch), zd, family = 'poisson')
#' fit_od <- glmer(count_od ~ treat + (1|batch), zd, family = 'poisson')
#' summary(fit)
#' summary(fit_od)
#' overdisp_fun(fit)
#' overdisp_fun(fit_od)
#' quasi_table(fit)
#' quasi_table(fit_od)
#' wald(fit)
#' 
#' wald(fit_od)
#' wald(fit_od, overdispersion = T)
#' @export   
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
#' @describeIn overdisp_fun display coefficient table with standard errors adjusted for overdispersion
#' @export   
quasi_table <- function(model,ctab=coef(summary(model)),
                        phi=overdisp_fun(model)["ratio"]) {
  qctab <- within(as.data.frame(ctab),
                  {   `Std. Error` <- `Std. Error`*sqrt(phi)
                  `z value` <- Estimate/`Std. Error`
                  `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
                  })
  return(qctab)
}