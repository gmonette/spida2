#' Bias-corrected and accelerated confidence intervals
#'
#' This function uses the method proposed by DiCiccio and Efron (1996)
#' to generate confidence intervals that produce more accurate coverage
#' rates when the distribution of bootstrap draws is non-normal.
#' This code is adapted from the \code{BC.CI()} function within the
#' \code{\link[mediation]{mediate}} function in the \code{mediation} package.
#'
#' @param theta a vector that contains draws of a quantity of interest using bootstrap samples.
#' The length of \code{theta} is equal to the number of iterations in the previously-run
#' bootstrap simulation.
#' @param conf.level the level of the desired confidence interval, as a proportion. Defaults to
#' .95 which returns the 95 percent confidence interval.
#' @param which any combination of 'bca','percentile' or 'normal'; or 'all' to get all.
#' Defaults to 'bca'.
#' @details \eqn{BC_a} confidence intervals are typically calculated using influence statistics
#' from jackknife simulations. For our purposes, however, running jackknife simulation in addition
#' to ordinary bootstrapping is too computationally expensive. This function follows the procedure
#' outlined by DiCiccio and Efron (1996, p. 201) to calculate the bias-correction and acceleration
#' parameters using only the draws from ordinary bootstrapping.
#' @return returns a vector of length 2 in which the first element is the lower bound and the
#' second element is the upper bound
#'
#' @author Jonathan Kropko <jkropko@@virginia.edu> and Jeffrey J. Harden <jharden@@nd.edu>, based
#' on the code for the \code{\link[mediation]{mediate}} function in the \code{mediation} package
#' by Dustin Tingley, Teppei Yamamoto, Kentaro Hirose, Luke Keele, and Kosuke Imai.
#' @references DiCiccio, T. J. and B. Efron. (1996). Bootstrap Confidence Intervals. \emph{Statistical Science}.
#' 11(3): 189–212. \url{https://doi.org/10.1214/ss/1032280214}
#' @seealso \code{\link[coxed]{coxed}}, \code{\link[rms]{bootcov}}, \code{\link[mediation]{mediate}}
#' @export
#' @examples
#' theta <- rnorm(1000, mean=3, sd=4)
#' bca(theta, conf.level = .95)
#' bca(theta, conf.level = .95, c('norm', 'perc'))
#' # 
#' # Why confidence intervals types matter:
#' #
#' bca(rnorm(1000), 'all')
#' bca(exp(rnorm(1000)), 'all')
bca <- function(theta, which = 'bca', conf.level = .95){
  
  if(is.numeric(which)) {
    conf.level <- which
    which <- 'bca'
  }

  if(length(which) == 1 && which == 'all') which <- c('bca','normal','percentile')
  if(var(theta)==0){
    lower <- mean(theta)
    upper <- mean(theta)
    return(c(lower, upper))
  }

  if(('bca' %in% which) & (max(theta)==Inf | min(theta)==-Inf)){
    warning("bca() function does not work when some values are infinite")
  }

  low <- (1 - conf.level)/2
  high <- 1 - low
  sims <- length(theta)
  z.inv <- length(theta[theta < mean(theta)])/sims
  z <- qnorm(z.inv)
  U <- (sims - 1) * (mean(theta, na.rm=TRUE) - theta)
  top <- sum(U^3)
  under <- 6 * (sum(U^2))^{3/2}
  a <- top / under
  lower.inv <-  pnorm(z + (z + qnorm(low))/(1 - a * (z + qnorm(low))))
  lower <- quantile(theta, lower.inv, names=FALSE)
  upper.inv <-  pnorm(z + (z + qnorm(high))/(1 - a * (z + qnorm(high))))
  upper <- quantile(theta, upper.inv, names=FALSE)
  if(length(which) == 1 && which == 'bca') return(c(lower,upper))
  
  ret <- data.frame(
    lower = c(lower, mean(theta) - qnorm(high)*sd(theta), quantile(theta,low)),
    upper = c(upper, mean(theta) + qnorm(high)*sd(theta), quantile(theta,high))
  )
  rownames(ret) <- c('bca','normal','percentile')
  attr(ret,'conf.level') <- c(conf.level=conf.level) 
  class(ret) <- c('conf_interval', class(ret))
  ret <- ret[which,]
  if(dim(ret)[1] == 1) ret <- unlist(ret)
  ret
}
#' @export
print.conf_interval <- function(x,...) {
  NextMethod()
  print(attr(x,'conf.level'))
}

