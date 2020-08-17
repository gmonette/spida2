#' Extract sample and sampling parameters from Bayesian posterior
#' 
#' @param x a fitted Bayesian regression object, for example,
#'        using \code{rstan}.  
#' @export
getS <- function(x,...) UseMethod('getS')
#' @export
getS.default <- function(x,...) {
  stop('No method implemented for object of class ',class(x),'\n')
}
#' @export
getS.stanfit <- function(sf) {
  # function to create a data frame
  # of samples from a stanfit object
  # together with sampler parameters
  sam <- rstan:::as.data.frame.stanfit(sf)
  pars <- as.data.frame(do.call(rbind,rstan::get_sampler_params(sf, inc_warmup = FALSE)))
  cbind(sam, pars)
}