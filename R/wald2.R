## wald2.R
##
## This a collection of functions designed to facilitate testing hypotheses
## with Wald tests.
## The methods are readily extended to any fitting method that returns a vector
## of estimated parameters and an estimated variance matrix.
## Extensions are implemented through the 'getFix' generic function.
##
##
## see also Leff
## Changes:
##
## 2020_02_15: adding LRT
## 2018 04 07: commented out isnarow code in getX since it seems in error
## 2014 06 04: changed fit@fixef to fixef(fit) in a number of 'getFix' methods
## October 2, 2011:  modified 'wald' to better handle rank-deficient models
##                   previously columns of L and columns and rows of vcov
##                   corresponding to NAs in beta were deleted ignoring the
##                   possibility that L beta is not estimable if any
##                   non-zero element of L multiplies an NA in beta.
##
## 2013 06 18: added df argument to wald to override denominator dfs. Useful
##             for saturated binomial fits or other binomials where residual
##             dfs are not appropriate.
##
## 2013 09 17: added getFix.multinom with df = Inf
##

# scraps:
# L <- list( 'marginal value of education' =Lform( fit, form = list(0, 1, 2 *
# education, 0, 0, type == 'prof', type == 'wc', 2 * education * (type ==
# 'prof'), 2 * education * (type == 'wc')), data = Prestige)) wald( fit, L )
# chat <- coef( wald( fit, L ), se = 2) xyplot( coef +coefp+coefm ~ education
# | type, cbind(Prestige,chat)[order(Prestige$education),], type = 'l')
# xyplot( chat~ education | type, Prestige)
#
# # This approach can be used to predict responses with a fitting method that
# has a # 'model.matrix' method but does not have a 'predict' method or does
# not return # estimated standard errors with its 'predict' method.
#
# datafit <- model.frame(fit) # or just the data frame used for the fit ww <-
# wald(fit, model.matrix(fit)) datafit <- cbind(datafit, coef(ww, se =2)) #
# ...etc as above

#' General Linear Hypothesis for the 'fixed' portion of a model with Wald test
#'
#' General Linear Hypothesis with Wald test for lm, glm, lme, nlme and
#' lmer objects.  Can be extended to other objects (e.g.) 'glm' by writing
#' 'getFix.glm'
#'
#' Tests a general linear hypothesis for the linear fixed portion of a model.
#' The hypothesis can be specified in a variety of ways such as a hypothesis
#' matrix or a pattern that is used as a regular expression to be matched with
#' the names of coefficients of the model. A number of tools are available to
#' facilitate the generation of hypothesis matrices.
#'
#' Usage:
#'
#' wald(fit, L) where L is a hypothesis matrix
#'
#' wald(fit, 'pat') where 'pat' a regular expression (see ?regex) used to
#' match names of coefficients of fixed effects.  e.g. wald( fit, ':.*:') tests
#' all 2nd and higher order interactions.
#'
#' wald(fit, c(2,5,6)) to test 2nd, 5th and 6th coefficients.
#'
#' wald(fit, list( hyp1= c(2,5,6), H2 = 'pat')) for more than one hypothesis
#' matrix
#'
#' There are number of functions to help construct hypothesis matrices. See in
#' particular \code{\link{Lfx}}.
#'
#' To extend the 'wald' function to a new class of objects, one needs to
#' write a 'getFix' method to extract estimated coefficients, their estimated
#' covariance matrix, and the denominator degrees of freedom for each
#' estimated coefficient. See the examples below for 
#' a \code{\link{getFix}} method for a \code{\link{glm}} object
#' and for 'mipo' objects in the packages 'mice':
#'
#' @param fit a model for which a \code{getFix} method exists.
#' @param Llist a hypothesis matrix or a pattern to be matched or a list of
#'        these
#' @param clevel level for confidence intervals. No confidence intervals if clevel is NULL
#' @param pred prediction data frame to evaluate fitted model using
#'        \code{getX(fit) %*% coef}
#' @param data data frame used as 'data' attribute fot list elements returned only if
#'        the corresponding element of \code{Llist} has a NULL data attribute
#' @param debug (default FALSE) produce verbose information
#' @param maxrows maximum number of rows of hypothesis matrix for which a full
#'        variance-covariance matrix is returned
#' @param full if TRUE, the hypothesis matrix is the model matrix for
#'        \code{fit} such that the estimated coefficients are the predicted values for
#'        the fixed portion of the model. This is designed to allow the calculation of
#'        standard errors for models for which the \code{predict} method does not
#'        provide them.
#' @param pred (default NULL) a data frame to use to create a model matrix. 
#'        This is an alternative to `full` when the model matrix needs to
#'        be based on data frame other than the data frame used for 
#'        fitting the model.
#' @param fixed if \code{Llist} is a character to be used a regular expression,
#'        if \code{fixed} is TRUE \code{Llist} is interpreted literally, i.e.
#'        characters that have a special meaning in regular expressions are
#'        interpreted literally.
#' @param invert if \code{Llist} is a character to be used a regular
#'        expression, \code{invert == TRUE} causes the matches to be inverted so that
#'        coefficients that do not match will be selected.
#' @param method 'svd' (current default) or 'qr' is the method used to find the
#'        full rank version of the hypothesis matrix.  'svd' has correctly identified
#'        the rank of a large hypothesis matrix where 'qr' has failed.
#' @param pars passed to \code{\link[rstan]{extract}} method for stanfit objects.
#' @param include passed to \code{\link[rstan]{extract}} method for stanfit objects.#' 
#' @param help obsolete
#' @return An object of class \code{wald}, with the following components:
#'       COMPLETE
#' @seealso \code{\link{Lform}},
#'          \code{\link{Lfx}}, \code{\link{getX}}, \code{\link{M}},
#'          \code{\link{Lall}},\code{\link{Lc}},\code{\link{Lequal}},
#'          \code{\link{Ldiff}},\code{\link{Lmu}},\code{\link{Lmat}},\code{\link{Lrm}},
#'          \code{\link{as.data.frame.wald}}. To extend to new
#'          models see \code{\link{getFix}}. To generate hypothesis matrices for general
#'          splines see \code{\link{gsp}} and \code{\link{sc}}.
#' @references REFERENCES HERE
#' @examples
#' data(hs)
#' library(nlme)
#' ###
#' ### Using wald to create and plot a data frame with predicted values
#' ###
#'   fit <- lme(mathach ~ (ses+I(ses^2)) * Sex * Sector, hs, random = ~ 1|school)
#'   summary(fit)
#'   pred <- expand.grid( ses = seq(-2,2,.1), Sex = levels(hs$Sex), Sector = levels(hs$Sector))
#'   pred
#'   w <- wald(fit, getX(fit,data=pred)) # attaches data to wald.object so it can be included in data frame
#'   w <- wald(fit, pred = pred)
#'   w <- as.data.frame(w)
#'   head(w)
#'   library(latticeExtra)
#'   xyplot(coef ~ ses | Sector, w, groups = Sex,
#'      auto.key = T, type = 'l',
#'      fit = w$coef,
#'      upper = with(w,coef+2*se),
#'      lower = with(w,coef-2*se),
#'      subscript = T) +
#'      glayer( gpanel.fit(...))
#'
#' wald( fit, 'Sex')  # sig. overall effect of Sex
#' wald( fit, ':Sex') # but no evidence of interaction with ses
#' wald( fit, '\\^2') # nor of curvature
#'
#' # but we continue for the sake of illustration
#'
#' L <- Lform( fit, list( 0, 1, 2*ses, 0, Sex == 'Male', (Sex == 'Male')*2*ses), hs)
#' L
#' (ww <- wald ( fit, L ))
#' wald.dd <- as.data.frame(ww, se = 2)
#' head( wald.dd )
#'
#' require(lattice)
#' xyplot( coef + U2 + L2 ~ ses | Sex, wald.dd,
#'  main= 'Increase in predicted mathach per unit increase in ses')
#' 
#' # Example of a getFix method for a glm oject:
#' 
#' getFix.glm <- function(fit,...) { 
#'   ss <- summary(fit) 
#'   ret <- list() 
#'   ret$fixed <- coef(fit) 
#'   ret$vcov <- vcov(fit) 
#'   ret$df <- rep(ss$df.residual,length(ret$fixed)) 
#'   ret
#' } 
#' 
#' # Example of a getFix method for a mipo object in the mice package:
#' 
#' getFix.mipo <- function( fit, ...){ 
#'   # pooled multiple imputation object in mice 
#'   # 'wald' will use the minimal df for components with non-zero weights 
#'   #  -- this is probably too conservative and should be improved 
#'   ret <- list()
#'   ret$fixed <- fit$qbar 
#'   ret$vcov <- fit$t 
#'   ret$df <- fit$df 
#'   ret 
#' }
#' @export
wald2 <- 
  function(fit, Llist = "", clevel = 0.95,
           pred = NULL,
           data = NULL, debug = FALSE , maxrows = 25,
           full = FALSE, fixed = FALSE,
           invert = FALSE, method = 'svd',
           df = NULL, pars = NULL,...) {
    # New version with support for stanfit
    if (full) return(wald(fit, getX(fit)))
    if(!is.null(pred)) return(wald(fit, getX(fit,pred)))
    dataf <- function(x,...) {
      x <- cbind(x)
      rn <- rownames(x)
      if(length(unique(rn)) < length(rn)) rownames(x) <- NULL
      data.frame(x, ...)
    }
    as.dataf <- function(x, ...) {
      x <- cbind(x)
      rn <- rownames(x)
      if(length(unique(rn)) < length(rn)) rownames(x) <- NULL
      as.data.frame(x, ...)
    }
    unique.rownames <- function(x) {
      ret <- c(tapply(1:length(x), x, function(xx) {
        if(length(xx) == 1) ""
        else 1:length(xx)
      })) [tapply(1:length(x), x)]
      ret <- paste(x, ret, sep="")
      ret
    }
    if(inherits(fit, 'stanfit')) {
      fix <- if(is.null(pars)) getFix(fit) else getFix(fit,pars=pars,...)
      if(!is.matrix(Llist)) stop(
        paste('Sorry: wald needs Llist to be a n x',
              length(fix$fixed),'matrix for this stanfit object'))
    } else {
      fix <- getFix(fit)
    }
    beta <- fix$fixed
    vc <- fix$vcov
    
    dfs <- if(is.null(df) ) fix$df else df + 0*fix$df

    if(is.character(Llist) ) Llist <- structure(list(Llist), names=Llist)
    if(!is.list(Llist)) Llist <- list(Llist)
    
    ret <- list()
    for (ii in 1:length(Llist)) {
      ret[[ii]] <- list()
      Larg <- Llist[[ii]]
      # Create hypothesis matrix: L
      L <- NULL
      if(is.character(Larg)) {
        L <- Lmat(fit,Larg, fixed = fixed, invert = invert)
      } else {
        if(is.numeric(Larg)) {   # indices for coefficients to test
          if(is.null(dim(Larg))) {
            if(debug) disp(dim(Larg))
            if((length(Larg) < length(beta)) && (all(Larg>0)||all(Larg<0)) ) {
              L <- diag(length(beta))[Larg,]
              dimnames(L) <- list( names(beta)[Larg], names(beta))
            } else L <- rbind( Larg )
          }
          else L <- Larg
        }
      }
      if (debug) {
        disp(Larg)
        disp(L)
      }
      # get data attribute, if any, in case it gets dropped
      Ldata <- attr( L , 'data')
      
      ## identify rows of L that are not estimable because they depend on betas that are NA
      Lna <- L[, is.na(beta), drop = FALSE]
      narows <- apply(Lna,1, function(x) sum(abs(x))) > 0
      
      L <- L[, !is.na(beta),drop = FALSE]
      ## restore the data attribute
      attr(L,'data') <- Ldata
      beta <- beta[ !is.na(beta) ]
      
      ## Anova
      if( method == 'qr' ) {
        qqr <- qr(t(na.omit(L)))
        # Qqr <- Q(t(L))
        L.rank <- qqr$rank
        # L.rank <- attr(Qqr,'rank')
        # L.miss <- attr(Qqr,'miss')
        if(debug)disp( t( qr.Q(qqr)))
        L.full <- t(qr.Q(qqr))[ 1:L.rank,,drop=FALSE]
        #L.full <- t(Qqr[!L.miss,])[ 1:L.rank,,drop=F]
      } else if ( method == 'svd' ) {
        if(debug) disp(L)
        #              if(debug)disp( t(na.omit(t(L))))
        #              sv <- svd( t(na.omit(t(L))) , nu = 0 )
        sv <- svd(na.omit(L) , nu = 0 )
        
        if(debug)disp( sv )
        tol.fac <- max( dim(L) ) * max( sv$d )
        if(debug)disp( tol.fac )
        if ( tol.fac > 1e6 ) warning( "Poorly conditioned L matrix, calculated numDF may be incorrect")
        tol <- tol.fac * .Machine$double.eps
        if(debug)disp( tol )
        L.rank <- sum( sv$d > tol )
        if(debug)disp( L.rank )
        if(debug)disp( t(sv$v))
        L.full <- t(sv$v)[seq_len(L.rank),,drop = FALSE]
      } else stop("method not implemented: choose 'svd' or 'qr'")
      
      # from package(corpcor)
      # Note that the definition tol= max(dim(m))*max(D)*.Machine$double.eps
      # is exactly compatible with the conventions used in "Octave" or "Matlab".
      
      if (debug && method == "qr") {
        disp(qqr)
        disp(dim(L.full))
        disp(dim(vc))
        disp(vc)
      }
      if (debug) disp(L.full)
      if (debug) disp(vc)
      
      vv <-  L.full %*% vc %*% t(L.full)
      eta.hat <- L.full %*% beta
      Fstat <- (t(eta.hat) %*% qr.solve(vv,eta.hat,tol=1e-10)) / L.rank
      included.effects <- apply(L,2,function(x) sum(abs(x),na.rm=TRUE)) != 0
      denDF <- min( dfs[included.effects])
      numDF <- L.rank
      ret[[ii]]$anova <- list(numDF = numDF, denDF = denDF,
                              "F-value" = Fstat,
                              "p-value" = pf(Fstat, numDF, denDF, lower.tail = FALSE))
      ## Estimate
      
      etahat <- L %*% beta
      
      # NAs if not estimable:
      
      etahat[narows] <- NA
      if( nrow(L) <= maxrows ) {
        etavar <- L %*% vc %*% t(L)
        etasd <- sqrt( diag( etavar ))
      } else {
        etavar <- NULL
        etasd <- sqrt( apply( L * (L%*%vc), 1, sum))
      }
      
      denDF <- apply( L , 1 , function(x,dfs) min( dfs[x!=0]), dfs = dfs)
      
      aod <- cbind(
        Estimate=c(etahat),
        Std.Error = etasd,
        DF = denDF,
        "t-value" = c(etahat/etasd),
        "p-value" = 2*pt(abs(etahat/etasd), denDF, lower.tail =FALSE))
      colnames(aod)[ncol(aod)] <- 'p-value'
      if (debug ) disp(aod)
      if ( !is.null(clevel) ) {
        #print(aod)
        #print(aod[,'DF'])
        #print(aod[,'etasd'])
        hw <- qt(1 - (1-clevel)/2, aod[,'DF']) * aod[,'Std.Error']
        #print(hw)
        aod <- cbind( aod, LL = aod[,"Estimate"] - hw, UL = aod[,"Estimate"] + hw)
        #print(aod)
        if (debug ) disp(colnames(aod))
        labs <- paste(c("Lower","Upper"), format(clevel))
        colnames(aod) [ ncol(aod) + c(-1,0)] <- labs
      }
      if (debug ) disp(rownames(aod))
      aod <- as.dataf(aod)
      
      rownames(aod) <- rownames(as.dataf(L))
      labs(aod) <- names(dimnames(L))[1]
      ret[[ii]]$estimate <- aod
      ret[[ii]]$coef <- c(etahat)
      ret[[ii]]$vcov <- etavar
      ret[[ii]]$L <- L
      ret[[ii]]$se <- etasd
      ret[[ii]]$L.full <- L.full
      ret[[ii]]$L.rank <- L.rank
      if( debug ) disp(attr(Larg,'data'))
      data.attr <- attr(Larg,'data')
      if(is.null(data.attr) && !(is.null(data))) data.attr <- data
      ret[[ii]]$data <- data.attr
    }
    names(ret) <- names(Llist)
    attr(ret,"class") <- "wald"
    ret
  }

# Test
if(FALSE){
library(nlme)
fit <- lme(mathach ~ ses * Sex * Sector, hs, random = ~ 1|school)
summary(fit)
pred <- expand.grid( ses = seq(-2,2,1), Sex = levels(hs$Sex), Sector = levels(hs$Sector))
pred
wald(fit,model.matrix(fit,data=pred))
model.matrix(fit,data = pred)
model.matrix(~ ses * Sex * Sector,data=pred)
}
