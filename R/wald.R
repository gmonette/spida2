## wald.R
## This a collection of functions designed to facilitate testing hypotheses
## with Wald tests.
## The methods are readily extended to any fitting method that returns a vector
## of estimated parameters and an estimated variance matrix.
## Extensions are implemented through the 'getFix' generic function.
##
##
## see also Leff
## Changes:
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
#' a \code{\link{gefFix}} method for a \code{\link{glm}} object
#' and for 'mipo' objects in the packages 'mice':
#'
#' @param fit a model for which a \code{getFix} method exists.
#' @param Llist a hypothesis matrix or a pattern to be matched or a list of
#'        these
#' @param clevel level for confidence intervals. No confidence intervals if clevel is NULL
#' @param pred prediction data frame to evaluate fitted model using '
#'        \code{getX(fit) %*% coef}
#' @param data data frame used as 'data' attribute fot list elements returned only if
#'        the corresonding element of \code{Llist} has a NULL data attribute
#' @param debug (default FALSE) produce verbose information
#' @param maxrows maximum number of rows of hypothesis matrix for which a full
#'        variance-covariance matrix is returned
#' @param full if TRUE, the hypothesis matrix is the model matrix for
#'        \code{fit} such that the estimated coefficients are the predicted values for
#'       the fixed portion of the model. This is designed to allow the calculation of
#'       standard errors for models for which the \code{predict} method does not
#'       provide them.
#' @param fixed if \code{Llist} is a character to be used a regular expression,
#'       if \code{fixed} is TRUE \code{Llist} is interpreted literally, i.e.
#'       characters that have a special meaning in regular expressions are
#'       interpreted literally.
#' @param invert if \code{Llist} is a character to be used a regular
#'       expression, \code{invert == TRUE} causes the matches to be inverted so that
#'       coefficients that do not match will be selected.
#' @param method 'svd' (current default) or 'qr' is the method used to find the
#'       full rank version of the hypothesis matrix.  'svd' has correctly identified
#'       the rank of a large hypothesis matrix where 'qr' has failed.
#' @param pars passed to \code{\link{rstan::extract}} method for stanfit objects.
#' @param include passed to \code{\link{rstan::extract}} method for stanfit objects.#' 
#' @param help obsolete
#' @return An object of class \code{wald}, with the following components:
#'       COMPLETE
#' @seealso \code{\link{Lform}}, \code{\link{xlevels}}, \code{\link{dlevels}},
#'          \code{\link{Lfx}}, \code{\link{getX}}, \code{\link{M}},
#'          \code{\link{Lall}},\code{\link{Lc}},\code{\link{Lequal}},
#'          \code{\link{Ldiff}},\code{\link{Lmu}},\code{\link{Lmat}},\code{\link{Lrm}},
#'          \code{\link{Leff}}, \code{\link{as.data.frame.wald}}. To extend to new
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
wald <- 
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
        sv <- svd( na.omit(L) , nu = 0 )
        
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

#' @describeIn wald experimental version with RHS?
#' @export
wald2 <- function(fit, Llist = "",clevel=0.95, data = NULL, debug = FALSE , maxrows = 25,
                 full = FALSE, fixed = FALSE, invert = FALSE, method = 'svd',df = NULL, RHS = 0) {
#' GM: 2015 08 11:  to do:
#'  Experimental version of wald with RHS
#' NEEDS to be restructured with
#' 1. printing must show RHS
#' 2. needs to work with a list as second argument
#' 3. should redo handling of list so RHS is in list and so
#'    list handing is outside main function
#'

    if (full ) return( wald ( fit, model.matrix(fit)))
  dataf <- function(x,...) {
    x <- cbind(x)
    rn <- rownames(x)
    if( length( unique(rn)) < length(rn)) rownames(x) <- NULL
    data.frame(x,...)
  }
  as.dataf <- function(x,...) {
    x <- cbind(x)
    rn <- rownames(x)
    if( length( unique(rn)) < length(rn)) rownames(x) <- NULL
    as.data.frame(x,...)
  }

  unique.rownames <- function(x) {
    ret <- c(tapply(1:length(x), x , function(xx) {
      if ( length(xx) == 1) ""
      else 1:length(xx)
    })) [tapply(1:length(x),x)]
    ret <- paste(x,ret,sep="")
    ret
  }
  #      if(debug) disp( Llist)
  if( is.character(Llist) ) Llist <- structure(list(Llist),names=Llist)
  if(!is.list(Llist)) Llist <- list(Llist)

  ret <- list()
  fix <- getFix(fit)
  #      if(debug) disp(fix)
  beta <- fix$fixed
  vc <- fix$vcov

  dfs <- if(is.null(df) ) fix$df else df + 0*fix$df
  #      if(debug) disp(Llist)
  for (ii in 1:length(Llist)) {
    ret[[ii]] <- list()
    Larg <- Llist[[ii]]
    #          if(debug) {
    #               disp(ii)
    #               disp(Larg)
    #          }
    L <- NULL
    if ( is.character(Larg)) {
      L <- Lmat(fit,Larg, fixed = fixed, invert = invert)
    } else {
      if ( is.numeric(Larg)) {
        if ( is.null(dim(Larg))) {
          if(debug) disp(dim(Larg))
          if ( (length(Larg) < length(beta)) && (all(Larg>0)||all(Larg<0)) ) {
            L <- diag(length(beta))[Larg,]
            dimnames(L) <- list( names(beta)[Larg], names(beta))
          } else L <- rbind( Larg )
        }
        else L <- Larg
      }
    }
    #          if (debug) {
    #             disp(Larg)
    #             disp(L)
    #          }
    ## Delete coefficients that are NA
    Ldata <- attr( L , 'data')

    ## identify rows of L that are not estimable because they depend on betas that are NA
    Lna <- L[, is.na(beta), drop = FALSE]
    narows <- apply(Lna,1, function(x) sum(abs(x))) > 0

    L <- L[, !is.na(beta),drop = FALSE]
    attr(L,'data') <- Ldata
    beta <- beta[ !is.na(beta) ]

    ## Anova
    if( method == 'qr' ) {
      qqr <- qr(t(na.omit(L)))
      #Qqr <- Q(t(L))
      L.rank <- qqr$rank
      #L.rank <- attr(Qqr,'rank')
      #L.miss <- attr(Qqr,'miss')
      if(debug)disp( t( qr.Q(qqr)))
      L.full <- t(qr.Q(qqr))[ 1:L.rank,,drop=FALSE]
      #L.full <- t(Qqr[!L.miss,])[ 1:L.rank,,drop=F]
    } else if ( method == 'svd' ) {
      if(debug) disp(L)
      #              if(debug)disp( t(na.omit(t(L))))
      #              sv <- svd( t(na.omit(t(L))) , nu = 0 )
      sv <- svd( na.omit(L) , nu = 0 )

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
    if (debug) disp( vc )

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
    etahat <- L %*% beta-RHS
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

    aod <- cbind( Estimate=c(etahat),
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
    ret[[ii]]$RHS <- RHS
    ret[[ii]]$L <- L
    ret[[ii]]$se <- etasd
    ret[[ii]]$L.full <- L.full
    ret[[ii]]$L.rank <- L.rank
    if( debug ) disp(attr(Larg,'data'))
    ret[[ii]]$data <- attr(Larg,'data')
  }
  names(ret) <- names(Llist)
  attr(ret,"class") <- "wald"
  ret
}
#' Print method for wald objects
#'
#' @param x
#' @param round
#' @param pround
#' @param \dots
#' @return called for
#' @note This is a note
#' @author GM
#' @seealso \code{\link{wald}}
#' @examples
#' # coming soon
#' @export
print.wald <- function(x,round = 6, pround = 5,...) {
  pformat <- function(x, digits = pround) {
      x <- format(xx <- round(x,digits))
      x[ as.double(xx) == 0 ] <- paste(c("<.",rep('0',digits-1),'1'),collapse="")
      x
  }
  rnd <- function(x,digits) {
      if (is.numeric(x)) x <- round(x,digits=digits)
      format(x)
  }
  for( ii in 1:length(x)) {
    nn <- names(x)[ii]
    tt <- x[[ii]]
    ta <- tt$anova

    ta[["p-value"]] <- pformat(ta[["p-value"]])
    print(as.data.frame(ta,row.names=nn))
    te <- tt$estimate
    rowlab <- attr(te,"labs")

    te[,'p-value'] <- pformat( te[,'p-value'])
    if ( !is.null(round)) {
       for ( ii in 1:length(te)) {
           te[[ii]] <- rnd(te[[ii]],digits=round)
       }
    }
    labs(te) <- rowlab
    print(te,digits=round,...)
    cat("\n")
  }
  invisible(x)
}
#' Transform output of a Wald test into a data frame
#'
#' The generic method is included in the possibly false hope that it will ensure that the method
#' in this package is used.
#'
#' @param obj a 'wald' object
#' @param se a vector with the multiples of standard error used to generate lower and upper limits. 'names(se)'
#'        is appended to 'L' and 'U' to label the variables.
#' @param which selects elements of 'obj' to turn to a data.frame.
#' @return A data frame with estimated coefficient, standard error, and, optionally, upper and lower limits and
#'         the variables included the 'data' element of 'obj' if present.
#'         If \code{length(which) > 1}, the returned object is a list of data frames.
#' @examples
#' # coming soon!
#' @export
as.data.frame.wald <- function(obj, se = 2, digits = 3, sep = "", which = 1) {
  # modified by GM 2010_09_20 to avoid problems with coefs with duplicate rownames
  dataf <- function(x, ...) {
    x <- cbind(x)
    rn <- rownames(x)
    if(length(unique(rn)) < length(rn)) rownames(x) <- NULL
    data.frame(x, ...)
  }
  obj = obj[which]
  if (length(obj) == 1) { # e.g. is length(which) > 1
    cf <- obj[[1]]$coef
    ret <- data.frame(coef = cf, se = obj[[1]]$se)
    if(is.null(names(se))) names(se) <-
                 sapply(se, function(x) as.character(round(x, digits)))
    SE <- obj[[1]]$se
    SEmat <- cbind(SE) %*% rbind(se)
    cplus <- cf + SEmat
    cminus <- cf - SEmat
    colnames(cplus) <- paste("U",colnames(cplus),sep=sep)
    colnames(cminus) <- paste("L",colnames(cminus),sep=sep)
    ret <- cbind(ret, cplus, cminus)

    if(is.null(dd <- obj[[1]]$data)) return(ret)
    else return(cbind(ret, dd))
  }
  else ret <- lapply( obj, as.data.frame.wald)
  ret
}
#' Wald function producing a data frame for graphing
#'
#' A version of the wald function that produces a data frame directly, analogously to \code{as.data.frame(wald(...))}
#'
#' @inheritParams wald
#' @param se a vector with the multiples of standard error used to generate lower and upper limits. 'names(se)'
#'        is appended to 'L' and 'U' to label the variables.
#' @param which selects elements of 'obj' to turn to a data.frame.
#' @return A data frame with estimated coefficient, standard error, and, optionally, upper and lower limits and
#'         the variables included the 'data' element of 'obj' if present.
#'         If \code{length(which) > 1}, the returned object is a list of data frames.
#' @examples
#' \dontrun{
#' ###
#' ### Using walddf to create and plot a data frame with predicted values
#' ###
#'   data(hs)
#'   library( nlme )
#'   fit <- lme(mathach ~ (ses+I(ses^2)) * Sex * Sector, hs, random = ~ 1|school)
#'   summary(fit)
#'   pred <- expand.grid( ses = seq(-2,2,.1), Sex = levels(hs$Sex), Sector = levels(hs$Sector))
#'   head(pred)
#'   w <- walddf(fit, getX(fit,data=pred)) # attaches data to wald.object so it can included in data frame
#'   head(w)
#'   library(latticeExtra)
#'   xyplot(coef ~ ses | Sector, w, groups = Sex,
#'      auto.key = T, type = 'l',
#'      fit = w$coef,
#'      upper = w$L,
#'      lower = w$U,
#'      xlim = c(0,2),
#'      subscript = T) +
#'      glayer( gpanel.fit(...))
#'   wald(fit, 'ses')
#'   wald( fit, 'Sex')  # sig. overall effect of Sex
#'   wald( fit, ':Sex') # but no evidence of interaction with ses
#'   wald( fit, '\\^2') # nor of curvature
#'
#'   # conditional effect of ses
#'   head(getX(fit))
#'
#'   ###
#'   ###  Effect of ses: Differentiating with respect to ses
#'   ###
#'
#'   L <- Lfx(fit, list(
#'          0,
#'          1,
#'          2 * ses,
#'          0 * M(Sex),
#'          0 * M(Sector),
#'          1 * M(Sex),
#'          2 * ses * M(Sex),
#'          1 * M(Sector),
#'          2 * ses * M(Sector),
#'          0 * M(Sex) * M(Sector),
#'          1 * M(Sex) * M(Sector),
#'          2 * ses * M(Sex) * M(Sector)),
#'          pred)
#'   head(wald(fit, L), 20)
#'   w <- walddf(fit, L)
#'   xyplot(coef ~ ses | Sector, w, groups = Sex,
#'      auto.key = list(columns = 2, lines = T),
#'      type = 'l',
#'      fit = w$coef,
#'      upper = w$L,
#'      lower = w$U,
#'      xlim = c(0,2),
#'      ylab = 'estimate change in mathach per unit increase in ses',
#'      subscripts = T) +
#'      glayer( gpanel.fit(...)) +
#'      layer(panel.abline(a=0,b=0,lwd = 1, color ='black'))
#'
#'   ###
#'   ###  Difference in effect of ses between Sectors
#'   ###
#'
#'   pred.d <- expand.grid( ses = seq(-2,2,.1), Sex = levels(hs$Sex), Sector = levels(hs$Sector), Sector0 = levels(hs$Sector))
#'   pred.d <- subset(pred.d, Sector > Sector0)
#'   head(pred.d)
#'   L <- Lfx(fit, list(
#'          0,
#'          0,
#'          0 * ses,
#'          0 * M(Sex),
#'          0 * M(Sector),
#'          0 * M(Sex),
#'          0 * ses * M(Sex),
#'          1 * M(Sector) - M(Sector0),
#'          2 * ses * (M(Sector) - M(Sector0)),
#'          0 * M(Sex) * M(Sector),
#'          1 * M(Sex) * (M(Sector) - M(Sector0)),
#'          2 * ses * M(Sex) * (M(Sector) - M(Sector0))),
#'          pred.d)
#'   w <- walddf(fit, L)
#'   w
#'   w <- sortdf(w, ~ Sex/ses)
#'   xyplot(coef  ~ ses | Sex, w,
#'      type = 'l',
#'      auto.key = T,
#'      fit = w$coef,
#'      lower = w$L2,
#'      upper = w$U2,
#'      xlim = c(0,2),
#'      ylab = 'effect of ses (Public minus Catholic)',
#'      subscripts = T) +
#'      layer(gpanel.fit(...)) +
#'      layer(panel.abline(a=0,b=0,lwd = 1, color ='black'))
#'
#'
#' }
#' @export
walddf <- function(fit, Llist = "", clevel = 0.95,
                 data = NULL, debug = FALSE ,
                 full = FALSE, fixed = FALSE,
                 invert = FALSE, method = 'svd',
                 df = NULL,
                 se = 2, digits = 3, sep = '') {
  obj <- wald(fit = fit, Llist = Llist, clevel = clevel,
              data = data, debug = debug ,
              full = FALSE, fixed = FALSE,
              invert = FALSE, method = 'svd',
              df = NULL)
  ret <- yscs:::as.data.frame.wald(obj,
              se = se, digits = digits,
              sep = sep)
  ret
}

#' Extract coefficients from Wald object
#'
#' @param obj wald object
#' @param se (default: FALSE) include standard errors
#' @return coefficients and, optionally, standard errors
#' @author GM
#' @examples
#' # coming soon
#' @export
coef.wald <- function( obj , se = FALSE ) {
 if ( length(obj) == 1) {
    ret <-
    ret <- obj[[1]]$coef
    if ( is.logical(se) && (se == TRUE) ) {
       ret <- cbind( coef = ret, se = obj[[1]]$se)

    } else if ( se > 0 ){
       ret <- cbind( coef = ret, coefp = ret+se*obj[[1]]$se,
              coefm = ret - se*obj[[1]]$se)
       attr(ret,'factor') <- se
    }
 }
 else ret <- sapply( obj, coef.wald )
 ret
}

##
##
##   Functions to perform a GLH on lm, lme or lmer models
##   August 13, 2005
##
##
##
##   Lmat: generate a hypothesis matrix based on a pattern
##
##   glh
##   Lmat
##   Ldiff
##   getFix
##
##   print.glh
##

#' Get information on fixed effects from a model object
#'
#' \code{getFix} extracts coefficients, covariance matrix and degrees of
#' freedom from model objects. Its main purpose is to extract information need
#' by the \code{wald} function. To extend the wald function to a new class of
#' objects, it is merely necessary to write a method for \code{getFix}.
#'
#' Extending the \code{\link{wald}} function to a new class of objects only
#' requires a \code{getFix} method for the new class.
#'
#' @param fit A fitted model object
#' @param \dots Other arguments [unused]
#' @return Returns a list with the following components:
#'      \itemize{
#'           \item{fixed}{Fixed effect parameter estimates}
#'           \item{vcov}{Covariance matrix of the parameters}
#'           \item{df}{denominator degrees of freedom for each effect}
#'           }
#' @author Georges Monette
#' @seealso \code{\link{wald}}
#' @examples
#' library(nlme)
#' fit <- lme( mathach ~ (ses + I(ses^2)) * Sex, hs, random = ~ 1 + ses| school)
#' getFix(fit)
#'
#' data(Prestige, package="car")
#' mod.prestige <- lm(prestige ~ education + income, data=Prestige)
#' getFix(mod.prestige)
#'
#' @export
getFix <- function(fit,...) UseMethod("getFix")
#' @describeIn getFix method for multinom objects in package nnet
#' @export
getFix.multinom <- function(fit,...) {
  ret <- list()
  ret$fixed <- c(t(coef(fit)))
  ret$vcov <- vcov(fit)
  names(ret$fixed) <- rownames(ret$vcov)
  df <- nrow(na.omit(cbind(fit$residuals))) - length(ret$fixed)
  ret$df <- rep( df, length(ret$fixed))
  ret
}

#' @describeIn getFix method for lm objects
#' @export
getFix.lm <- function(fit,...) {
       ss <- summary(fit)
       ret <- list()
       ret$fixed <- coef(fit)
       ret$vcov <- ss$sigma^2 * ss$cov.unscaled
       ret$df <- rep(ss$df[2], length(ret$fixed))
       ret
}

#' @describeIn getFix method for glm objects
#' @export
getFix.glm <- function(fit,...) {
       ss <- summary(fit)
       ret <- list()
       ret$fixed <- coef(fit)
       ret$vcov <- vcov(fit)
       ret$df <- rep(ss$df.residual, length(ret$fixed))
       ret
}
#' @describeIn getFix method for lme objects in the nlme package
#' @export
getFix.lme <- function(fit,...) {
       require(nlme)
       ret <- list()
       ret$fixed <- nlme::fixef(fit)
       ret$vcov <- fit$varFix
       ret$df <- fit$fixDF$X
       ret
}
#' @describeIn getFix method for gls objects in the nlme package
#' @export
getFix.gls <- function(fit,...) {
  require(nlme)
  ret <- list()
  ret$fixed <-coef(fit)
  ret$vcov <- vcov(fit)
  ds <- fit$dims
  df <- ds[[1]] - sum( unlist( ds[-1]))
  ret$df <- rep(df, length(coef(fit)))
  ret
}

#' @describeIn getFix method for lmer objects in the lme4 package
#' @export
getFix.lmer <- function(fit,...) {
  # 2014 06 04: changed fit@fixef to fixef(fit)
  ret <- list()
       ret$fixed <- fixef(fit)
       ret$vcov <- as.matrix( vcov(fit) )
       # ret$df <- Matrix:::getFixDF(fit)
       ret$df <- rep( Inf, length(ret$fixed))
       ret
}
#' @describeIn getFix method for glmer objects in the lme4 package
#' @export
getFix.glmer <- function(fit,...) {
  # 2014 06 04: changed fit@fixef to fixef(fit)

  ret <- list()
       ret$fixed <- fixef(fit)
       ret$vcov <- as.matrix(vcov(fit))
       # ret$df <- Matrix:::getFixDF(fit)
       ret$df <- rep( Inf, length(ret$fixed))
       ret
}

#' @describeIn getFix method for mer objects in the lme4 package
#' @export
getFix.mer <- function(fit,...) {
  # 2014 06 04: changed fit@fixef to fixef(fit)

       ret <- list()
       ret$fixed <- fixef(fit)
       ret$vcov <- as.matrix(vcov(fit))
       # ret$df <- Matrix:::getFixDF(fit)
       ret$df <- rep( Inf, length(ret$fixed))
       ret
}
#' @describeIn getFix method for zeroinfl objects in the pscl?? package
#' @export
getFix.zeroinfl <- function(fit,...){
       ret <- list()
       ret$fixed <- coef(fit)
       ret$vcov <- as.matrix(vcov(fit))
       # ret$df <- Matrix:::getFixDF(fit)
       ret$df <- rep( Inf, length(ret$fixed))
       ret
}
#' @describeIn getFix method for mipo objects in the mice package
#' @export
getFix.mipo <- function( fit, ...){
  # pooled multiple imputation object in mice
  # uses the minimal df for components with non-zero weights
  # -- this is probably too conservative and should
  # improved
  ret <- list()
  ret$fixed <- fit$qbar
  ret$vcov <- fit$t
  ret$df <- fit$df
  ret
}
#' @describeIn getFix method for MCMCglmm objects in the MCMCglmm package
#' @export
getFix.MCMCglmm <- function(fit,...) {
  ret <- list()
  ret$fixed <- apply(fit$Sol, 2, mean)
  ret$vcov <- var( fit $ Sol)
  ret$df <- rep(Inf, length(ret$fixed))
  ret
}
#' @describeIn getFix method for `stanfit` objects in the `rstan` package
#' @export
getFix.stanfit <-
function(fit, pars, include = TRUE, ...) {
  if(missing(pars)) pars <- dimnames(fit)$parameter
  sam <- as.matrix(fit, pars = pars , include = include)
  ret <- list()
  ret$fixed <- apply(sam, 2, mean)
  ret$vcov <- var(sam)
  ret$df <- rep(Inf, length(ret$fixed))
  ret
}
#' @describeIn getFix print message if getFix id used for a class for which a method has not been written
#' @export
getFix.default <- function(fit, ...) stop(paste("Write a 'getFix' method for class",class(fit)))

#' Generic 'vcov' extended to objects with a \code{getFix} method
#'
#' @param fit
#' @return estimated variance covariance matrix of fixed effects
#' @author GM
#' @export
Vcov <- function(fit) {
     getFix(fit)$vcov
}

#' Correlation matrix of fixed effects
#'
#' @describeIn Vcov correlation matrix of fixed effects
#' @export
Vcor <- function(fit) {
     vc <- cov2cor(getFix(fit)$vcov)
     svds <- svd(vc)$d
     attribute(vc,'conditionNumber') <- svds[1]/svds[length(svds)]
     vc
}

#
# getFix.rdc <- function(fit, ...) UseMethod("getFix")
#
# getFix.rdc.lmer <- function(fit, ...) {
#             ret <- list()
#             ret$fixed <- fixef(fit)
#             ret$vcov <- vcov(fit)
#             ret$df <- Matrix:::getFixDF(fit)
#             ret
# }
#
# getFix.rdc.lm <- function(fit, ...) {
#           ret <- list()
#           ret$fixed <- coef(fit)
#           ret$vcov <- vcov(fit)
#           ret$df <- fit$df.residuals
#           ret
# }



#' Extract the model data frame for various fitting methods
#'
#' getData is to lme what model.frame is to lm
#' getData is implemented as a method in nlme
#' so we just add a method for 'lm' and other objects
#'
#' @param fit a fitted object
#' @export
getData <- function(x,...) UseMethod("getData")
#' @describeIn getData method for lmer objects
#' @export
getData.lmer <- function(x,...) slot(x,'frame')
#' @describeIn getData method for lme objects
#' @export
getData.lme <- function(x,...) nlme:::getData.lme(x,...)
#' @describeIn getData method for gls objects
#' @export
getData.gls <- function(x,...) nlme:::getData.gls(x,...)
#' @describeIn getData method for lm objects
#' @export
getData.lm <- function(x,...) model.frame(x,...)


#' get the names of variables that are factors
#'
#' @param object a data frame or object with a \code{\link{getData}} method that produces a data frame or a list
#' @param ... other orguments (currently not used)
#' @export
getFactorNames <- function(object, ...) UseMethod("getFactorNames")
#' @describeIn getFactorNames
#' @export
getFactorNames.data.frame <- function(object,...) {
     names(object)[ sapply(object, is.factor)  ]
}
#' @describeIn getFactorNames
#' @export
getFactorNames.default <- function(object,...) getFactorNames( getData(object))


#' Print method for 'cat' objects
#'
#' @param x
#' @param \dots
#' @return invisible(x)
#' @export
print.cat <- function(object,...) {
    cat(object)
    invisible(object)
}
#' Hypothesis matrix generated by expressions
#'
#' Creates an L matrix using expressions evaluated in 'data' for each column of
#' the L matrix
#'
#' If \code{Lform} is called with only a \code{fit} argument, it outputs code
#' consisting of an expression that would, if used as the \code{fmla} argument
#' to \code{Lform} would generate the full design matrix for the linear model.
#'
#' If \code{Lform} is called with two or three arguments, it generates a
#' hypothesis matrix by evaluating the expressions in \code{form} in the
#' environment \code{data}. The function \code{M} is designed to facilitate the
#' generation of blocks of the hypothesis matrix corresponding to main effects
#' or interaction effects of factors.
#'
#' \verb{ Creates a linear hypothesis
#' matrix, i.e. an L matrix, using formulas evaluated in 'data' for each column
#' of the L matrix. This approach lends itself to creating hypotheses and
#' estimates based on data such as partial derivatives with respect to
#' variables evaluated at each data point.
#'
#' An example is the estimation of growth rates in longitudinal models.
#'
#' library(car) library(spida) fit <- lm( income ~ (education + I(education^2)
#' )* type, Prestige) summary(fit)
#'
#' . . .  Coefficients: Estimate Std. Error t value Pr(>|t|) (Intercept) 891.3
#' 23889.1 0.037 0.97032 education 210.0 5638.8 0.037 0.97037 I(education^2)
#' 38.3 328.3 0.117 0.90740 typeprof 191523.2 63022.0 3.039 0.00312 ** typewc
#' 25692.3 73888.0 0.348 0.72887 education:typeprof -28133.0 10236.0 -2.748
#' 0.00725 ** education:typewc -4485.4 14007.9 -0.320 0.74956
#' I(education^2):typeprof 1017.5 451.8 2.252 0.02679 * I(education^2):typewc
#' 170.9 671.6 0.255 0.79967 . . .
#'
#' # estimate the marginal value of occupation for each occupation in the data
#' set
#'
#' L <- list( 'marginal value of education' =Lform( fit, form = list(
#' 0,1,2*education, 0,0, type == 'prof', type == 'wc', 2*education
#' *(type=='prof'), 2* education * (type == 'wc')), data = Prestige)) wald(
#' fit, L ) chat <- coef( wald( fit, L ), se = 2) xyplot( coef +coefp+coefm ~
#' education | type, cbind(Prestige,chat)[order(Prestige$education),], type =
#' 'l') xyplot( chat~ education | type, Prestige) }
#'
#' @param fit a fitted model with a 'getFix' method.
#' @param expr.list a list of expressions with one component for each column
#'    (or groups of columns) of the hypothesis matrix corresponding to each term
#'    of the model. A term with multiple degrees of freedom can either be
#'    generated as separate single terms or with an expression that evaluates to a
#'    suitable matrix.
#' @param data the data frame in which expressions are evaluated.
#' @param formula as an argument of \code{M}, a one-sided formula defining a
#' main effect or an interaction term involving factors.
#' @param expr as an argument of \code{M}, an expression that can be evaluated
#' in each row of \code{data} to form element(s) of the corresponding row of
#' \code{L}.  formula defining a main effect or an interaction term involving
#' factors.
#' @return hypothesis matrix
#' @seealso \code{\link{wald},\link{Leff}} and \code{\link{Lfx}} for a new
#' improved but experimental version.
#' @examples
#' \dontrun{
#'       library(car)
#'       mod <- lm( income ~ (education + I(education^2) )* type, Prestige)
#'       summary(mod)
#'
#'       # estimate the marginal value of an extra year of education for a
#'       # range of years for each type
#'
#'       years.type <- expand.grid( education = seq(6,18,2), type = levels(Prestige$type))
#'       Lf <- Lform( mod,
#'          list( 0, 1, 2*education, 0, 0, type =="prof", type =="wc",
#'             2*education*(type =="prof"), 2*education*(type =="wc")),
#'          years.type)
#'       Lf
#'       ww <- wald( mod, Lf)
#'       ww
#'       ytderiv <- as.data.frame( ww, se = 2)
#'       head( ytderiv )
#'       xyplot( coef ~ education, ytderiv, groups = type, type = 'l',
#'               auto.key = list(columns = 3, lines = TRUE, points = FALSE)
#' }
#' @export
Lform <- function( fit, form, data = getData(fit)) {
 # 2011-12-01: replaced with version below
  # 2012 12 04
  # Plan for Lform
  #

  # 2012 12 05: Lform becomes Lex to acknowledge the fact that it uses
  # expressions instead of formulas
    if (missing(form)) return ( Lcall(fit))
      gg <- getFix(fit)
      Lsub <- do.call(cbind,eval( substitute( form ), data))
      if( (nrow(Lsub) != nrow( data))) {
          if ((nrow(Lsub)==1)) Lsub <- Lsub[rep(1,nrow(data)),]
          else stop('nrow(Lsub) != nrow(data)')
      }
      if( is.null( colnames(Lsub))) colnames(Lsub) <- rep('',ncol(Lsub))
      L <- matrix( 0, nrow = nrow(Lsub), ncol = length( gg$fixed))
      rownames(L) <- rownames(data)
      colnames(L) <- names( gg$fixed)
      Lpos <- Lsub[, colnames(Lsub) == '', drop = FALSE]
      # disp(Lpos)
      Lnamed <- Lsub[ , colnames(Lsub) !='', drop  = FALSE]
      # disp(Lnamed)
      for ( ip in seq_len( ncol( Lpos ))) L[,ip] <- Lpos[,ip]
      if ( ncol( Lnamed ) > 0 ) {
         if ( length( unknown <- setdiff( colnames(Lnamed) , colnames(L)))) {
              stop( paste("Unknown effect(s):" , unknown, collapse = " "))
         }
         for ( nn in colnames(Lnamed)) L[,nn] <- Lnamed[,nn]
      }
      attr(L,"data") <- data
      L
}
#' Generate a hypothesis matrix
#'
#' Generates a hypothesis matrix to test whether a group of coefficients in a
#' linear model are jointly zero, selecting the coefficients that match a
#' regular expression.
#'
#' @param fit a fitted object with a \code{\link{getFix}} method.
#' @param pattern a regular expression that matches names of coefficients
#' @return hypothesis matrix to test the hypothesis that the true values of
#'        matched coefficients are zero
#' @export
Lmat <- function(fit, pattern, fixed = FALSE, invert = FALSE, debug = FALSE) {
     # pattern can be a character used as a regular expression in grep
     # or a list with each component generating  a row of the matrix
     umatch <- function( pat, x, fixed , invert ) {
            ret <- rep(0,length(pat))
            for ( ii in 1:length(pat)) {
                imatch <- grep(pat[ii], x, fixed= fixed, invert = invert)
                if ( length(imatch) != 1) {
                   cat("\nBad match of:\n")
                   print(pat)
                   cat("in:\n")
                   print(x)
                   stop("Bad match")
                }
                ret[ii] <- imatch
            }
            ret
     }
     if ( is.character(fit)) {
        x <- pattern
        pattern <- fit
        fit <- x
     }
     fe <- getFix(fit)$fixed
     ne <- names(fe)
     if (is.character(pattern)) {
        L.indices <- grep(pattern,names(fe), fixed = fixed, invert = invert)
        ret <- diag( length(fe)) [L.indices,,drop = FALSE]
        if (debug) disp(ret)
        rownames(ret) <- names(fe) [L.indices]
        labs(ret) <- "Coefficients"
     } else if (is.list(pattern)){
        ret <- matrix(0, nrow = length(pattern), ncol = length(fe))
        colnames(ret) <- ne
        for ( ii in 1:length(pattern)) {
            Lcoefs <- pattern[[ii]]
            pos <- umatch(names(Lcoefs), ne, fixed = fixed, invert = invert)
            if ( any( is.na(pos))) stop("Names of L coefs not matched in fixed effects")
            ret[ii, pos] <- Lcoefs
        }
        rownames(ret) <- names(pattern)
      }
      labs(ret) <- "Coefficients"
      ret
}
#' Older version of Ldiff
#'
#' @param fit
#' @param pat
#' @param levnames
#' @param reflevel
#' @param cut
#' @return hypothesis matrix
#' @export
Ldiff.old <- function(fit, pat, levnames = c(reflevel,substring(rownames(L),cut+1)),
                       reflevel = "<ref>", cut = nchar(pat)) {
         L <- Lmat(fit, pat)
         nam <- rownames(L)
         n <- nrow(L)
         if(n < 2) return(L)
         plus <- unlist( apply( rbind( 2:n), 2, seq, n))
         minus <- rep(1:(n-1), (n-1):1)
         Lp <- L[ plus, ]
         Lm <- L[ minus, ]
         Lret <- rbind( L, Lp - Lm)
         rn <- paste( levnames[ c(1:n,plus) + 1], levnames[ c(rep(0,n), minus)+1], sep = " - ")
         rownames(Lret) <- rn
         Lret
}
#' Version of Ldiff used in RDC
#'
#' @param fit
#' @param nam
#' @param ref
#' @return hypothesis matrix
#' @export
Ldiff.rdc <- function( fit, nam , ref = "no longer used") {
      # based on Lmu
      # Tests all pairwise difference in factor with model with Intercept term
      Lm <- Lmu(fit, nam)
      levs <- rownames(Lm)
      n <- nrow(Lm)
      if (n < 2) return (Lm)
      plus <- unlist( apply ( rbind(2:n), 2, seq, n))
      minus <- rep(1:(n-1), (n-1):1)
      Lret <- Lm[plus,] - Lm[minus,]
      rn <- paste( levs [plus], levs[minus] , sep = " - ")
      rownames(Lret) <- rn
      Lret
}

#' Hypothesis matrix to test differencs in factor levels
#'
#' @param fit
#' @param pat
#' @param levnames
#' @param reflevel
#' @param cut
#' @param verbose
#' @return hypothesis matrix
#' @export
Ldiff <- function( fit, pat, levnames = c(reflevel,substring(rownames(L),cut+1)),
         reflevel ="<ref>", cut=nchar(pat),verbose=F) {
      L <- Lmat(fit, paste("^",pat,sep=""))
      nam <- rownames(L)
      n <- nrow(L)
      zm <- matrix(1:n,nrow=n,ncol=n)
      plus <- zm[ col(zm) < row(zm)]
      minus <- rep(1:(n-1), (n-1):1)
      Lp <- L[plus,]
      Lm <- L[minus,]
      Lret <- rbind( L, Lp - Lm)
         pnames <- levnames [ c(1:n, plus) +1]
      mnames <- levnames [ c(rep(0,n), minus) + 1]
      if (verbose) {
        print(levnames)
        print(plus)
        print(minus)
        print(Lp)
        print(Lm)
        print(L)
        print(Lret)
        print(pnames)
        print(mnames)
      }
      rn <- paste( levnames[ c(1:n,plus)+1], levnames[ c(rep(0,n),minus) + 1], sep = " - ")
      rownames(Lret) <- rn
      Lret
}

#' Estimate predicted response for a factor level.
#'
#' @param fit
#' @param nam
#' @param verbose
#' @export
Lmu <- function(fit, nam, verbose = 0) {
       ## "Works only if 'nam' is a factor and a main effect and model has Intercept")
       if ( class(fit) != 'lmer' ) stop( "only implemented for lmer")
       v <- fit@frame[[nam]]
       if( !is.factor(v)) stop ("nam needs to specify the name of a factor")
       levs <- levels(v)
       if( verbose > 0) print(levs)
       cmat <- contrasts(v)
       if( verbose > 0) print(cmat)
       #  print(cmat)
       fe <- getFix(fit)$fixed
       if( verbose > 0) print(fe)
       if ( substring(nam,1,1) != '^') nam <- paste("^",nam,sep="")
       L.indices <- grep(nam,names(fe))
       if( verbose > 0) print(L.indices)
       L <- matrix(0,nrow=length(levs), ncol = length(fe))

       colnames(L) <- names(fe)
       if( verbose > 0) print(L)
          rownames(L) <- levs
       L[,L.indices] <- cmat
       if('(Intercept)' %in% colnames(L)) L[,'(Intercept)'] <- 1
       L
}

#' Hypothesis matrix for lmer objects: comparisons with reference level
#'
#' @param fit
#' @param nam
#' @param ref
#' @param verbose
#'
#' @export
Lc <- function(fit, nam, ref = 1, verbose = 0) {
       ## Comparisons with one level
       ## Use Lmu
       ## "Works only if 'nam' is a factor and a main effect and model has Intercept?")
       if ( class(fit) != 'lmer' ) stop( "only implemented for lmer")
       L <- Lmu( fit, nam)
       Lref <- L[ ref,,drop = FALSE]
       index <- 1:nrow(L)
       names(index) <- rownames(L)
       refind <- index[ref]
       if (length(refind) != 1) stop( paste( ref, "does not refer to a single level"))
       Lret <- L[-refind,]
       Lret <- Lret - cbind( rep(1,nrow(Lret))) %*% Lref
       attr(Lret,"heading") <- paste("Comparisons with reference level:", rownames(L)[refind])
       Lret
}

#' Construct hypothesis matrix to test repeated measures factor effects
#'
#' @param fit
#' @param nam
#' @param vals
#' @return hypothesis matrix
#' @export
Lrm <- function(fit, nam, vals = 1:nrow(L.mu)) {
    ## Repeated measures polynomial contrasts
    ## Uses Lmu
    ## "Works only if 'nam' is a factor and a main effect and model has Intercept?")
    ##
    L.mu <- Lmu(fit, nam)
    # print(L.mu)
    pp <- cbind( 1, Poly(vals, nrow(L.mu) -1))
    # print(pp)
    ortho <- Q(pp)[,-1] # (linear, quad, etc.)
    # print(ortho)
    ortho <- ortho[,-1]
    maxp <- max( 5, nrow(L.mu))
    colnames(ortho) <- c('linear','quadratic','cubic', paste("^",4:maxp,sep=''))[1:ncol(ortho)]
    L <- t(ortho) %*% L.mu
    L
}
# Lrm(fit, "SRH94")

#' Construct hypothesis matrix to test ????
#'
#' @param fit
#' @param factors
#' @return a hypothesis matrix
#' @export
Lcall <- function( fit , factors = getFactorNames(fit), debug = F){

      nams <- names(getFix(fit)$fixed)

      nams <- gsub( "^", ":", nams)   # delineate terms
      nams <- gsub( "$", ":", nams)   # delineate terms
      for ( ff in factors)   {
          ff.string <- paste( ff, "([^:]*)" , sep = '')
          if(debug) disp( ff.string)
          ff.rep <- paste(ff, " == \\'\\1\\'", sep = '')
          if(debug) disp(ff.rep)
          nams <- gsub( ff.string, ff.rep, nams)
      }
     # for ( ii in seq_along(matrix)) {
     #     mm.all   <- paste( "(:",names(matrix)[ii], "[^\\)]*\\))",sep='')
     #     mm.match <- paste( "(",names(matrix)[ii], "[^\\)]*\\))",matrix[ii], sep ='')
     #     mm.rep   <- paste( "\\1")
     #     which.null <- grepl( mm.all, nams) mm.null  <-
     #
     # }
      nams <- sub("(Intercept)", 1, nams)
      nams <- gsub( "^:","(",nams)
      nams <- gsub( ":$",")",nams)
      nams <- gsub( ":", ") * (", nams)
      #if(comment) nams <- paste( nams, "  #",nams)
      nams <- paste( "with (data, \n cbind(", paste( nams, collapse = ",\n"), ")\n)\n", collapse = "")
      class(nams) <- 'cat'
      nams
}

#' Hypothesis matrix to test equality of factor level effects
#'
#' @param fit
#' @param pat
#' @return hypothesis matrix
#' @export
Lequal <- function(fit, pat) {
       # Test for equality within levels of pat using all differences
         L <- Lmat(fit, pat)
         nam <- rownames(L)
         n <- nrow(L)
         if(n < 2) return(L)
         plus <- unlist( apply( rbind( 2:n), 2, seq, n))
         minus <- rep(1:(n-1), (n-1):1)
         Lp <- L[ plus, ]
         Lm <- L[ minus, ]
         Lret <- rbind( Lp - Lm)
         rn <- paste( nam[plus], nam[minus], sep = " - ")
         rownames(Lret) <- rn
         Lret
}




#Lc <- function(fit, vec ){
#   fe <- getFix(fit)$fixed
#   ret <- 0 * fe
#   if ( is.null(names(vec))) ret[] <-
#}
# Lmu(fit,"SRH")

#' Hypothesis matrix to test for lmer objects
#'
#' @param fit
#' @param pat
#' @return hypothesis matrix
#' @export
Lall <- function( fit , nam ) {
        if ( class(fit) != 'lmer' ) stop( "only implemented for lmer")
        v <- fit@frame[[nam]]
        if( !is.factor(v)) stop ("nam needs to specify the name of a factor")
        lev0 <- levels(v)[1]
        ret <-list()
        namf <- nam
        if ( substring(namf,1,1) != "^") namf <- paste("^", namf, sep ="")
        ret[[ nam ]] <- Lmat( fit, namf)
        ret[[ paste(nam,"mu",sep = '.') ]] <- Lmu(fit,nam)
        ret[[ paste(nam,"diff",sep = '.') ]] <- Ldiff( fit , nam)
        ret

}


#' @describeIn wald transforms estimates canfidence intervals using delta method
#' @export
wald.transform <- function (x, fun, label = "Transformed coefficients") {
  # transforms estimates and confidence intervals for wald, uses delta method
  # for se
  ant <- x[[1]]$estimate
  coefs <- as.double(ant$Estimate)
#   disp(coefs)
  trs <- numericDeriv( quote(fun(coefs)), 'coefs')
#   disp(trs)
  ant$Estimate <- c(trs)
  derivs <- abs( diag(as.matrix(attr(trs,'gradient'))))
  ant[["Std.Error"]] <- derivs * ant[["Std.Error"]]
  low.ind <- grep("^Lower ", colnames(ant))
  up.ind <- grep("^Upper ", colnames(ant))
  ant[[low.ind]] <- fun(ant[[low.ind]])
  ant[[up.ind]] <- fun(ant[[up.ind]])
  attr(ant,'labs') <- label
  x <- x[1]
  x[[1]]$estimate <- ant
  class(x) <- 'wald'
  x
}
#' @describeIn wald same as \code{wald}, kept for backward compatibility
#' @export
glh <- function( ...) wald( ...)    # previous name for 'wald' function

#' Design matrix for a fit possibly on a new data frame
#'
#' This function return the X matrix for a fit possibly based on a different data frame than the model.
#' It performs a function very close to that of \code{\link{model.matrix}} except that \code{\link{model.matrix}}
#' expects the variable for the LHS of the formula to be in the data set, ostensibly in order to remove rows
#' for which the LHS variable(s) are NA. In addition, \code{getX} attaches the argument data set as an attribute.
#'
#' Extending \code{getX} to new classes merely requires a \code{\link{getData}} method. The \code{\link{formula}}
#' method is also used but usually already exitsts.
#'
#' @param fit a fitted object with formula method
#' @param data (default NULL) a data frame on which to evaluate the design matrix
#' @return a design matrix
#' @export
getX <- function(fit, data = getData(fit)) {
  f <- formula(fit)
  if(length(f) == 3) f <- f[-2]
  ret <- model.matrix(f, data = data)
  attr(ret,'data') <- data #include data as attribute
  ret
}
# if(FALSE){ #TESTS:
#   library(nlme)
#   fit <- lme(mathach ~ ses * Sex * Sector, hs, random = ~ 1|school)
#   summary(fit)
#   getX(fit)
#   pred <- expand.grid( ses = seq(-2,2,1), Sex = levels(hs$Sex), Sector = levels(hs$Sector))
#   pred
#   getX(fit,pred)
#   w <- wald(fit,getX(fit,data=pred))
#   w
#   as.data.frame(w)
#   g <- getX(fit,pred)
#   attr(g,'data') <- pred
#   w <- wald(fit,g)
#   w
#   as.data.frame(w)
# }
