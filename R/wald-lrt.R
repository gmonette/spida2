## see wald-lrt in R2 for version saved before removing all the
## duplicated function that are in wald.R
## 
## 
## 
## Combines wald-lrt.R by Tino Ntentes
## with some new functions to compute kernels 
## conjugate complements
##
## wald-lrt.R
#' Kernel (null space) of a linear transformation
#' 
#' @param L a matrix
#' @param tol the smallest singular value that is considered to
#'        indicate a non-zero dimension.
#' @export
ker <- function(L, tol = 1e-14) {
  sv <- svd(t(L), nu = NCOL(L), nv = 0)
  ret <- sv$u[,sv$d < tol]
  attr(ret, 'd') <- sv$d  # diagnostic for close callss
  if(ncol(ret) > 0) ret else matrix(0,NCOL(L),1)
}
#' Canonical correlation on uncentered matrices
#'
#' The vector of canonical correlations of uncentered
#' matrices helps to identify whether the column
#' spaces of two matrices are identical, or whether one
#' is a subspace of the other, or whether the spaces
#' are orthogonal. 
#' 
#' @param x,y two matrices with the same number of rows
#' @export  
ccor <- function(x,y) {
  # canonical correlations without centering
  # to test whether two
  # matrices span the same or orthogonal column spaces,
  # also the dimension of the intersection
  # is the number of 1's
  cancor(x,y,F,F)$cor
}
#' SVD to generate an othonormal basis
#'
#' @param x a matrix 
#' @param tol the smallest singular value that is considered to
#'        indicate a non-zero dimension.
#' @return an orthonormal basis for the column span of x
#' @export
ob <- function(x,  tol = 1e-14) {
  # orthogonal basis
  sv <- svd(x, nv = 0)
  disp(sum(sv$d > tol))
  ret <- sv$u[,seq_len(sum(sv$d > tol))]
  attr(ret, 'd') <- sv$d
  ret
} 
#' Conjugate complement of span(X) in span(Z) with respect 
#' to inner product ip
#' 
#' @param X	matrix defining space whose complement is computed. 
#'          Not necessarily of full column rank
#' @param Z	matrix defining space within which the complement 
#'          is computed. 
#'          Should be of full column rank. Default: diag(nrow(X))
#' @param ip positive definite matrix defining inner product with 
#'        respect to which complement is computed. Default: diag(nrow(X))
#' @param tol tolerance (default 1e-14)
#' @examples
#' #' cc_svd(cbind(1:3))
#' cc_svd(cbind(1:3), cbind(1:3,1))
#' cc_svd(cbind(1:3), cbind(3*c(1,-2,1), 1:3,1))
#' # challenge with large numbers
#' K <- cc_svd(cbind(1,1:10,(1:10)^2))
#' K
#' ret <- list()
#' for(i in 1:20){
#'   cbind(1,(1:10-10^i)^2,1:10) %>% 
#'     { list(svd = ccor(cc_svd(.), K) -1, qr =ccor(cc_qr(.), K) -1)
#'     }  %>% lapply(abs) %>% 
#'     {list(med=lapply(.,median), ave = lapply(.,mean))} %>% unlist->ret[[i]]
#' }
#' ret <- do.call(rbind, ret)
#' ret
#' library(lattice)
#' library(latticeExtra)
#' xyplot(ts(log2(ret)),  scales = list(y = 'same'), ylab = 'cancorr: log base 2 of deviation from 1')
#' #
#' # challenge with small numbers
#' #
#' i <- 1
#' K <- cc_svd(cbind(1,1:10,(1:10 == 1)* 10^{-i}))
#' K
#' ret <- list()
#' for(i in 1:30){
#'   cbind(1,1:10,(1:10 == 1)* 10^{-i}) %>% 
#'     { list(svd = ccor(cc_svd(.), K) -1, qr =ccor(cc_qr(.), K) -1)
#'     }  %>% lapply(abs) %>% 
#'     {list(med=lapply(.,median), ave = lapply(.,mean))} %>% unlist->ret[[i]]
#' }
#' ret <- do.call(rbind, ret)
#' ret
#' 
#' xyplot(ts(log2(ret)),  scales = list(y = 'same'), ylab = 'cancorr: log base 2 of deviation from 1')
#' 
#' #
#' # challenge with small numbers added to 1
#' #
#' i <- 1
#' K <- cc_svd(cbind(1,1:10,1+(1:10 == 1)* 10^{-i}, (1:10 != 1)))
#' K
#' ret <- list()
#' for(i in 1:30){
#'   cbind(1,1:10, 1 + (1:10 == 1)* 10^{-i}) %>% 
#'     { list(svd = ccor(cc_svd(.), K) -1, qr =ccor(cc_qr(.), K) -1)
#'     }  %>% lapply(abs) %>% 
#'     {list(med=lapply(.,median), ave = lapply(.,mean))} %>% unlist->ret[[i]]
#' }
#' ret <- do.call(rbind, ret)
#' ret
#' @export
cc_svd <- function( X , Z = NULL , ip = NULL, tol = 1e-14 ) {
  #
  # conjugate complement using the SVD
  tX <- t(X)
  if(!is.null(ip)) tX <- tX %*% ip
  if(!is.null(Z)) tX <- tX %*% Z
  K <- ker(tX, tol = tol)
  if(!is.null(Z)) Z %*% K else K
}
#' 
#' @describeIn cc_svd Conjugate complement using QR algorithm
#' @export
cc_qr <- function( X , Z = diag( NROW(X)) , ip = diag( NROW(X)), tol = 1e-14 ) {
  # conjugate complement using QR
  # keep versions consistent
  # ConjComp returns a basis for the conjugate complement of the
  # conjugate projection of X into span(Z) with respect to inner product with
  # matrix ip.
  # Note: Z is assumed to be of full column rank but not necessarily X.
  nr <- NROW(X)
  if(!is.null(ip)) X <- ip %*% X
  if(!is.null(Z)) X <- t(Z) %*% X
  xq <- qr(X, tol = tol)
  if ( xq$rank == 0 ) return(if(!is.null(Z)) Z else diag(nr)) 
  a <- qr.Q( xq, complete = TRUE ) [ ,-seq_len(xq$rank)]
  if(!is.null(Z)) Z %*% a else a
}


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
#' @param LRT if TRUE, provide likelihood ratio test statistic and p-value
#' @param pred prediction data frame to evaluate fitted model using
#'        `getX(fit) %*% coef`
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
wald_lrt <- 
  function(fit, Llist = "", clevel = 0.95,
           LRT = TRUE, pred = NULL,
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
        if ( tol.fac > 1e6 ) warning("Poorly conditioned L matrix, calculated numDF may be incorrect")
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
                              "Wald F-value" = Fstat,
                              "Wald p-value" = pf(Fstat, numDF, denDF, lower.tail = FALSE))
      ## LRT
      if (LRT) {
        
        if ( class(fit) %in% c("lm", "glm") ) {
          model_mat <- getX(fit)
          sv2 <- svd(na.omit(L) , nu = 0, nv = NCOL(L))
          rmv <- if (numDF == 0) 1:NCOL(L) else -(1:numDF)
          constrainedX <- as.matrix(model_mat %*% sv2$v[ , rmv, drop=FALSE])
          numCol <- NCOL(constrainedX)
          names1 <- as.character(1:numCol)
          names1 <- sapply(names1, function(x) paste0("XConCol", x))
          newData <- data.frame(cbind(constrainedX, getData(fit)))
          names(newData) <- c(names1, names(getData(fit)))
          constrained_fit <- try(update(fit, as.formula(paste0(". ~ ", paste(names1, collapse = "+"),  "- 1")), 
                                        data = newData))
          if (is.null(constrained_fit) | class(constrained_fit) == "try-error") {
            warning(paste0(class(fit), " object from original fit threw an error with ML method, LRT will be suppressed."))
          } else {
            changeDF <- NCOL(L) - numDF
            lrt_stat <- 2 * ( logLik(mlfit) - logLik(constrained_fit) )
            if (lrt_stat < 0) warning("LRT stat is negative, original fit may have convergence problems.")
            ret[[ii]]$anova[["changeDF"]] <- changeDF
            ret[[ii]]$anova[["LRT ChiSq"]] <- lrt_stat
            ret[[ii]]$anova[["LRT p-value"]] <- pchisq(lrt_stat, changeDF, lower.tail = FALSE)
          }
          
        } else if ( class(fit) %in% c("lme", "gls") ) {
          mlfit <- try(update(fit, method = "ML"))   # refit both models using 'ML' for proper comparison
          if (!(is.null(mlfit) | class(mlfit) == "try-error")) {
            model_mat <- getX(fit)
            sv2 <- svd(na.omit(L) , nu = 0, nv = NCOL(L))
            rmv <- if (numDF == 0) 1:NCOL(L) else -(1:numDF)
            constrainedX <- as.matrix(model_mat %*% sv2$v[ , rmv, drop=FALSE])
            numCol <- NCOL(constrainedX)
            names1 <- as.character(1:numCol)
            names1 <- sapply(names1, function(x) paste0("XConCol", x))
            newData <- data.frame(cbind(constrainedX, getData(fit)))
            names(newData) <- c(names1, names(getData(fit)))
            constrained_fit <- try(update(fit, fixed. = as.formula(paste0(". ~ ", paste(names1, collapse = "+"),  "- 1")), 
                                          data = newData, method = "ML"))
            if (is.null(constrained_fit) | class(constrained_fit) == "try-error") {
              warning(paste0(class(fit), " object from original fit threw an error with ML method, LRT will be suppressed."))
            } else {
              changeDF <- NCOL(L) - numDF
              lrt_stat <- 2 * ( logLik(mlfit) - logLik(constrained_fit) )
              if (lrt_stat < 0) warning("LRT stat is negative, original fit may have convergence problems.")
              ret[[ii]]$anova[["changeDF"]] <- changeDF
              ret[[ii]]$anova[["LRT ChiSq"]] <- lrt_stat
              ret[[ii]]$anova[["LRT p-value"]] <- pchisq(lrt_stat, changeDF, lower.tail = FALSE)
            }
          } else {
            warning(paste0(class(fit), " object from original fit threw an error with ML method, LRT will be suppressed."))
          }
          
        } else if ( class(fit) %in% c("lmer", "lmerMod", "glmer", "glmerMod", "lmerModLmerTest") ) {
          mlfit <- try(update(fit, REML = FALSE))    # refit both models using 'ML' for proper comparison
          if (!(is.null(mlfit) | class(mlfit) == "try-error")) {
            model_mat <- model.matrix(fit)
            sv2 <- svd(na.omit(L) , nu = 0, nv = NCOL(L))
            rmv <- if (numDF == 0) 1:NCOL(L) else -(1:numDF)
            constrainedX <- as.matrix(model_mat %*% sv2$v[ , rmv, drop=FALSE])
            numCol <- NCOL(constrainedX)
            names1 <- as.character(1:numCol)
            names1 <- sapply(names1, function(x) paste0("XConCol", x))
            newData <- data.frame(cbind(constrainedX, getData(fit)))
            names(newData) <- c(names1, names(getData(fit)))
            formulalmer <- formula(fit)
            formulalmer <- paste(deparse(formulalmer, width.cutoff = 500), collapse="")
            formind <- gregexpr("\\({1}?[^\\(|\\||\\)]*\\|{1,2}[^//)]*[\\)]{1}?", formulalmer)[[1]]
            lengs <- attr(formind, which = "match.length")
            ranEffects <- substring(formulalmer, first = formind, last = formind + lengs - 1)
            constrained_fit <- try(update(fit, as.formula(paste0(". ~ ", paste(names1, collapse = "+"),  
                                                                 "- 1 + ", paste(ranEffects, collapse = "+"))),
                                          data = newData, REML = FALSE))
            if (is.null(constrained_fit) | class(constrained_fit) == "try-error") {
              warning(paste0(class(fit), " object from original fit threw an error with ML method, LRT will be suppressed."))
            } else {
              changeDF <- NCOL(L) - numDF
              lrt_stat <- 2 * ( logLik(mlfit) - logLik(constrained_fit) )
              if (lrt_stat < 0) warning("LRT stat is negative, original fit may have convergence problems.")
              ret[[ii]]$anova[["changeDF"]] <- changeDF
              ret[[ii]]$anova[["LRT ChiSq"]] <- lrt_stat
              ret[[ii]]$anova[["LRT p-value"]] <- pchisq(lrt_stat, changeDF, lower.tail = FALSE)
            }
          } else {
            warning(paste0(class(fit), " object from original fit threw an error with ML method, LRT will be suppressed."))
          }
        } else {
          message(paste0("LRT not yet tested with class ", class(fit)))
        }
      }
      
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
    attr(ret,"class") <- c("wald_lrt", 'wald')
    ret
  }

# Test
if(FALSE){
  library(nlme)
  fit <- lmer(mathach ~ ses * Sex * Sector + (1 | school), hs)
  summary(fit)
  pred <- expand.grid( ses = seq(-2,2,1), Sex = levels(hs$Sex), Sector = levels(hs$Sector))
  pred
  wald(fit, 'Sector')
  model.matrix(fit,data = pred)
  test1 <- model.matrix(~ ses * Sex * Sector,data=pred)
}

#' Print method for wald objects
#'
#' @param x wald object
#' @param round FIXME
#' @param pround FIXME
#' @param \dots FIXME
#' @return called for  FIXME
#' @note This is a note
#' @author GM
#' @seealso \code{\link{wald}}
#' @examples
#' # coming soon FIXME
#' @export
print.wald_lrt <- function(x, round = 6, pround = 5,...) {
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
    
    ta[["Wald p-value"]] <- pformat(ta[["Wald p-value"]])
    if (!is.null(ta[["LRT p-value"]])) {
      ta[["LRT p-value"]] <- pformat(ta[["LRT p-value"]])
    }
    print(as.data.frame(ta, row.names = nn))
    te <- tt$estimate
    te[,'p-value'] <- pformat( te[,'p-value'])
    if ( !is.null(round)) {
      for ( ii in 1:length(te)) {
        te[[ii]] <- rnd(te[[ii]],digits=round)
      }
    }
    temat <- as.matrix(te)
    rownames(temat) <- rownames(tt$L)
    print(temat, quote = F, justify = 'right')
  }
  invisible(x)
}
