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

{
  waldx <- function(fit, Llist = "", clevel = 0.95,
             pred = NULL,
             data = NULL, debug = FALSE , maxrows = 25,
             full = FALSE, fixed = FALSE,
             invert = FALSE, method = 'svd',
             overdispersion = FALSE,
             df = NULL, pars = NULL,...) {
      # New version with support for stanfit
      if (full) return(waldx(fit, getX(fit)))                                     # call waldx
      if(!is.null(pred)) return(waldx(fit, getX(fit,pred)))                       # call waldx
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
      vc <- if(isTRUE(overdispersion)) overdisp_fun(fit)["ratio"] * fix$vcov 
      else if(isFALSE(overdispersion)) fix$vcov
      else fix$vcov
      
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
        
        ## Is L already reduced for a singular model?
        ## 
        Lna <- L[, is.na(beta), drop = FALSE]
        narows <- apply(Lna,1, function(x) sum(abs(x))) > 0
        # give a warning if L has non-zero coefficients for a beta that is NA        
        if(sum(narows) >0) warning("L matrix has non-zero coefficients for estimates that are NA")
        if(sum(is.na(beta)) > 0) {    # singular model
          
          if(ncol(L) != sum(!is.na(beta))) {   # L and non.missing beta of different dimensions
            L <- L[, !is.na(beta),drop = FALSE]                                       # drop NAs from L
            ## restore the data attribute
            attr(L,'data') <- Ldata
            
          }
        }
        isnabeta <- is.na(beta)
        beta <- beta[ !isnabeta ]                                              # always do this
        vc <- vc[!isnabeta, !isnabeta, drop = FALSE]                        # drop NAs from vc
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
        if(!isFALSE(overdispersion)) {
          ret[[ii]]$anova <- c(
            ret[[ii]]$anova, 
            `overdispersion variance factor` =
              if(isTRUE(overdispersion)) overdisp_fun(fit)[["ratio"]]
            else overdispersion)
        }
        if(isTRUE(overdispersion)) {
          ret[[ii]]$anova <- c(
            ret[[ii]]$anova, overdisp_fun(fit)[c(1,3,4)])
        }
        # disp(ret[[ii]]$anova)
        
        ## Estimate
        
        etahat <- L %*% beta
        
        # NAs if not estimable:
        
        # etahat[narows] <- NA
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
}
{
  waldx <- 
    function(fit, Llist = "", clevel = 0.95,
             pred = NULL,
             data = NULL, debug = FALSE , maxrows = 25,
             full = FALSE, fixed = FALSE,
             invert = FALSE, method = 'svd',
             overdispersion = FALSE,
             df = NULL, pars = NULL,...) {
      # New version with support for stanfit
      if (full) return(waldx(fit, getX(fit)))
      if(!is.null(pred)) return(waldx(fit, getX(fit,pred)))
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
      vc <- if(isTRUE(overdispersion)) overdisp_fun(fit)["ratio"] * fix$vcov 
      else if(isFALSE(overdispersion)) fix$vcov
      else fix$vcov
      
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
        
        ## 
        ## 1. Is the model singular? i.e. Are there NA in beta
        ##    If so:
        ##      Determine whether L conforms with non-singular portion of model or full model
        ##      - if full:
        ##        - examine L for non-zero coefficients for NA and report L not estimable
        ##        - truncate L to singular portion
        ##      Truncate beta and vc 
        ## 
        narows <- NULL
        if(sum(is.na(beta)) > 0) {  # singular model
          if(!(ncol(L) %in% c(length(beta), length(na.omit(beta))) )) stop("ncol(L) incorrect for singular or non-singular model")
          if(ncol(L) == length(beta)) {  # L for singular model
            # Is L estimable?
            Lna <- L[, is.na(beta), drop = FALSE]
            narows <- apply(Lna,1, function(x) sum(abs(x))) > 0
            warning("Row(s): ", paste(which(narows),collapse = ' '), " of L not estimable. L coefficients set to 0.")
            L <- L[, !is.na(beta), drop = FALSE]
            attr(L,'data') <- Ldata
          }
          # fix beta and vc
          nas <- is.na(beta)
          beta <- beta[!nas]
          vc <- vc[!nas, !nas, drop= FALSE]
        }


        ## Anova
        if( method == 'qr' ) {
          qqr <- qr(t(na.omit(L))) # omit rows with NAs
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
        if(!isFALSE(overdispersion)) {
          ret[[ii]]$anova <- c(
            ret[[ii]]$anova, 
            `overdispersion variance factor` =
              if(isTRUE(overdispersion)) overdisp_fun(fit)[["ratio"]]
            else overdispersion)
        }
        if(isTRUE(overdispersion)) {
          ret[[ii]]$anova <- c(
            ret[[ii]]$anova, overdisp_fun(fit)[c(1,3,4)])
        }
        # disp(ret[[ii]]$anova)
        
        ## Estimate
        
        etahat <- L %*% beta
        
        # NAs if not estimable:
        if(!is.null(narows)) etahat[narows] <- NA
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
}

undebug(waldx)


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
pred <- merge(pred, list(x=1:10))
ww <- as.data.frame(waldx(fit, pred = pred))
ww$L
waldx(fit, ww$L)

# Workflow for singular models
# 


dd
# differences in each cell from a reference level
# 
pred.grid
pred <- with(dd, pred.grid(a,b, x = 10:12))
pred
waldx(fit, pred = pred)
             