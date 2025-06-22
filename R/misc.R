#' Nickname for getAnywhere
#' 
#' See [utils::getAnywhere()].
#' 
#' @export
ga <- function (x) {
  if (tryCatch(!is.character(x), error = function(e) TRUE)) 
    x <- as.character(substitute(x))
  objs <- list()
  where <- character()
  visible <- logical()
  if (length(pos <- find(x, numeric = TRUE))) {
    objs <- lapply(pos, function(pos, x) get(x, pos = pos), 
                   x = x)
    where <- names(pos)
    visible <- rep.int(TRUE, length(pos))
  }
  if (length(grep(".", x, fixed = TRUE))) {
    np <- length(parts <- strsplit(x, ".", fixed = TRUE)[[1L]])
    for (i in 2:np) {
      gen <- paste(parts[1L:(i - 1)], collapse = ".")
      cl <- paste(parts[i:np], collapse = ".")
      if (gen == "" || cl == "") 
        next
      Call <- substitute(getS3method(gen, cl, TRUE), list(gen = gen, 
                                                          cl = cl))
      f <- eval.parent(Call)
      if (!is.null(f) && !is.null(environment(f))) {
        ev <- topenv(environment(f), baseenv())
        nmev <- if (isNamespace(ev)) 
          getNamespaceName(ev)
        else NULL
        objs <- c(objs, list(f))
        msg <- paste("registered S3 method for", gen)
        if (!is.null(nmev)) 
          msg <- paste(msg, "from namespace", nmev)
        where <- c(where, msg)
        visible <- c(visible, FALSE)
      }
    }
  }
  for (i in loadedNamespaces()) {
    ns <- asNamespace(i)
    if (exists(x, envir = ns, inherits = FALSE)) {
      f <- get(x, envir = ns, inherits = FALSE)
      objs <- c(objs, list(f))
      where <- c(where, paste0("namespace:", i))
      visible <- c(visible, FALSE)
    }
  }
  ln <- length(objs)
  dups <- rep.int(FALSE, ln)
  if (ln > 1L) 
    for (i in 2L:ln) for (j in 1L:(i - 1L)) if (identical(objs[[i]], 
                                                          objs[[j]], ignore.environment = TRUE)) {
      dups[i] <- TRUE
      break
    }
  structure(list(name = x, objs = objs, where = where, visible = visible, 
                 dups = dups), class = "getAnywhere")
}
#' Special contrasts
#' 
#' Normalized Helmert contrasts for parsimonious random effects models
#' 
#' - [R Library Contrast Coding Systems for categorical variables](https://stats.idre.ucla.edu/r/library/r-library-contrast-coding-systems-for-categorical-variables/)
#' 
#' @param n a vector of levels for a factor or the number of levels.
#' @seealso [C_()] for weighted average
#' @examples
#' zf <- factor(c("A","B","B","C"))
#' contrasts(zf) <- contr.nhelmert(3)
#' contrasts(zf)
#' library(nlme)
#' set.seed(123)
#' n_within <- 10
#' n_schools <- 20
#' n_treatments <- 3
#' dd <- expand.grid(id = 1:n_within, schoolid = 1:n_schools)
#' dd <- within(
#'    dd,
#'    {
#'       method <- factor(sample(LETTERS[1:n_treatments], nrow(dd), TRUE))
#'       methodn <- C(method, contr.nhelmert)
#'       methodmat <- contrasts(methodn)[methodn,]
#'       grade <- 10 + seq_len(n_treatments)[method] + rnorm(n_schools)[schoolid] + rnorm(nrow(dd))
#'    })
#' fitfull <- lme(grade ~ method, dd, random = ~ 1 + method | schoolid,
#'  control = list(returnObject = T))     
#' summary(fitfull)
#' getG(fitfull) %>% svd(nu=0,nv=0)
#' fitfulln <- lme(grade ~ method, dd, random = ~ 1 + methodn | schoolid,
#'  control = list(returnObject = T))     
#' summary(fitfulln)
#' getG(fitfulln)
#' getG(fitfulln) %>% svd(nu=0,nv=0)
#' fitpars1 <- lme(grade ~ method, dd, 
#'    random = list( schoolid = pdBlocked(list(pdDiag(~1), pdIdent(~ methodn - 1))) ))
#' fitpars2 <- lme(grade ~ method, dd, 
#'    random = list( schoolid = pdBlocked(list(pdDiag(~1), pdIdent(~ methodmat - 1))) ))
#' summary(fitpars1)
#' summary(fitpars2)
#' ranef(fitpars2) 
#' @export
contr.nhelmert <- function(n, ...) {
  ret <- contr.helmert(n,...)
  ret/sqrt(colSums(ret*ret))
}
#' Contrasts to make intercept a weighted average 
#' 
#' A generalization of treatment contrasts that make
#' the intercept a weighted average of factor levels
#' (instead of placing all the weight on the reference level)
#'  
#' @param x a factor
#' @param weights a vector fo relative weights applied to 
#'        levels of `x`  
#' @seealso [contr.nhelmert()] for weighted average
#' @examples
#' x <- factor(rep(letters[1:4],4))
#' y <- 1:length(x)
#' lm(y ~ x)
#' xx <- C_(x)
#' lm(y ~ xx)
C_ <- function(x, weights = rep(1,length(levels(x)))) {
  x <- as.factor(x)
  if(is.null(dim(weights))){
    nams <- levels(x)
    weights <- rbind(weights, cbind(-1, diag(length(weights)-1)))
    rownames(weights) <- c('all', paste0(' ',nams[-1], '-', nams[1]) )
  }
  C(x, solve(weights)[, -1])
}
#' 
#' Substitute last occurrence of a pattern in a string
#' 
#' @param pattern character string containing a regular expression whose
#'        whose last match in 'x' will be replaced
#' @param replacement character string that is a replacement for the
#'        matched 'pattern' in 'x'
#' @param x a character vector where matches are sought.
#' @param ... arguments passed to \code{\link{gsub}}.
#' 
#' @examples
#' string <- 'a_b_c_D'
#' sublast('_','__',string)
#' dd <- data.frame(id = 1:3, X_a_1 = 1:3, X_a_2 = 1:3, X_b_1 = 1:3, X_b_2 = 1:3)
#' dd
#' names(dd) <- sublast('_','__',names(dd))
#' tolong(dd, sep = '__')
#' tolong(dd, sep = '__') %>% tolong(sep = '_', idvar = 'id2', timevar = 'time2')
#' 
#' @export
sublast <- function(pattern, replacement, x, ...) {
  pat <- paste0('(', pattern, ')(?!.*\\1)')
  # disp(pat)
  sub(pat, replacement, x, ..., perl = TRUE)
}
#' Change NAs to FALSE
#' 
#' @param x vector, possibly with NAs
#' @export
na2f <- function(x) {
  x[is.na(x)] <- FALSE
  x
}
#' Change NAs to TRUE
#' 
#' @param x vector, possibly with NAs
#' @export
na2t <- function(x) {
  x[is.na(x)] <- TRUE
  x
}
# Pipe from magrittr
#
# Removed because of conflict with pipe defined in tidyverse
# @importFrom magrittr %>%
# @export
# magrittr::`%>%`

#' Transform NAs to 0
#'
#' @param x vector, possibly with NAs
#' @export
na20 <- function(x) {
  x[is.na(x)] <- 0
  x
}
#' Directory or filename of active R script
#' 
#' Uses the `this.path` package return the file name
#' of the current R script or its enclosing directory. 
#' 
#' @param dir (default: TRUE) return the directory of the current R script, otherwise the full file name
#' 
#' @export
here <- function(dir = FALSE) {
  if(dir) {
    this.path::this.dir()
  } else {
    this.path::this.path()
  }
}
#' Set working directory to directory of active R script
#' 
#' Set working directory to directory of active R script
#' if script is not being knitted
#' 
#' @export
setwd_here <- function() {
  path <- this.path::this.dir()
  setwd(path)
  invisible(path)
}
#' Vectorized ifelse with multiple conditions
#' 
#' Avoids nested ifelse statements when the action depends on
#' the value of a variable
#' 
#' @param .select. a variable whose values determine the 
#'        argument to be used
#' @param \dots named arguments and one possibly unnamed argument. 
#' 
#' 
#' @details 
#' Each argument in \dots evaluates
#' to a vector whose value is returned where the name of the 
#' argument matches a value of \code{.select.}. 
#' 
#' The vectors in \dots are combined into a matrix with
#' \code{\link{cbind}} and the names of the arguments
#' are used as values of \code{.select.} to select which
#' vector value is returned.  See the examples. 
#' 
#' If there is an unnamed argument, its value is used
#' as a value in \code{.select.} is not matched by
#' an argument name.
#' 
#' See an alternative: \code{\link[dplyr]{case_when}}
#'
#' @examples
#' x <- c(letters[1:4],NA)
#' case(x, a = 'was an a', b = 'was a b', z = 'was a z')
#' case(x, a = 'was an a', x) # x is returned as default
#' # numerical 'select' is coerced to character
#' case(1:4, '1' = 'was a 1', '2' = 'was a 2')
#' 
#' location <- c('England','England','France','France',
#'      'Germany','Spain')
#' xvar <- c('yes','no','non','oui','nein','si')
#' case(location,
#'    'England' = tr(xvar, c('yes','no'), c(1,0)),
#'    'France'  = tr(xvar, c('oui','non'), c(1,0)),
#'    'Germany' = tr(xvar, c('nein','ja'), c(0,1)))
#' case(location,
#'    'England' = tr(xvar, c('yes','no'), c(1,0)),
#'    'France'  = tr(xvar, c('oui','non'), c(1,0)),
#'    'Germany' = tr(xvar, c('nein','ja'), c(0,1)),
#'    xvar) 
#' case(location,
#'    'England' = tr(xvar, c('yes','no'), c(1,0)),
#'    'France'  = tr(xvar, c('oui','non'), c(1,0)),
#'    'Germany' = tr(xvar, c('nein','ja'), c(0,1)),
#'    'no match')
#' @export
case <- function(.select., ...) {
  nas <- is.na(.select.)
  replace <- list(...)
  levels <- names(replace)
  # if "" is in levels, i.e. if there is an unnamed argument
  # then this is the default for non-matches
  # otherwise non-matches return NA
  which <- match(as.character(.select.), levels)
  if(length(default <- grep("^$", levels))) which[is.na(which)] <- default
  # But NAs in select nevertheless return NAs
  which[nas] <- NA
  what <- do.call(cbind, replace)
  what[cbind(1:nrow(what), which)]
}
#' 
#' Sequential ifelse with paired arguments
#' 
#' Equivalent of nested ifelse with alternating conditions and values
#' 
#' @param ... a sequence of alternating arguments with each pair consisting
#'        of a vector logical argument followed by a vector of values to be returned in
#'        the positions in which the logical argument is TRUE. Each pair
#'        corresponds to the first two arguments of a \code{\link{ifelse}}.
#'        To get a default value, use a condition 'TRUE' for the last pair.
#'        You may also need to specify is.na(x) if there are NAs.
#' @return a vector consisting of the value vector corresponding to the first
#'        logical vector that evaluates to TRUE in a position.
#' @export
esac <- function(...) {
# esac <- function(...) {
  # Equivalent of nested ifelse
  # Arguments are alternating (unnamed pairs) giving:
  # condition followed by the value if the condition is satisfied
  # To specify a default: end with: TRUE, default value
  # or use argument OTHER
  a <- list(...)
  sel <- do.call(cbind, a[seq(1,length(a), 2)])
  # repl needs to have the right shape without messing up its type
  repl <- a[seq(2,length(a),2)]
  repl[[1]] <- rep(repl[[1]], length.out = nrow(sel))
  repl <- do.call(cbind, repl)
  # disp(sel)
  first_col <- apply(sel, 1, function(x) min(which(x)))
  repl[cbind(seq_len(nrow(sel)), first_col)]
}
if(FALSE) {
  esac(c(T,F,T), 'a', c(T,F,F), 'b', TRUE, 'c')
}
#'
#' Left Cholesky factor
#' 
#' Decomposes positive-definite G = L'L where L is lower-triangular.
#' 
#' In R, \code{\link{chol}} returns a upper-triangular matrix \code{R}
#' such that G = R'R. \code{lchol} return a lower-triangular matrix.
#' 
#' @param x a positive-definite matrix
#' @examples
#' mm <- cbind( c(8,2,1), c(2,10,2), c(1,2,5))
#' mm
#' chol(mm)
#' lchol(mm)
#' crossprod(chol(mm))
#' t(chol(mm)) %*% chol(mm)
#' crossprod(lchol(mm))
#' t(lchol(mm)) %*% lchol(mm)
#' @export
lchol <- function(x) {
  rind <- rev(1:nrow(x))
  xret <- x[rind,][,rind]
  ret <- chol(xret)
  ret[rind,][,rind]
}
#' Vovk-Sellke Maximum p-Ratio
#' 
#' Calculates the Vovk-Sellke Maximum p-Ratio
#' 
#' @param p p-values
#' @aliases vovk sellke
#' @export 
vs <- function(p) {
  -1/(exp(1) * p * log(p))
}
#' Decomposes positive-definite G = L'L where L is lower-triangular.
#' 
#' In R, \code{\link{chol}} returns a upper-triangular matrix \code{R}
#' such that G = R'R. \code{lchol} return a lower-triangular matrix.
#' 
#' @param x a positive-definite matrix
#' @examples
#' mm <- cbind( c(8,2,1), c(2,10,2), c(1,2,5))
#' mm
#' chol(mm)
#' lchol(mm)
#' crossprod(chol(mm))
#' t(chol(mm)) %*% chol(mm)
#' crossprod(lchol(mm))
#' t(lchol(mm)) %*% lchol(mm)
#' @export
lchol <- function(x) {
  rind <- rev(1:nrow(x))
  xret <- x[rind,][,rind]
  ret <- chol(xret)
  ret[rind,][,rind]
}
#' Cartesian product of variable values for prediction
#'
#' Prediction data frame to caompute predicted values
#' 
#' Plotting predicted values for a model often requires computing
#' predicted values on a grid of predictor values other than the
#' original data set. Categorical variables (character or factor)
#' must be of the same form (same set of character values in the case
#' of character variables and the same levels in the same order in the case of 
#' factor variables) as they appear in the model data frame on which the 
#' model was fitted. \code{pred.grid} facilitates the process, in comparison
#' with \code{\link{expand.grid}}, by allowing references to variables
#' in the model data frame by name or by specifying a vector of values
#' in the case of numeric predictors.
#' 
#' @param ... names of variables or named arguments with values. The 
#'            arguments that are names are evaluated in the environment, which will
#'            typically
#'            be a data frame supplied as an environment with the \code{\link{with}} 
#'            function. Named arguments are evaluted in the usual way.  The unique values 
#'            of each argument are supplied to \code{\link{expand.grid}} to 
#'            create a data frame whose rows are the Cartesian product of the
#'            unique values of the input. See comments on usage in the examples.
#' @examples
#' 
#' hs <- within(
#'    hs,
#'    {
#'       id <- paste(Sector, school) %>% 
#'          as.factor %>% 
#'          reorder(ses + I(Sector == 'Catholic')*1000) 
#'       ses_mean <- capply(ses, id, mean, na.rm = TRUE)  
#'       mathach_mean <- capply(mathach, id, mean, na.rm = TRUE)  
#'    })
#' 
#' fit1 <- lm(mathach ~ (ses + I(ses^2)) * id, hs)
#' pred1 <- with(hs, pred.grid(id, ses = seq(-3,3, .1))) 
#' pred1$fit1 <- predict(fit1, newdata = pred1) 
#' 
#' fit2 <- lm(mathach ~ (ses + I(ses^2))* Sector, hs)
#' # add Sector to pred1
#' head(pred1); dim(pred1)
#' pred2 <- merge(
#'            pred1, 
#'            up(hs, ~id), 
#'            by = 'id', 
#'            all.x = TRUE)
#' head(pred2); dim(pred2)
#' pred2$fit2 <- predict(fit2, newdata = pred2) 
#'
#' # Existing methods allow you to graph lines fitted
#' # within each panel
#' 
#' library(lattice)
#' library(latticeExtra)
#' p <- xyplot(mathach ~ ses | id, hs, groups = Sex, 
#'    layout = c(7,6),
#'    alpha = c(.8, .5),
#'    auto.key = list(space = 'right'),
#'    between = list(y = rep(c(0,.4,0), c(2,1,2))),
#'    skip = rep(c(F,T,F), c(19,2,20)))
#' p
#' p + glayer(panel.smoother(...)) # 
#' 
#' # Using 'pred.grid' and 'expand.grid' makes it easier
#' # to graph lines fitted with models
#' 
#' td(pch = 1, cex = .5)
#' p
#' p + 
#' xyplot(fit1 ~ ses | id, pred1, type = 'l')
#' 
#' p + 
#' xyplot(fit1 ~ ses | id, pred1, type = 'l') +
#' xyplot(fit2 ~ ses | id, pred2, type = 'l', col = 'black') # ?????
#' 
#' pred2 <- sortdf(pred2, ~ ses)
#' p + 
#' xyplot(fit1 ~ ses | id, pred1, type = 'l') +
#' xyplot(fit2 ~ ses | id, pred2, type = 'l', col = 'black') # zig-zagging gone
#'
#' 
#' # Referring to other variables in panel functions
#' 
#' xyplot(mathach ~ ses | id, hs, groups = Sex, 
#'    layout = c(7,6),
#'    cex = .4,
#'    ses_mean = hs$ses_mean,
#'    mathach_mean = hs$mathach_mean,
#'    par.strip.text = list(cex = .7),
#'    auto.key = list(space = 'right'),
#'    between = list(y = rep(c(0,.4,0), c(2,1,2))),
#'    subscripts = TRUE,
#'    skip = rep(c(F,T,F), c(19,2,20))) +
#'  glayer(panel.smoother(..., se = FALSE, lwd = 2, lty =1)) +
#'  layer(panel.abline(v=ses_mean[subscripts],...,col = 'gray')) +  
#'  layer(panel.abline(h=mathach_mean[subscripts],...,col = 'gray'))  
#'  
#' p + glayer(panel.abline(v=ses_mean))
#' p + xyplot(fit1 ~ ses | id, pred1, type = 'l')    
#' @export
pred.grid <-function (...) {
  nams <- as.character(as.list(substitute(list(...)))[-1L])
  x <- list(...)
  if (is.null(names(x))) 
    names(x) <- nams
  else if (any(names(x) == "")) 
    names(x)[names(x) == ""] <- nams[names(x) == ""]
  ret <- lapply(x, unique)
  ret <- lapply(ret, sort)
  ret <- do.call(expand.grid, c(ret, stringsAsFactors = FALSE))
  for(nn in names(ret)){
    if(is.factor(ret[[nn]])) contrasts(ret[[nn]]) <- contrasts(x[[nn]])
  }
  ret
}
#' List debugged functions
#' 
#' @param environments default: search()
#' @param all list all functions including those not debugged, default: FALSE
#' @returns a data frame listing functions and whether they are currently debugged
#' @export
debugged <- function(environments=search(), all = FALSE) {
  r <- do.call("rbind", lapply(environments, function(environment.name) {
    return(do.call("rbind", lapply(ls(environment.name), function(x) {
      if(is.function(get(x))) {
        is.d <- try(isdebugged(get(x)))
        if(!(inherits(is.d,"try-error"))) {
          return(data.frame(function.name=x, debugged=is.d))
        } else { return(NULL) }
      }
    })))
  }))
  if(all) return(r) else subset(r, debugged == TRUE)
} 
#' Load package or install
#' 
#' @param package
#  @examples
#  lib("spida2")
#' @export
# lib <- function(package){
#   package <- as.character(substitute(package))
#   # print(package)
#   if(require(package, character.only = TRUE)) {
#     help(p=package)
#   } else {
#     install.packages(package)
#     library(package)
#     help(p=package)
#   }
# }
