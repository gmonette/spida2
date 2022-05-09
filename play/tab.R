#' Display matrix
#'
#' Transform a matrix of test results for display
#' 
#' @param x matrix
#' @export
.mat2arr <- function(x) {
      ret <- as.list(x)
      dim(ret) <- NULL
      nams <- expand.grid( dimnames(x))
      for ( ii in 1:length(nams)) {
          nams[[ii]] <- paste( names(nams)[ii] , nams[[ii]], sep = " = ")
      }
      nams <- c( nams, sep = ", ")
      nams <- do.call( paste, nams)
      names(ret) <- nams
      ret
}
#' Drop last facets of array
#' 
#' Primarily used to strips totals from a table bordered by totals 
#' by dropping the last facet. 
#'
#' @param arr array
#' @param drop drop parameter in subsetting, default FALSE
#' @param keep names of facets to be kept, default NULL
#' 
#' @examples
#' (arr <- array(1:24,2:4))
#' arr %>% dropLast
#' arr %>% dropLast(drop = TRUE)
#' tab(iris, ~ Species)
#' tab(iris, ~ Species) %>% dropLast 
#' tab(iris, ~ Species + I(Sepal.Length > 5)) %>% dropLast 
#' # row percentages:
#' tab(iris, ~ Species + I(Sepal.Length > 5), pct= 1)
#' tab(iris, ~ Species + I(Sepal.Length > 5), pct= 1) %>% dropLast
#' tab(iris, ~ Species + I(Sepal.Length > 5), pct= 1) %>% dropLast(keep = "All")
#' @export
dropLast <- function(arr, drop = FALSE, keep = NULL) {
  # NEW: 2013-08-20, G. Monette
  # utility function to drop totals from tab
  ind.last <- dim(arr)
  ind.keep <- as.list(- ind.last) # drop last index of each dim unless name in keep
  last <- function(x) x[length(x)]
  names.last <- sapply( dimnames(arr), last)
  disp(names.last)
  if (!is.null(keep)) {
    no.drop <- names.last %in% keep
    ind.keep[no.drop] <- lapply( ind.last[no.drop], seq_len)
  }
  call <- c(list(arr), ind.keep, drop = drop)
  disp(call)
  do.call( `[`,call)
}
#' Drop last elements of an array if it is a "Total"
#'
#' Used to drop "Total" rows and columns after using \code{\link{tab}}. Synonyms for
#' legacy: \code{Tab} and \{pab}.
#'
#' @param mat a matrix, array or table
#' @param names_to_drop (default "Total")
#' @param drop (default FALSE) should one-element dimension be dropped as dimensions
#' @seealso tab, Tab, dropLast
#' @export
dropLastTotal <- function (mat, names_to_drop = "Total", drop = FALSE) {
  cl <- class(mat)
  last <- function(x) x[length(x)]
  cutlast <- function(x) if(length(x)>0) x[-length(x)] else x
  ind.last <- dim(mat)
  ns <- dimnames(mat)
  # disp(ns)
  keep <- lapply(ind.last, seq_len)
  # disp(keep)
  for(i in seq_along(keep)) {
    # disp(i)
    # disp(keep[[i]])
    if( last(ns[[i]]) %in% names_to_drop) keep[[i]] <- cutlast(keep[[i]])
  }
  # disp(keep)
  call <- c(list(mat), keep, drop = drop)
  ret <- do.call(`[`, call)
  class(ret) <- cl
  ret
}

#' Table of frequencies or relative frequencies bordered with totals and
#' including NAs
#'
#' Generates a table of frequencies or relative frequencies or relative
#' percentages
#'
#' @param \dots as with \code{table}, one or more objects which can be interpreted as factors (including character strings), or a list (or data frame) whose components can be so interpreted.
#' @param data a data frame in which formula are interpreted
#' @param fmla a formula whose right-hand side names the variables to be used for tabulation. The optional left-hand side specifies a variable to be used for weights.
#' @param useNA whether to include NA levels. The default is "ifany". Can also
#' be set to "no" or "always".
#' @param pct margins to be scaled to sum to 100. This is the vector of margin indices on which percentages are conditioned. For example, with the call \code{tab(~ A + B + C, data, pct = 2:3)}, the table will contain conditional percentages for variable
#' \code{A} conditional of combinations of \code{B} and \code{C}.
#' @param pr margins to be scaled to sum to 1. This is the vector of margin indices on which percentages are conditioned.
#' @param total.margins if FALSE, generate table without margins
#' @param weights (not working temporarily) instead of generating a frequency table, generate a table
#' with the sum of the weights
#' @param keep names of margins to keep with 'Tab', default = "All". To drop
#' all margins, use default = "".
#' @param test (default FALSE) use \code{\link{chisq.test}}
#' @param simulate (default FALSE) simulate p-value with  \code{\link{chisq.test}}
#' @param B (default 2000) number of replications for simulation
#' @return An object of class 'table' of dimension equal to the number of
#' variables, with optional margins showing totals. Elements of the matrix can
#' be frequencies, relative frequencies, percentages or sums of weights.
#' @aliases tab tab.formula tab.data.frame tab.default Tab
#' @seealso \code{\link{tab_}} to drop "Total" margins and
#' \code{\link{tab__}} to drop "Total" and "All" margins.
#' @author Georges Monette
#' @examples
#' titanic <- as.data.frame(Titanic)
#' head(titanic)
#' tab(titanic, Freq ~ Sex + Survived + Age , test = T)
#' tab(titanic, Freq ~ Sex + Age)
#' tab(titanic, Freq ~ Sex + Survived + Age)
#' round(tab(titanic, Freq ~ Sex + Survived + Age,
#'     pct = c(1,3)),2)
#' round(Tab(titanic, Freq ~ Sex + Survived + Age,
#'     pct = c(1,3)),2)
#' round(Tab(titanic, Freq ~ Sex + Survived + Age,
#'     pct = c(1,3), keep = ""),2)
#' @export
tab <- function(x,...) UseMethod("tab")
#' @describeIn tab method class table
#' @export
tab.table <- function(x,...) {
    data <- as.data.frame(x)
    wt <- data$Freq
    data$Freq <- NULL
    tab.data.frame(data, ..., weights = wt)
}
#' @describeIn tab method for matrices
#' @export
tab.matrix <- function(x,...) tab(as.table(x),...)
#' @describeIn tab method for arrays
#' @export
tab.array <- function(x,...) tab(as.table(x),...)
#' @describeIn tab method for formulas
#' @export
tab.formula <- function( fmla, data = sys.frame(sys.parent()), ... ) tab.data.frame(data,fmla,...)
#' @describeIn tab method for data frames
#' @export
tab.data.frame <-
  function (dd, fmla,
        total.margins = TRUE,
        useNA = "ifany",
        pct = NULL, pr = NULL,
        test = FALSE,
        weights = NULL,
        na.rm = NULL,
        all.label = "All", 
        simulate = FALSE,
        B = 2000)
{
  # GM: 2014 08 22: modified handling of lhs to fix bug when variable in lhs
  #                 also appears in rhs
  if (missing(fmla)) {
    fmla <- parse(text = paste(c("~",paste(names(dd),collapse="+"))))
    fmla <- eval(fmla)
    environment(fmla) <- parent.frame()
  }
  if (is.null(weights) && (length(fmla) >2 )) {
    weights <- model.frame(fmla[-3], dd, na.action = NULL)[[1]]
    xx <- model.frame(fmla[-2], dd, na.action = NULL)
  } else {
    xx = model.frame(fmla, dd, na.action = NULL)
    #weights <- eval(substitute(weights), dd, environment(fmla))  # so weights is evaluated in parent.frame
  }
  if(!is.null(weights) && any(is.na(weights))) {
    warning("NAs in weights set to -Inf")
    weights[ is.na(weights) ] <- -Inf
  }
  #  disp(weights)
  if(simulate ==TRUE) test <- TRUE
  xx = c(xx, list( total.margins = total.margins, useNA = useNA,
                   pct=pct, pr = pr, test=test, na.rm = na.rm,
                   weights = weights,all.label=all.label,
                   simulate=simulate) )
  do.call("tab", xx)
}

#' @describeIn tab default method
#' @export
tab.default <- function (..., total.margins = TRUE,
                         pct = NULL, pr = NULL,
                         useNA = "ifany",
                         test = simulate,
                         weights = NULL,
                         na.rm = NULL,
                         all.label = "All",
                         simulate = FALSE, 
                         B = 2000)
{
  if(!is.null(na.rm)) useNA <- if(na.rm) "no" else "ifany"
  aa <- list(...)
  if (length(aa) == 1 && is.list(aa[[1]])) {
    aa <- c(aa[[1]],list( total.margins = total.margins,
                          useNA = useNA, pr = pr, pct = pct, test=test,
                          na.rm = na.rm, weights = weights))
    disp(aa[[1]])
    return(do.call("tab", aa[[1]]))
  }
  if (is.null(names(aa))) {
    nns = names(match.call())
    names(aa) = nns[2:(1 + length(aa))]
  }
  if(useNA=="ifany") for (ii in 1:length(aa)) aa[[ii]] <- factor(aa[[ii]], exclude = NULL)  # table uses 'ifany' correctly when a number but not for factors
  if( is.null(weights)){
    aa[["useNA"]] <- useNA
    aa[["na.rm"]] <- NULL
    ret <- do.call("table", aa)
  } else {
    ret <- as.table(tapply(weights, aa, sum ) )
    ret[is.na(ret)] <- 0
  }
  if (test) {
    if( length(dim(ret)) < 3) {
      test.out <- chisq.test(ret)
    } else if ( length(dim(ret)) == 3)  {
      test.out <- apply( ret, 3:length(dim(ret)), chisq.test, 
                         simulate.p.value = simulate,
                         B= B)

    } else if ( length(dim(ret)) > 3)  {
      test.out <- .mat2arr(apply( ret, 3:length(dim(ret)), chisq.test, 
                                  simulate.p.value = simulate,
                                  B= B))

    }
  }
  if ( !is.null(pr)) ret <- acond( ret, MARGIN = pr, all.label = all.label)
  else if ( !is.null(pct)) ret <- 100* acond( ret, MARGIN = pct, all.label = all.label)
  else if (total.margins) ret = atotal(ret)
  if( test ) attr(ret,'test') <- test.out
  ret <- as.table(ret)
  if( test ) class(ret) <- c('tab', class(ret))
  ret
}
#' Print method for 'tab' object
#' 
#' @export
print.tab <- function(x,...) {
  NextMethod()
  print(attr(x,'test'))
  invisible(x)
}
#' Transform a frequency table into relative frequencies relative to a margin.
#'
#' @param x a table
#' @param MARGIN on which to condition
#' @param total.margin logical whether to include a total margin (default TRUE)
#' @param all.label label for marginal distribution
#' @export
acond <- function (x, MARGIN = NULL, total.margins = TRUE, all.label = "All")
{
    debug <- F
    x <- if( total.margins ) atotal(x)
    if (is.null(MARGIN)|| max(MARGIN) < 1)
        return(x/x[length(x)])
    d <- dim(x)
    if (debug) {
        cat("dim(x)\n")
        print(d)
    }
    dn <- dimnames(x)
    if (debug) {
        cat("dimnames(x)\n")
        print(dn)
    }
    n <- length(d)
    m <- length(MARGIN)
    if ( length(setdiff( MARGIN, 1:n)) > 0) stop("MARGIN must select dimensions of x, or be equal to 0")
    perm <- c(MARGIN, setdiff(1:n, MARGIN))
    perm.inverse <- order(perm)
    if (debug)
        disp(perm)
    if (debug)
        disp(perm.inverse)
    x <- aperm(x, perm)
    d <- dim(x)
    zl <- list(x)
    for (ii in 1:length(d)) zl[[ii + 1]] <- seq(if (ii <= m)
        1
    else d[ii], d[ii])
    tots <- do.call("[", zl)
    ret <- array(c(x)/c(tots), dim = d)
    ret <- aperm(ret, perm.inverse)
    if (debug)
        disp(dim(ret))
    if (debug)
        disp(dn)
    for ( ii in MARGIN) dn[[ii]][length(dn[[ii]])] <- all.label
    dimnames(ret) <- dn
    ret
}
#' @rdname acond
#' @export
aprop <- acond    # older name


#' Add all marginal sums to an array
#'
#' atotal adds by default a border of sums to an array. The function FUN may be used instead of 'sum'. Additional
#' arguments to FUN can also be given.
#'
#' @param arr array for which margins need to be added
#' @param FUN function to compute margins, default: sum 
#' @param label for margin, default: 'Total'
#' @param \dots additional arguments to 'FUN'
#' @author Georges Monette
#' @return an array with dimension dim(arr) + 1
#' @export
atotal <- function(arr, FUN = sum, label = 'Total', ...) {
	d <- dim(arr)
	cls <- class(arr)
	dim1 <- FALSE  # to handle 1-dimensional tables
	if (length(d) == 1) {
    dim1 <- TRUE
    dn <- dimnames(arr)
		arr <- c(arr)
		d <- dim(arr)
	}
	if(is.character(FUN))
		FUN <- get(FUN, mode = "function")
	else if(mode(FUN) != "function") {
		farg <- substitute(FUN)
		if(mode(farg) == "name")
			FUN <- get(farg, mode = "function")
		else stop(paste("\"", farg, "\" is not a function", sep = ""))
	}

	if (is.null(d)) {
        ret <- structure(c(arr,FUN(arr,...)), names = c(names(arr), label), class = cls)
        if (dim1) {
            dn[[1]] <- c(dn[[1]],label)
            ret <- structure( ret, dim = length(ret), dimnames = dn)
        }
        return (ret)
    }
	n <- length(d)
	ret <- arr
	ind <- 1:n
	for ( i in n:1) {
		new <- apply(ret,ind[-i],FUN,...)
		ret <- abind( ret, new, i, label)
	}
	class( ret ) <- cls
    ret
}

###
###  abind
###
#' Bind comformable arrays
#'
#' \code{abind} binds two conformable arrays along a dimension.
#'
#' dim( arr1 ) and dim( arr2 ) must be equal except in the dth dimension. If
#' the length of dim( arr2 ) is 1 less than that of dim( arr1 ), then 'arr2' is
#' treated as if it had dropped the dth dimension with size 1.
#'
#' If length(dim(arr2)) == length(dim(arr1)) - 1, then arr2 is treated as if it
#' dropped a dimension of size 1 in the dth dimension. 'facename' is used as
#' the name of the dimnames list for this dimension.

#' @param arr1 first array
#' @param arr2 second array
#' @param d dimension along which arr1 and arr2 are joined, or, if the
#' \code{length(dim(arr1)) == length(dim(arr2)) + 1} the dimension in arr1 that
#' is extended by arr2.
#' @param facename Name for the new array dimension
#' @return The returned value, ret, is an array with dimension dim(arr1) except
#' for the dth dimension where dim(ret)[d] == dim(arr1)[d] + dim(arr2)[d].
#'
#' @author Georges Monette
#' @seealso \code{\link[base]{aperm}}, to permute arrays
#' @export
abind <- function(arr1,arr2,d,facename="") {
	d1 <- dim(arr1)
	n1 <- length(d1)
	dn1 <- dimnames(arr1)
	d2 <- dim(arr2)
	n2 <- length(d2)
	dn2 <- dimnames(arr2)

	arenull <- is.null(dn1) & is.null(dn2)
	if (is.null(dn1)){
		dn1 <- lapply( as.list(d1), function(x) seq(1,x))
		dimnames(arr1) <- dn1
	}

	if ( n1 != n2 ) {
		d2 <- d1
		d2[d] <- 1
		dn2 <- dn1
		dn2[[d]] <- facename
		dimnames(arr2) <- NULL
		dim(arr2) <- d2
		dimnames(arr2) <- dn2
		n2 <- n1
	}
	if (is.null(dn2)){
		dn2 <- lapply( as.list(d2), function(x) seq(1,x))
		dimnames(arr2) <- dn2
	}

	# check input for consistency
	# ... later
	#

	perm <- 1:n1
	perm[c(d,n1)] <- c(n1,d)	# perm is an involution

	#
	# permute arr1
	#

	arr.p1 <- aperm(arr1,perm)

	#
	# permute arr2 if length of dimension same as arr1
	#

	arr.p2 <- aperm(arr2,perm)
	dret <- d1[perm]
	dret[n1] <- dret[n1] + d2[d]
	dnret <- dn1
	dnret[[d]] <- c(dnret[[d]],dn2[[d]])

	ret <- c(arr.p1, arr.p2)
	dim(ret) <-  dret

	#
	# permute response back
	#

	ret <- aperm(ret, perm)

	dimnames(ret) <- dnret
	ret
}


#' tab without marginal totals
#'
#' Version of \code{\link{tab}}, with the option to drop selected margins without
#' necessarily dropping
#' the marginal average proportions denoted by "All".
#'
#' @param \dots arguments to the \code{\link{tab}} function.
#' @param names_to_drop (default "Total") names of margins to drop
#' @aliases Tab pab
#' @seealso \code{\link{tab}}
#' @export
tab_ <- function(..., names_to_drop = "Total") {
  # New version of Tab that handles pct and pr
  # To keep the "All" and not the "Total" rows,
  # specify keep = "All"
  # BUGS: would be more efficient if it
  #       called tab(...,total.margins=FALSE)
  #       when pct or pr arguments are not given
  as.table(dropLastTotal(
    tab(..., total.margins = FALSE),
    names_to_drop = names_to_drop))
}
#' tab without marginal totals and 'All'
#'
#' Version of \code{\link{tab}}, with the option to drop selected margins without
#' necessarily dropping
#' the marginal average proportions denoted by "All".
#'
#' @param \dots arguments to the \code{\link{tab}} function.
#' @param names_to_drop (default "Total") names of margins to drop
#' @aliases Tab pab
#' @seealso \code{\link{tab}}
#' @export
tab__ <- function(..., names_to_drop = c("All","Total")) {
  # New version of Tab that handles pct and pr
  # To keep the "All" and not the "Total" rows,
  # specify keep = "All"
  # BUGS: would be more efficient if it
  #       called tab(...,total.margins=FALSE)
  #       when pct or pr arguments are not given
  as.table(dropLastTotal(
    tab(..., total.margins = FALSE),
    names_to_drop = names_to_drop))
}

#' @export
pab <- tab_     # legacy
#' @export
Tab <- tab_    # future?
#' Return a table as as data frame
#' 
#' @param data a data frame
#' @param fmla a one-sided or two-sided formula
#' @param \dots other arguments to \code{\link{tab}}
#' 
#' @details The right hand side of 'fmla' is used for the dimensions of the table. The 
#' left-hand side, if present, is summed to generate the values in the cells of the table. If there is no
#' left-hand side, cell frequencies are shown.   
#' 
#' \code{tab_df(data, fmla)} is similar in effect to \code{as.data.frame(tab_(data, fmla))} except
#' that the former uses the name in the LHS of \code{fmla} as a variable name for counts while the latter
#' uses 'Freq'. Also, \code{tab_df} preserved ordered factors.
#'  
#' @export
tab_df <- function(data, fmla, ...){
  if(missing(fmla)) ret <- as.data.frame( tab_(data, ...))
  else if (length(fmla) >2 ) {
    nam <- as.character(fmla)[2]
    ret <- as.data.frame(tab_(data, fmla, ...))
    ret[[nam]] <- ret$Freq
    ret$Freq <- NULL
  } else ret <- as.data.frame(tab_(data, fmla, ...))
  nams.fac <- intersect(names(data)[sapply(data,is.factor)], names(ret))
  for( nn in nams.fac) {
    if(is.ordered(data[[nn]])) ret[[nn]] <- ordered(ret[[nn]],levels = levels(data[[nn]]))
    else ret[[nn]] <- factor(ret[[nn]], levels = levels(data[[nn]]))
  }
  ret
}