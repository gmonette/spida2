#' @name Deriv
#' @title Symbolic differentiation of an expression or function
#' @description Symbollic differentiation of an expression or function
#' @aliases Deriv drule
#' @concept symbollic differentiation
#'
#' @param f An expression or function to be differentiated.
#'  f can be \itemize{
#'   \item a user defined function: \code{function(x) x**n}
#'   \item a string: \code{"x**n"}
#'   \item an expression: \code{expression(x**n)}
#'   \item a call: \code{call("^", quote(x), quote(n))}
#'   \item a language: \code{quote(x**n)}
#'   \item a right hand side of a formula: \code{~ x**n} or \code{y ~ x**n}
#'  }
#' @param x An optional character vector with variable name(s) with resptect to which
#'  \code{f} must be differentiated. If not provided (i.e. x=NULL), x is
#'  guessed either from\ code{names(formals(f))} (if \code{f} is a function)
#'  or from all variables in f in other cases.
#'  To differentiate expressions including components of lists or vectors, i.e. by expressions like
#'  \code{p[1]}, \code{theta[["alpha"]]} or \code{theta$beta}, the vector of
#'  variables \code{x}
#'  must be a named vector. For the cited examples, \code{x} must be given
#'  as follows \code{c(p="1", theta="alpha", theta="beta")}. Note the repeated name \code{theta} which must be provided for every component of the list \code{theta} by which a
#'  differerentiation is required.
#' @param env An environment where the symbols and functions are searched for.
#'  Defaults to \code{parent.frame()} for \code{f} expression and to
#'  \code{environment(f)} if \code{f} is a function. For primitive function,
#'  it is set by default to .GlobalEnv
#' @param use.D An optional logical (default FALSE), indicates if base::D()
#'  must be used for differentiation of basic expressions.
#' @param cache.exp An optional logical (default TRUE), indicates if
#'  final expression must be optimized with cached subexpressions.
#'  If enabled, repeated calculations are made only once and their
#'  results stored in cache variables which are then reused.
#' @param nderiv An optional integer vector of derivative orders to calculate.
#'  Default NULL value correspond to one differentiation. If length(nderiv)>1,
#'  the resulting expression is a list where each component corresponds to derivative order
#'  given in nderiv. Value 0 corresponds to the original function or expression  non
#'  differentiated. All values must be non negative. If the entries in nderiv
#'  are named, their names are used as names in the returned list. Otherwise
#'  the value of nderiv component is used as a name in the resulting list.
#' @param combine An optional character scalar, it names a function to combine
#'  partial derivatives. Default value is "c" but other functions can be
#'  used, e.g. "cbind" (cf. Details, NB3), "list" or user defined ones. It must
#'  accept any number of arguments or at least the same number of arguments as
#'  there are items in \code{x}.
#'
#' @return \itemize{
#'  \item a function if \code{f} is a function
#'  \item an expression if \code{f} is an expression
#'  \item a character string if \code{f} is a character string
#'  \item a language (usually a so called 'call' but may be also a symbol or just a numeric) for other types of \code{f}
#' }
#'
#' @details
#' R already contains two differentiation functions: D and deriv. D does
#' simple univariate differentiation.  "deriv" uses D to do multivariate
#' differentiation.  The output of "D" is an expression, whereas the output of
#' "deriv" can be an executable function.
#'
#' R's existing functions have several limitations.  They can probably be fixed,
#' but since they are written in C, this would probably require a lot of work.
#' Limitations include:
#' \itemize{
#'  \item The derivatives table can't be modified at runtime, and is only available
#' in C.
#'  \item Function cannot substitute function calls.  eg:
#'	f <- function(x, y) x + y; deriv(~f(x, x^2), "x")
#' }
#'
#' So, here are the advantages of this implementation:
#'
#' \itemize{
#'  \item It is entirely written in R, so would be easier to maintain.
#'  \item Can do multi-variate differentiation.
#'  \item Can differentiate function calls:
#'  \itemize{
#'	   \item if the function is in the derivative table, then the chain rule
#'	is applied.  For example, if you declared that the derivative of
#'	sin is cos, then it would figure out how to call cos correctly.
#'	   \item if the function is not in the derivative table (or it is anonymous),
#'	then the function body is substituted in.
#'	   \item these two methods can be mixed.  An entry in the derivative table
#'	need not be self-contained -- you don't need to provide an infinite
#'	chain of derivatives.
#'  }
#'  \item It's easy to add custom entries to the derivatives table, e.g.
#'
#'   \code{drule[["cos"]] <- alist(x=-sin(x))}
#'
#'   The chain rule will be automatically applied if needed.
#'  \item The output is an executable function, which makes it suitable
#'      for use in optimization problems.
#'  \item Compound functions (i.e. piece-wise functions based on if-else operator) can
#'      be differentiated (cf. examples section).
#'  \item in case of multiple derivatives (e.g. gradient and hessian calculation),
#'      caching can make calculation economies for both
#' }
#'
#' Two work environments \code{drule} and \code{simplifications} are
#' exported in the package namescape.
#' As their names indicate, they contain tables of derivative and
#' simplification rules.
#' To see the list of defined rules do \code{ls(drule)}.
#' To add your own derivative rule for a function called say \code{sinpi(x)} calculating sin(pi*x), do \code{drule[["sinpi"]] <- alist(x=pi*cospi(x))}.
#' Here, "x" stands for the first and unique argument in \code{sinpi()} definition. For a function that might have more than one argument,
#' e.g. \code{log(x, base=exp(1))}, the drule entry must be a list with a named rule
#' per argument. See \code{drule$log} for an example to follow.
#' After adding \code{sinpi} you can differentiate expressions like
#' \code{Deriv(~ sinpi(x^2), "x")}. The chain rule will automatically apply.
#'
#' NB. In \code{abs()} and \code{sign()} function, singularity treatment
#'     at point 0 is left to user's care.
#'     For example, if you need NA at singular points, you can define the following:
#'     \code{drule[["abs"]] <- alist(x=ifelse(x==0, NA, sign(x)))}
#'     \code{drule[["sign"]] <- alist(x=ifelse(x==0, NA, 0))}
#'
#' NB2. In Bessel functions, derivatives are calculated only by the first argument,
#'      not by the \code{nu} argument which is supposed to be constant.
#'
#' NB3. There is a side effect with vector length. E.g. in
#'      \code{Deriv(~a+b*x, c("a", "b"))} the result is \code{c(a = 1, b = x)}.
#'      To avoid the difference in lengths of a and b components (when x is a vector),
#'      one can use an optional parameter \code{combine}
#'      \code{Deriv(~a+b*x, c("a", "b"), combine="cbind")} which gives
#'      \code{cbind(a = 1, b = x)} producing a two column matrix which is
#'      probably the desired result here.
#'      \cr Another example illustrating a side effect is a plain linear
#'      regression case and its Hessian:
#'      \code{Deriv(~sum((a+b*x - y)**2), c("a", "b"), n=c(hessian=2)}
#'      producing just a constant \code{2} for double differentiation by \code{a}
#'      instead of expected result \code{2*length(x)}. It comes from a simplification of
#'      an expression \code{sum(2)} where the constant is not repeated as many times
#'      as length(x) would require it. Here, using the same trick
#'      with \code{combine="cbind"} would not help as all 4 derivatives are just scalars.
#'      Instead, one should modify the previous call to explicitly use a constant vector
#'      of appropriate length:
#'      \code{Deriv(~sum((rep(a, length(x))+b*x - y)**2), c("a", "b"), n=2)}
#'
#' @author Andrew Clausen (original version) and Serguei Sokol (actual version and maintainer)
#' @examples
#'
#' \dontrun{f <- function(x) x^2}
#' \dontrun{Deriv(f)}
#' # function (x)
#' # 2 * x
#'
#' \dontrun{f <- function(x, y) sin(x) * cos(y)}
#' \dontrun{Deriv(f)}
#' # function (x, y)
#' # c(x = cos(x) * cos(y), y = -(sin(x) * sin(y)))
#'
#' \dontrun{f_ <- Deriv(f)}
#' \dontrun{f_(3, 4)}
#' #              x         y
#' # [1,] 0.6471023 0.1068000
#'
#' \dontrun{Deriv(~ f(x, y^2), "y")}
#' # -(2 * (y * sin(x) * sin(y^2)))
#'
#' \dontrun{Deriv(quote(f(x, y^2)), c("x", "y"), cache.exp=FALSE)}
#' # c(x = cos(x) * cos(y^2), y = -(2 * (y * sin(x) * sin(y^2))))
#'
#' \dontrun{Deriv(expression(sin(x^2) * y), "x")}
#' # expression(2*(x*y*cos(x^2)))
#'
#' Deriv("sin(x^2) * y", "x") # differentiate only by x
#' "2 * (x * y * cos(x^2))"
#'
#' Deriv("sin(x^2) * y", cache.exp=FALSE) # differentiate by all variables (here by x and y)
#' "c(x = 2 * (x * y * cos(x^2)), y = sin(x^2))"
#'
#' # Compound function example (here abs(x) smoothed near 0)
#' fc <- function(x, h=0.1) if (abs(x) < h) 0.5*h*(x/h)**2 else abs(x)-0.5*h
#' Deriv("fc(x)", "x", cache.exp=FALSE)
#' "if (abs(x) < h) x/h else sign(x)"
#'
#' # Example of a first argument that cannot be evaluated in the current environment:
#' \dontrun{
#'   suppressWarnings(rm("xx", "yy"))
#'   Deriv(xx^2+yy^2)
#' }
#' # c(xx = 2 * xx, yy = 2 * yy)
#'
#' # Automatic differentiation (AD), note itermediate variable 'd' assignment
#' \dontrun{Deriv(~{d <- ((x-m)/s)^2; exp(-0.5*d)}, "x")}
#' #{
#' #   d <- ((x - m)/s)^2
#' #   .d_x <- 2 * ((x - m)/s^2)
#' #   -(0.5 * (.d_x * exp(-(0.5 * d))))
#' #}
#'
#' # Custom derivation rule
#' \dontrun{
#'   myfun <- function(x, y=TRUE) NULL # do something usefull
#'   dmyfun <- function(x, y=TRUE) NULL # myfun derivative by x.
#'   drule[["myfun"]] <- alist(x=dmyfun(x, y), y=NULL) # y is just a logical
#'   Deriv(myfun(z^2, FALSE), "z")
#'   # 2 * (z * dmyfun(z^2, FALSE))
#' }
#' # Differentiantion by list components
#' \dontrun{
#'   theta <- list(m=0.1, sd=2.)
#'   x <- names(theta)
#'   names(x)=rep("theta", length(theta))
#'   Deriv(~exp(-(x-theta$m)**2/(2*theta$sd)), x, cache.exp=FALSE)
#' # c(theta_m = exp(-((x - theta$m)^2/(2 * theta$sd))) *
#' #  (x - theta$m)/theta$sd, theta_sd = 2 * (exp(-((x - theta$m)^2/
#' #  (2 * theta$sd))) * (x - theta$m)^2/(2 * theta$sd)^2))
#' }

Deriv <- function(f, x=if (is.function(f)) NULL else all.vars(if (is.character(f)) parse(text=f) else f), env=if (is.function(f)) environment(f) else parent.frame(), use.D=FALSE, cache.exp=TRUE, nderiv=NULL, combine="c") {
	tf <- try(f, silent=TRUE)
	fch <- deparse(substitute(f))
	if (is.primitive(f)) {
		# get the true function name (may be after renaming in caller env like f=cos)
		fch=sub('^\\.Primitive\\("(.+)"\\)', "\\1", format1(f))
	}
	if (inherits(tf, "try-error")) {
		f <- substitute(f)
	}
	# create dsym and scache local envs (to keep clean nested calls)
	dsym <- new.env()
	dsym$l <- list()
	scache <- new.env()
	scache$l <- list()

	if (is.null(env))
		env <- .GlobalEnv
	if (is.null(x)) {
		# primitive function or function given by a list membr or alike
		af <- formals(args(f))
		x <- names(af)
		rule <- drule[[fch]]
		if (!is.null(rule)) {
			# exclude arguments by which we cannot not differentiate from x
			x=as.list(x)
			x[sapply(rule, is.null)] <- NULL
			if (length(x) == 0) {
				stop(sprintf("There is no differentiable argument in the function %s", fch))
			}
			x=unlist(x)
		}
		if (is.function(f)) {
			rmget=mget(fch, mode="function", envir=env, inherits=TRUE, ifnotfound=NA)
			if (!is.function(rmget[[fch]])) {
				# no function with name stored in fch => replace it by its body
				fd=body(f)
			} else {
				fd <- as.call(c(as.symbol(fch), lapply(names(af), as.symbol)))
			}
		}
		pack_res <- as.call(alist(as.function, c(af, list(res)), envir=env))
	} else {
		x[] <- as.character(x)
		if (any(nchar(x) == 0)) {
			stop("Names in the second argument must not be empty")
		}
		fd <- NULL
	}
	# prepare fd (a call to differentiate)
	# and pack_res (a cal to evaluate and return as result)
	if (!is.null(fd)) {
		; # we are already set
	} else if (is.character(f)) {
		# f is to parse
		fd <- parse(text=f)[[1]]
		pack_res <- as.call(alist(format1, res))
	} else if (is.function(f)) {
#browser()
		b <- body(f)
		if ((is.call(b) && (b[[1]] == as.symbol(".Internal") || b[[1]] == as.symbol(".External") || b[[1]] == as.symbol(".Call"))) || (is.null(b) && (is.primitive(f)) || !is.null(drule[[fch]]))) {
			if (fch %in% dlin || !is.null(drule[[fch]])) {
				arg <- lapply(names(formals(args(f))), as.symbol)
				fd <- as.call(c(as.symbol(fch), arg))
				pack_res <- as.call(alist(as.function, c(formals(args(f)), list(res)), envir=env))
			} else {
				stop(sprintf("Internal or external function '%s()' is not in derivative table.", fch))
			}
		} else {
			fd <- b
			pack_res <- as.call(alist(as.function, c(formals(f), list(res)), envir=env))
		}
	} else if (is.expression(f)) {
		fd <- f[[1]]
		pack_res <- as.call(alist(as.expression, res))
	} else if (is.language(f)) {
		if (is.call(f) && f[[1]] == as.symbol("~")) {
			# rhs of the formula
			fd <- f[[length(f)]]
			pack_res <- quote(res)
		} else {
			# plain call derivation
			fd <- f
			pack_res <- quote(res)
		}
	} else {
		fd <- substitute(f)
		pack_res <- quote(res)
		#stop("Invalid type of 'f' for differentiation")
	}
	res <- Deriv_(fd, x, env, use.D, dsym, scache, combine)
	if (!is.null(nderiv)) {
		# multiple derivatives
		# prepare their names
		if (any(nderiv < 0)) {
			stop("All entries in nderiv must be non negative")
		}
		nm_deriv <- names(nderiv)
		nderiv <- as.integer(nderiv)
		if (is.null(nm_deriv))
			nm_deriv <- nderiv
		iempt <- nchar(nm_deriv)==0
		nm_deriv[iempt] <- seq_along(nderiv)[iempt]
		# prepare list of repeated derivatives
		lrep <- as.list(nderiv)
		names(lrep) <- nm_deriv
		# check if 0 is nderiv
		iz <- nderiv==0
		lrep[iz] <- list(fd)
		# set first derivative
		i <- nderiv==1
		lrep[i] <- list(res)
		maxd <- max(nderiv)
		for (ider in seq_len(maxd)) {
			if (ider < 2)
				next
			res <- Deriv_(res, x, env, use.D, dsym, scache, combine)
			i <- ider == nderiv
			lrep[i] <- list(res)
		}
		if (length(lrep) == 1) {
			res <- lrep[[1]]
		} else {
			res <- as.call(c(quote(list), lrep))
		}
	}
#browser()
	# This subs in sub expressions
	#print(res)
	#if (cache.exp)
	#	res <- Cache(Simplify(deCache(res), scache=scache))
	#print(res)
	eval(pack_res)
}

# workhorse function doing the main work of differentiation
Deriv_ <- function(st, x, env, use.D, dsym, scache, combine="c") {
	stch <- format1(if (is.call(st)) st[[1]] else st)
	# Make x scalar and wrap results in a c() call if length(x) > 1
	iel=which("..." == x)
	if (length(iel) > 0) {
		# remove '...' from derivable arguments
		x=as.list(x)
		x[iel]=NULL
		x=unlist(x)
	}
	nm_x <- names(x)
	if (!is.null(nm_x))
		nm_x[is.na(nm_x)] <- ""
	else
		nm_x <- rep("", length(x))
	if (length(x) > 1 && stch != "{") {
#browser()
		# many variables => recursive call on single name
		# we exclude the case '{' as we put partial derivs inside of '{'
		# so it can be well optimized by Cache()
		res <- lapply(seq_along(x), function(ix) Deriv_(st, x[ix], env, use.D, dsym, scache, combine))
		names(res) <- if (is.null(nm_x)) x else ifelse(is.na(nm_x) | nchar(nm_x) == 0, x, paste(nm_x, x, sep="_"));
		return(as.call(c(as.symbol(combine), res)))
	}
	# differentiate R statement 'st' (a call, or a symbol or numeric) by a name in 'x'
	get_sub_x <- !(is.null(nm_x) | nchar(nm_x) == 0 | is.na(nm_x))
	is_index_expr <- is.call(st) && any(format1(st[[1]]) == c("$", "[", "[["))
	is_sub_x <- is_index_expr &&
				format1(st[[2]]) == nm_x && format1(st[[3]]) == x
	if (is.conuloch(st) || (is_index_expr && !is_sub_x)) {
		return(0)
	} else if (is.symbol(st) || (get_sub_x && is_index_expr)) {
#browser()
		stch <- format1(st)
		if ((stch == x && !get_sub_x) || (get_sub_x && is_sub_x)) {
			return(1)
		} else if ((get_sub_x && is_index_expr && !is_sub_x) ||
				(if (get_sub_x) is.null(dsym$l[[nm_x]][[x]][[stch]]) else
				is.null(dsym$l[[x]][[stch]]))) {
			return(0)
		} else {
			return(if (get_sub_x) dsym$l[[nm_x]][[x]][[stch]] else dsym$l[[x]][[stch]])
		}
	} else if (is.call(st)) {
#browser()
		stch <- format1(st[[1]])
		args <- as.list(st)[-1]
		if (stch %in% dlin) {
			# linear case
			# differentiate all arguments then pass them to the function
			dargs <- lapply(args, Deriv_, x, env, use.D, dsym, scache)
			return(Simplify_(as.call(c(st[[1]], dargs)), scache))
		}
		nb_args=length(st)-1
		# special cases: out of rule table or args(stch) -> NULL
		if (stch == "{") {
#browser()
			# AD differentiation (may be with many x)
			res <- list(st[[1]])
			# initiate dsym[[x[ix]]] or dsym[[nm_x[ix]}}[[x[ix]]]
			for (ix in seq_along(x)) {
				if (get_sub_x[ix]) {
					if (is.null(dsym$l[[nm_x[ix]]][[x[ix]]]))
						dsym$l[[nm_x[ix]]][[x[ix]]] <- list()
				} else {
					if (is.null(dsym$l[[x[ix]]]))
						dsym$l[[x[ix]]] <- list()
				}
			}
			# collect defined var names (to avoid redifferentiation)
			defs <- sapply(args, function(e) if (is.assign(e)) format1(e[[2]]) else "")
#			alva <- list()
			last_res <- list()
			for (iarg in seq_along(args)) {
#browser()
				a <- args[[iarg]]
				if (is.assign(a)) {
					if (!is.symbol(a[[2]]))
						stop(sprintf("In AD mode, don't know how to deal with a non symbol '%s' at lhs", format1(a[[2]])))
					# put in scache the assignement
					ach <- format1(a[[2]])
					for (ix in seq_along(x)) {
						d_ach <- paste(".", ach, "_", x[ix], sep="")
						d_a <- as.symbol(d_ach)
						if (any(d_ach == defs)) {
							# already differentiated in previous calls
							if (get_sub_x[ix])
								dsym$l[[nm_x[ix]]][[x[ix]]][[ach]] <- d_a
							else
								dsym$l[[x[ix]]][[ach]] <- d_a
							next
						}
						de_a <- Deriv_(a[[3]], x[ix], env, use.D, dsym, scache)
						if (get_sub_x[ix])
							dsym$l[[nm_x[ix]]][[x[ix]]][[ach]] <- de_a
						else
							dsym$l[[x[ix]]][[ach]] <- de_a
						if (is.numconst(de_a, 0)) {
							if (iarg < length(args))
								next
						} else if (!is.call(de_a)) {
							if (iarg < length(args))
								next
						}
						if (get_sub_x[ix])
							dsym$l[[nm_x[ix]]][[x[ix]]][[ach]] <- d_a
						else
							dsym$l[[x[ix]]][[ach]] <- d_a
						res <- append(res, call("<-", d_a, de_a))
#						alva <- append(alva, list(c(d_ach, all.vars(de_a))))
						# store it in scache too
						#scache$l[[format1(de_a)]] <- as.symbol(d_a)
						if (iarg == length(args))
							last_res[[ix]] <- d_a
					}
					Simplify_(a, scache)
					res <- append(res, a)
#					alva <- append(alva, list(all.vars(a)))
					if (iarg == length(args)) {
						names(last_res) <- ifelse(get_sub_x, paste(nm_x, x, sep="_"), x)
						res <- append(res, as.call(c(as.symbol(combine), last_res)))
					}
				} else {
					de_a <- lapply(seq_along(x), function(ix) Deriv_(a, x[ix], env, use.D, dsym, scache))
					if (length(x) > 1) {
						names(de_a) <- ifelse(get_sub_x, paste(nm_x, x, sep="_"), x)
						res <- append(res, as.call(c(as.symbol(combine), de_a)))
					} else {
						res <- append(res, de_a)
					}
				}
			}
#browser()
#			if (length(alva) == length(res)) {
#				i <- toporder(alva[-length(alva)]) # the last expression must stay the last
#			} else {
#				i <- toporder(alva)
#			}
#			res[-c(1, length(res))] <- res[-c(1, length(res))][i]
			return(Simplify(as.call(res)))
		} else if (is.uminus(st)) {
			return(Simplify(call("-", Deriv_(st[[2]], x, env, use.D, dsym, scache)), scache=scache))
		} else if (stch == "(") {
#browser()
			return(Simplify(Deriv_(st[[2]], x, env, use.D, dsym, scache), scache=scache))
		} else if(stch == "ifelse") {
			return(Simplify(call("ifelse", st[[2]], Deriv_(st[[3]], x, env, use.D, dsym, scache),
				Deriv_(st[[4]], x, env, use.D, dsym, scache)), scache=scache))
		} else if (stch == "rep") {
#browser()
			# 'x' argument is named or positional?
			i=if ("x" %in% names(st)) "x" else 2
			dst=st
			dst[[i]]=Simplify(Deriv_(st[[i]], x, env, use.D, dsym, scache), scache=scache)
			return(dst)
		} else if (stch == "if") {
			return(if (nb_args == 2)
				Simplify(call("if", st[[2]], Deriv_(st[[3]], x, env, use.D, dsym, scache)), scache=scache)
				else
				Simplify(call("if", st[[2]], Deriv_(st[[3]], x, env, use.D, dsym, scache),
					Deriv_(st[[4]], x, env, use.D, dsym, scache)), scache=scache))
		}
		rule <- drule[[stch]]
		if (is.null(rule)) {
#browser()
			# no derivative rule for this function
			# see if its arguments depend on x. If not, just send 0
			dargs <- lapply(args, Deriv_, x, env, use.D, dsym, scache)
			if (all(sapply(dargs, identical, 0))) {
				return(0)
			}
			else if(stch == "M" & grepl(x,args))
			{
			  dargs <- call("M",dargs[[1]])
			  return(dargs)
			}
			else if(stch == "I" & grepl(x,args))
			{
			  dargs <- call("I",dargs[[1]])
			  return(dargs)
			}
			else
			{
  			# otherwise try to get the body and differentiate it
  			ff <- get(stch, mode="function", envir=env)
  			
  			bf <- body(ff)
  			if (is.null(bf)) {
  				stop(sprintf("Could not retrieve body of '%s()'", stch))
  			}
  			if (is.call(bf) && (bf[[1]] == as.symbol(".External") || bf[[1]] == as.symbol(".Internal") || bf[[1]] == as.symbol(".Call"))) {
  #cat("aha\n")
  				stop(sprintf("Function '%s()' is not in derivative table", stch))
  			}
  			mc <- match.call(ff, st)
  			st <- Simplify_(do.call("substitute", list(bf, as.list(mc)[-1])), scache)
  			return(Deriv_(st, x, env, use.D, dsym, scache))
			}
		}
		# there is a rule!
		if (use.D) {
			return(Simplify(D(st, x), scache=scache))
		}
#if (stch == "myfun")
#browser()
		# prepare replacement list
		da <- try(args(stch), silent=TRUE)
		if (inherits(da, "try-error")) {
			# last chance to get unknown function definition
			# may be it is a user defined one?
			da <- args(get(stch, mode="function", envir=env))
		}
		mc <- as.list(match.call(definition=da, call=st, expand.dots=FALSE))[-1]
		da <- as.list(da)
		da <- da[-length(da)] # all declared arguments with default values
		aa <- modifyList(da, mc) # all arguments with actual values
		# actualize the rule with actual arguments
		rule <- lapply(rule, function(r) do.call("substitute", list(r, aa)))
#browser()
		# which arguments can be differentiated?
		iad <- which(!sapply(rule, is.null))
		rule <- rule[iad]
		lsy <- unlist(lapply(dsym$l, function(it) if (get_sub_x && is.list(it)) unlist(lapply(it, ls, all.names=TRUE)) else ls(it, all.names=TRUE)))
		if (!any(names(which(sapply(mc, function(it) {av <- all.vars(it); (if (get_sub_x) any(nm_x == av) else any(x == av)) || any(av %in% lsy)}))) %in% names(rule))) {
			#warning(sprintf("A call %s cannot be differentiated by the argument '%s'", format1(st), x))
			return(0)
		}
		n <- length(rule)
		
		dargs <- lapply(names(rule), function(nm_a) if (is.null(mc[[nm_a]])) 0 else Deriv_(mc[[nm_a]], x, env, use.D, dsym, scache))
		ize <- sapply(dargs, identical, 0)
		dargs <- dargs[!ize]
		rule <- rule[!ize]
		if (length(rule) == 0) {
			return(0)
		}

		# apply chain rule where needed
		ione <- sapply(dargs, identical, 1)
		imone <- sapply(dargs, identical, -1)
		#print(rule)
		for (i in seq_along(rule)[!(ione|imone)]) {
			rule[[i]] <- Simplify(call("*", dargs[[i]], rule[[i]]), scache=scache)
		}
		for (i in seq_along(rule)[imone]) {
			rule[[i]] <- Simplify(call("-", rule[[i]]), scache=scache)
		}
		#print(rule)
		return(Simplify(li2sum(rule), scache=scache))
	} else if (is.function(st)) {
#browser()
		# differentiate its body if can get it
		args <- as.list(st)[-1]
		names(args)=names(formals(ff))
		if (is.null(names(args))) {
			stop(sprintf("Could not retrieve arguments of '%s()'", stch))
		}
		st <- do.call("substitute", list(body(ff), args))
		Deriv_(st, x, env, use.D, dsym, scache)
	} else {
		stop("Invalid type of 'st' argument. It must be constant, symbol or a call.")
	}
}

drule <- new.env()

# linear functions, i.e. d(f(x))/dx == f(d(arg)/dx)
dlin=c("+", "-", "c", "t", "sum", "cbind", "rbind", "list")

# rule table
# arithmetics
drule[["*"]] <- alist(e1=e2, e2=e1)
drule[["^"]] <- alist(e1=e2*e1^(e2-1), e2=e1^e2*log(e1))
drule[["/"]] <- alist(e1=1/e2, e2=-e1/e2^2)
# log, exp, sqrt
drule[["sqrt"]] <- alist(x=0.5/sqrt(x))
drule[["log"]] <- alist(x=1/(x*log(base)), base=-log(x, base)/(base*log(base)))
drule[["logb"]] <- drule[["log"]]
drule[["log2"]] <- alist(x=1/(x*log(2)))
drule[["log10"]] <- alist(x=1/(x*log(10)))
drule[["log1p"]] <- alist(x=1/(x+1))
drule[["exp"]] <- alist(x=exp(x))
drule[["expm1"]] <- alist(x=exp(x))
# trigonometric
drule[["sin"]] <- alist(x=cos(x))
drule[["cos"]] <- alist(x=-sin(x))
drule[["tan"]] <- alist(x=1/cos(x)^2)
drule[["asin"]] <- alist(x=1/sqrt(1-x^2))
drule[["acos"]] <- alist(x=-1/sqrt(1-x^2))
drule[["atan"]] <- alist(x=1/(1+x^2))
drule[["atan2"]] <- alist(y=x/(x^2+y^2), x=-y/(x^2+y^2))
if (getRversion() >= "3.1.0") {
	drule[["sinpi"]] <- alist(x=pi*cospi(x))
	drule[["cospi"]] <- alist(x=-pi*sinpi(x))
	drule[["tanpi"]] <- alist(x=pi/cospi(x)^2)
}
# hyperbolic
drule[["sinh"]] <- alist(x=cosh(x))
drule[["cosh"]] <- alist(x=sinh(x))
drule[["tanh"]] <- alist(x=(1-tanh(x)^2))
drule[["asinh"]] <- alist(x=1/sqrt(x^2+1))
drule[["acosh"]] <- alist(x=1/sqrt(x^2-1))
drule[["atanh"]] <- alist(x=1/(1-x^2))
# sign depending functions
drule[["abs"]] <- alist(x=sign(x))
drule[["sign"]] <- alist(x=0)
#drule[["abs"]] <- alist(x=ifelse(x==0, NA, sign(x)))
#drule[["sign"]] <- alist(x=ifelse(x==0, NA, 0))
# special functions
drule[["besselI"]] <- alist(x=(if (nu == 0) besselI(x, 1, expon.scaled) else 0.5*(besselI(x, nu-1, expon.scaled) + besselI(x, nu+1, expon.scaled)))-if (expon.scaled) besselI(x, nu, TRUE) else 0, nu=NULL, expon.scaled=NULL)
drule[["besselK"]] <- alist(x=(if (nu == 0) -besselK(x, 1, expon.scaled) else -0.5*(besselK(x, nu-1, expon.scaled) + besselK(x, nu+1, expon.scaled)))+if (expon.scaled) besselK(x, nu, TRUE) else 0, nu=NULL, expon.scaled=NULL)
drule[["besselJ"]] <- alist(x=if (nu == 0) -besselJ(x, 1) else 0.5*(besselJ(x, nu-1) - besselJ(x, nu+1)), nu=NULL)
drule[["besselY"]] <- alist(x=if (nu == 0) -besselY(x, 1) else 0.5*(besselY(x, nu-1) - besselY(x, nu+1)), nu=NULL)
drule[["gamma"]] <- alist(x=gamma(x)*digamma(x))
drule[["lgamma"]] <- alist(x=digamma(x))
drule[["digamma"]] <- alist(x=trigamma(x))
drule[["trigamma"]] <- alist(x=psigamma(x, 2L))
drule[["psigamma"]] <- alist(x=psigamma(x, deriv+1L), deriv=NULL)
drule[["beta"]] <- alist(a=beta(a, b)*(digamma(a)-digamma(a+b)), b=beta(a, b)*(digamma(b)-digamma(a+b)))
drule[["lbeta"]] <- alist(a=digamma(a)-digamma(a+b), b=digamma(b)-digamma(a+b))
# probability densities
drule[["dbinom"]] <- alist(x=NULL, size=NULL, prob=if (size == 0) -x*(1-prob)^(x-1) else if (x == size) size*prob^(size-1) else (size-x*prob)*(x-size+1)*dbinom(x, size-1, prob)/(1-prob)^2/(if (log) dbinom(x, size, prob) else 1), log=NULL)
drule[["dnorm"]] <- alist(x=-(x-mean)/sd^2*if (log) 1 else dnorm(x, mean, sd),
	mean=(x-mean)/sd^2*if (log) 1 else dnorm(x, mean, sd),
	sd=(((x - mean)/sd)^2 - 1)/sd * if (log) 1 else dnorm(x, mean, sd),
	log=NULL)
drule[["pnorm"]] <- alist(q=dnorm(q, mean, sd)*(if (lower.tail) 1 else -1)/(if (log.p) pnorm(q, mean, sd, lower.tail) else 1), mean=dnorm(q, mean, sd)*(if (lower.tail) -1 else 1)/(if (log.p) pnorm(q, mean, sd, lower.tail) else 1), sd=dnorm(q, mean, sd)*(mean-q)/sd*(if (lower.tail) 1 else -1)/(if (log.p) pnorm(q, mean, sd, lower.tail) else 1), lower.tail=NULL, log.p=NULL)
# data mangling
#drule[["rep"]] <- alist(x=rep(1, ...)) # cannot handle '...' yet
drule[["rep.int"]] <- alist(x=rep.int(1, times), times=NULL)
drule[["rep_len"]] <- alist(x=rep_len(1, length.out), length.out=NULL)

#' @name Simplify
#' @title Symbollic simplification of an expression or function
#' @description Symbollic simplification of an expression or function
#' @aliases Simplify simplifications Cache deCache
#' @concept symbolic simplification
# \usage{
# Simplify(expr, env=parent.frame(), scache=new.env())
# }
#' 
#' 
#' @param expr An expression to be simplified, expr can be
#' \itemize{
#'    \item an expression: \code{expression(x+x)}
#'    \item a string: \code{"x+x"}
#'    \item a function: \code{function(x) x+x}
#'    \item a right hand side of a formula: \code{~x+x}
#'    \item a language: \code{quote(x+x)}
#' }
#' @param env An environment in which a simplified function is created
#'  if \code{expr} is a function. This argument is ignored in all other cases.
#' @param scache An environment where there is a list in which simplified expression are cached
#' @param st A language expression to be cached
#' @param prefix A string to start the names of the cache variables
#' @return A simplified expression. The result is of the same type as
#'  \code{expr} except for formula, where a language is returned.
#' @details An environment \code{simplifications} containing simplification rules, is exported in the namespace accessible by the user.
#'  Cache() is used to remove redundunt calculations by storing them in
#'  cache variables. Default parameters to Cache() does not have to be provided
#'  by user. deCache() makes the inverse job -- a series of assignements
#'  are replaced by only one big expression without assignement.
#'  Sometimes it is usefull to
#'  apply deChache() and only then pass its result to Cache().
Simplify <- function(expr, env=parent.frame(), scache=new.env()) {
  if (is.null(scache$l))
    scache$l <- list() # for stand alone use of Simplify
  if (is.expression(expr)) {
    as.expression(Simplify_(expr[[1]], scache))
  } else if (is.function(expr)) {
    as.function(c(as.list(formals(expr)),
                  Simplify_(body(expr), scache)),
                envir=env)
  } else if (is.call(expr) && expr[[1]] == as.symbol("~")) {
    Simplify_(expr[[length(expr)]], scache)
  } else if (is.character(expr)) {
    format1(Simplify_(parse(text=expr)[[1]], scache))
  } else {
    Simplify_(expr, scache)
  }
}

#' @name format1
#' @title Wrapper for base::format() function
#' @description Wrapper for base::format() function
# \usage{
# format1(expr)
# }
#' 
#' 
#' @param expr An expression or symbol or language to be converted to a string.
#' @return A character vector of length 1 contrary to base::format() which
#'  can split its output over several lines.
format1 <- function(expr) {
  res <- if (is.symbol(expr)) as.character(expr) else if (is.call(expr) && expr[[1]]==as.symbol("{")) sapply(as.list(expr), format1) else format(expr)
  n <- length(res)
  if (n > 1) {
    if (res[1] == "{" && n > 2) {
      b <- paste0(res[-1], collapse="; ")
      res <- paste0("{", b, "}", collapse="")
    } else {
      res <- paste0(res, collapse=" ")
    }
  }
  return(res)
}

Simplify_ <- function(expr, scache) {
  if (is.call(expr)) {
    che <- format1(expr)
    res <- scache$l[[che]]
    if (!is.null(res)) {
      if (typeof(res) == "logical" && is.na(res)) {
        # recursive infinite call
        scache$l[[che]] <- expr
        return(expr)
      } else {
        return(res)
      }
    }
    scache$l[[che]] <- NA # token holder
    #cat("simp expr=", format1(expr), "\n", sep="")
    args <- lapply(as.list(expr)[-1], Simplify_, scache)
    expr[-1]=args
    if (all(sapply(args, is.conuloch))) {
      # if all arguments are like numeric, evaluate them
      res <- eval(expr)
      scache$l[[che]] <- res
      return(res)
    } else {
      # is there a rule in the table?
      sym.name <- format1(expr[[1]])
      Simplify.rule <- simplifications[[sym.name]]
      res <- if (!is.null(Simplify.rule)) Simplify.rule(expr, scache=scache) else expr
      scache$l[[che]] <- res
      return(res)
    }
  } else {
    expr
  }
}

# in what follows no need to Simplify_ args neither to check if
# all arguments are numeric. It is done in the upper Simplify_()
`Simplify.(` <- function(expr, scache=NULL) {
  expr[[2]]
}
`Simplify.+` <- function(expr, add=TRUE, scache=NULL) {
  if (length(expr) == 2) {
    if (add)
      return(expr[[2]])
    else if (is.uminus(expr[[2]]))
      return(expr[[2]][[2]])
    else if (is.uplus(expr[[2]]))
      return(call("-", expr[[2]][[2]]))
    else
      return(expr)
  }
  a <- expr[[2]]
  b <- expr[[3]]
  
  if (is.numconst(a, 0) || (is.call(a) && format1(a[[1]]) %in% c("rep", "rep.int", "rep_len") && is.numconst(a[[2]], 0))) {
    #browser()
    return(if (add) b else call("-", b))
  } else if (is.numconst(b, 0) || (is.call(b) && format1(b[[1]]) %in% c("rep", "rep.int", "rep_len") && is.numconst(b[[2]], 0))) {
    #browser()
    return(a)
  } else if (add && is.uminus(a) && !is.uminus(b)) {
    a <- b
    b <- expr[[2]][[2]]
    add <- FALSE
    expr <- call("-", a, b)
  } else if (identical(a, b)) {
    return(if (add) Simplify_(call("*", 2, a), scache) else 0)
  } else if (!is.call(a) && !is.call(b)) {
    if (add) {
      # just reorder
      expr[-1] <- expr[1+order(sapply(expr[-1], format1))]
    }
    return(expr)
  }
  # factorise most repeated terms
  alc <- Lincomb(a)
  blc <- Lincomb(b)
  if (add) {
    lc <- c(alc, blc)
  } else {
    # inverse sminus in b
    blc <- lapply(blc, function(it) {it$sminus <- !it$sminus; it})
    lc <- c(alc, blc)
  }
  #browser()
  # sum purely numeric terms
  inum <- which(sapply(lc, function(it) length(it$num)==0 && length(it$den)==0))
  if (length(inum) > 1) {
    term <- sum(sapply(lc[inum], function(it) (if (it$sminus) -1 else 1)*it$fa$num/it$fa$den))
    lc[[inum[1]]] <- list(fa=list(num=abs(term), den=1), sminus=term<0)
    lc <- lc[-inum[-1]]
  }
  bch <- ta <- tsim <- po <- ilc <- ind <- list()
  for (cnd in c("num", "den")) {
    # character bases in num/den
    bch[[cnd]] <- unlist(lapply(lc, function(it) {lapply(it[[cnd]]$b, format1)}))
    # powers
    po[[cnd]] <- do.call(c, lapply(lc, function(it) it[[cnd]]$p), quote=TRUE)
    # index of the lc term for each bnch
    ta[[cnd]] <- table(bch[[cnd]])
    ta[[cnd]] <- ta[[cnd]][ta[[cnd]] > 1] # keep only repeated bases
    tsim[[cnd]] <- outer(bch[[cnd]], names(ta[[cnd]]), `==`)
    ilc[[cnd]] <- unlist(lapply(seq_along(lc), function(i) {rep(i, length(lc[[i]][[cnd]]$b))}))
    # index of the base in a given term (nd) for each bnch
    ind[[cnd]] <- unlist(lapply(seq_along(lc), function(i) {seq_along(lc[[i]][[cnd]]$b)}))
  }
  #browser()
  # fnd will be the name "num" or "den" where the first factor
  # will be taken. ond is the "other" name (if fnd=="num", then ond == "den")
  # we select the candidate which is most repeated provided that it
  # has at least one numeric power occurance.
  taa <- unlist(ta)
  ota <- order(taa, decreasing=TRUE)
  ntan <- length(ta$num)
  fnd <- NA
  for (i in ota) {
    cnd <- if (i > ntan) "den" else "num"
    ita <- i - if (i > ntan) ntan else 0
    ib <- bch[[cnd]] == names(ta[[cnd]])[ita]
    if (any(sapply(po[[cnd]][ib], is.numeric))) {
      fnd <- cnd
      iit <- which(ib) # the bases equal to factor
      p_fa <- min(sapply(po[[cnd]][ib], function(p) if (is.numeric(p)) p else NA), na.rm=TRUE)
      i_lc <- ilc[[cnd]][iit]
      i_nd <- ind[[cnd]][iit]
      break
    }
  }
  #browser()
  if (is.na(fnd))
    return(lc2expr(lc, scache)) # nothing to factorize, just order terms
  ond <- if (fnd == "num") "den" else "num"
  # create nd with the first factor
  fa_nd <- list(num=list(b=list(), p=list()),
                den=list(b=list(), p=list()),
                sminus=FALSE, fa=list(num=1, den=1))
  fa_nd[[fnd]]$b <- lc[[i_lc[1]]][[fnd]]$b[i_nd[1]]
  fa_nd[[fnd]]$p <- list(p_fa)
  # decrease p in the lc terms
  for (i in seq_along(i_lc)) {
    lc[[i_lc[i]]][[fnd]]$p[[i_nd[i]]] <- Simplify_(call("-", lc[[i_lc[i]]][[fnd]]$p[[i_nd[i]]], p_fa), scache)
  }
  
  for (cnd in c(fnd, ond)) {
    # see if other side can provide factors
    for (i in seq_along(ta[[cnd]])) {
      if ((cnd == fnd && i == ita) || ta[[fnd]][ita] != ta[[cnd]][i] || any(ilc[[cnd]][tsim[[cnd]][,i]] != i_lc)) {
        next # no common layout with the factor
      }
      ib <- bch[[cnd]] == names(ta[[cnd]])[i]
      # see if it has numeric power
      if (!any(sapply(po[[cnd]][ib], is.numeric))) {
        next
      }
      iit <- which(ib) # the bases equal to factor
      p_fa <- min(sapply(po[[cnd]][ib], function(p) if (is.numeric(p)) p else NA), na.rm=TRUE)
      i_lc <- ilc[[cnd]][iit]
      i_nd <- ind[[cnd]][iit]
      fa_nd[[cnd]]$b <- append(fa_nd[[cnd]]$b, lc[[i_lc[1]]][[cnd]]$b[i_nd[1]])
      fa_nd[[cnd]]$p <- append(fa_nd[[cnd]]$p, p_fa)
      # decrease p in the lc terms
      for (i in seq_along(i_lc)) {
        lc[[i_lc[i]]][[cnd]]$p[[i_nd[i]]] <- Simplify_(call("-", lc[[i_lc[i]]][[cnd]]$p[[i_nd[i]]], p_fa), scache)
      }
    }
  }
  #browser()
  # form final symbolic expression
  # replace all i_lc by one product of fa_nd and lincomb of the reduced nds
  rest <- Simplify_(lc2expr(lc[i_lc], scache), scache)
  if (is.neg.expr(rest)) {
    rest <- negate.expr(rest)
    fa_nd$sminus <- !fa_nd$sminus
  }
  fa_nd$num$b <- append(fa_nd$num$b, rest)
  fa_nd$num$p <- append(fa_nd$num$p, 1)
  lc <- c(list(fa_nd), lc[-i_lc])
  return(lc2expr(lc, scache))
}

`Simplify.-` <- function(expr, scache=NULL)
{
  `Simplify.+`(expr, add=FALSE, scache=scache)
}

`Simplify.*` <- function(expr, div=FALSE, scache=NULL)
{
  #print(expr)
  #browser()
  a <- expr[[2]]
  b <- expr[[3]]
  if (is.uminus(a)) {
    sminus <- TRUE
    a <- a[[2]]
  } else {
    sminus <- FALSE
  }
  if (is.uminus(b)) {
    sminus <- !sminus
    b <- b[[2]]
  }
  #browser()
  if (is.numconst(a, 0) || (is.call(a) && format1(a[[1]]) %in% c("rep", "rep.int", "rep_len") && is.numconst(a[[2]], 0)) || (is.numconst(b, 0) || (is.call(b) && format1(b[[1]]) %in% c("rep", "rep.int", "rep_len") && is.numconst(b[[2]], 0)) && !div)) {
    #	if (a == 0 || (b == 0 && !div)) {
    #browser()
    0
  } else if (is.numconst(a, 1) && !div) {
    if (sminus) call("-", b) else b
  } else if (is.numconst(b, 1)) {
    if (sminus) call("-", a) else a
  } else if (div && identical(a, b)) {
    if (sminus) -1 else 1
  } else {
    #browser()
    # get numerator and denominator for a and b then combine them
    nd_a <- Numden(a)
    nd_b <- Numden(b)
    if (div) {
      nd <- list(
        num=list(b=c(nd_a$num$b, nd_b$den$b),
                 p=c(nd_a$num$p, nd_b$den$p)),
        den=list(b=c(nd_a$den$b, nd_b$num$b),
                 p=c(nd_a$den$p, nd_b$num$p))
      )
      sminus=xor(sminus, xor(nd_a$sminus, nd_b$sminus))
    } else {
      nd <- list(
        num=list(b=c(nd_a$num$b, nd_b$num$b),
                 p=c(nd_a$num$p, nd_b$num$p)),
        den=list(b=c(nd_a$den$b, nd_b$den$b),
                 p=c(nd_a$den$p, nd_b$den$p))
      )
      sminus=xor(sminus, xor(nd_a$sminus, nd_b$sminus))
    }
    # reduce numerics to only one factor
    fa=list()
    if (div) {
      fa$num <- nd_a$fa$num*nd_b$fa$den
      fa$den <- nd_a$fa$den*nd_b$fa$num
    } else {
      fa$num <- nd_a$fa$num*nd_b$fa$num
      fa$den <- nd_a$fa$den*nd_b$fa$den
    }
    res <- fa$num/fa$den
    if (all(as.integer(res) == res)) {
      fa$num <- res
      fa$den <- 1
    } else if (fa$den != 1) {
      res <- fa$den/fa$num
      if (all(as.integer(res) == res)) {
        fa$num <- 1
        fa$den <- res
      }
    }
    # group identical bases by adding their powers
    #browser()
    bch=list()
    for (na in c("num", "den")) {
      bch[[na]] <- sapply(nd[[na]]$b, format1)
      if (length(nd[[na]]$b) <= 1)
        next
      ta <- table(bch[[na]])
      ta <- ta[ta > 1]
      if (length(ta) == 0)
        next
      nd_eq <- outer(bch[[na]], names(ta), `==`)
      for (inum in seq(len=ncol(nd_eq))) {
        isim <- which(nd_eq[,inum])
        if (length(isim)) {
          # add powers for this base
          nd[[na]]$p[[isim[1]]] <- Simplify_(li2sum(nd[[na]]$p[isim]), scache)
          # set grouped powers to 0
          nd[[na]]$p[isim[-1]] <- 0
        }
      }
      # remove power==0 terms
      ize <- isim[-1]
      if (length(ize)) {
        nd[[na]]$b <- nd[[na]]$b[-ize]
        nd[[na]]$p <- nd[[na]]$p[-ize]
        bch[[na]] <- bch[[na]][-ize]
      }
    }
    # simplify identical terms in num and denum by subtracting powers
    nd_eq <- outer(bch$den, bch$num, `==`)
    ipair <- matrix(0, nrow=2, ncol=0)
    for (inum in seq(len=ncol(nd_eq))) {
      iden <- which(nd_eq[,inum]) # of length at most 1 as terms are already grouped
      if (length(iden)) {
        # simplify power for this pair
        ipair <- cbind(ipair, c(inum, iden))
        res <- Simplify_(call("-", nd$num$p[[inum]], nd$den$p[[iden]]), scache)
        if (is.neg.expr(res)) {
          nd$num$p[[inum]] <- 0
          nd$den$p[[iden]] <- negate.expr(res)
        } else {
          nd$num$p[[inum]] <- res
          nd$den$p[[iden]] <- 0
        }
      }
    }
    #browser()
    # remove power==0 terms
    for (na in c("num", "den")) {
      if (length(nd[[na]]$b) == 0)
        next
      ize=sapply(nd[[na]]$p, `==`, 0)
      nd[[na]]$b <- nd[[na]]$b[!ize]
      nd[[na]]$p <- nd[[na]]$p[!ize]
    }
    nd[["fa"]] <- fa
    nd[["sminus"]] <- sminus
    expr <- nd2expr(nd, scache)
    expr
  }
}
`Simplify./` <- function(expr, scache=NULL)
{
  `Simplify.*`(expr, div=TRUE, scache=scache)
}
`Simplify.^` <- function(expr, scache=NULL)
{
  a <- expr[[2]]
  b <- expr[[3]]
  
  if (is.numconst(a, 0)) {
    0
  } else if (is.numconst(b, 0) || is.numconst(a, 1)) {
    1
  } else if (is.numconst(b, 1)) {
    a
  } else if (identical(b, 0.5)) {
    substitute(sqrt(a))
  } else if (b == -0.5) {
    substitute(1/sqrt(a))
  } else if (is.call(a)) {
    if (a[[1]] == as.symbol("^")) {
      # product of exponents
      b <- Simplify_(call("*", a[[3]], b), scache)
      a <- a[[2]]
    } else if (a[[1]] == as.symbol("sqrt")) {
      # divide by 2
      b <- Simplify_(call("/", b, 2), scache)
      a <- a[[2]]
    } else if (a[[1]] == as.symbol("abs") && is.numeric(b) && b%%2 == 0) {
      # remove abs() for even power
      a <- a[[2]]
    }
    expr[[2]] <- a
    expr[[3]] <- b
    expr
  } else {
    expr
  }
}
Simplify.log <- function(expr, scache=NULL) {
  if (is.call(expr[[2]])) {
    # the argument of log is a function
    subf <- format1(expr[[2]][[1]])
    if (subf == "^") {
      p <- expr[[2]][[3]]
      expr[[2]] <- expr[[2]][[2]]
      expr <- Simplify_(call("*", p, expr), scache)
    } else if (subf == "exp") {
      if (length(expr) == 2)
        expr <- expr[[2]][[2]]
      else
        expr <- Simplify_(call("/", expr[[2]][[2]], call("log", expr[[3]])), scache)
    } else if (subf == "sqrt") {
      expr[[2]] <- expr[[2]][[2]]
      expr <- Simplify_(call("*", 0.5, expr), scache)
    } else if (subf == "*") {
      a <- expr
      a[[2]] <- expr[[2]][[2]]
      expr[[2]] <- expr[[2]][[3]] # unitary "+" cannot appear here
      expr <- Simplify_(call("+", a, expr), scache)
    } else if (subf == "/") {
      a <- expr
      a[[2]] <- expr[[2]][[2]]
      expr[[2]] <- expr[[2]][[3]] # unitary "+" cannot appear here
      expr <- Simplify_(call("-", a, expr), scache)
    } else if (subf == "+") {
      # replace log(1+x) by log1p(x)
      if (expr[[2]][[2]] == 1) {
        expr <- call("log1p", expr[[2]][[3]])
      } else if (expr[[2]][[3]] == 1) {
        expr <- call("log1p", expr[[2]][[2]])
      }
    }
  }
  if (length(expr) == 3 && identical(expr[[2]], expr[[3]])) {
    1
  } else {
    expr
  }
}
Simplify.sqrt <- function(expr, scache=NULL) {
  if (is.call(expr[[2]])) {
    # the argument of sqrt is a function
    subf <- format1(expr[[2]][[1]])
    if (subf == "^") {
      p <- expr[[2]][[3]]
      Simplify_(call("^",  call("abs", expr[[2]][[2]]), call("/", p, 2)), scache)
    } else if (subf == "exp") {
      expr[[2]][[2]] <- Simplify_(call("/", expr[[2]][[2]], 2), scache)
      expr[[2]]
    } else if (subf == "sqrt") {
      Simplify_(call("^", expr[[2]][[2]], 0.25), scache)
    } else if (subf == "*" && identical(expr[[2]][[2]], expr[[2]][[3]])) {
      Simplify_(call("abs", expr[[2]][[2]]), scache)
    } else {
      expr
    }
  } else {
    expr
  }
}
Simplify.abs <- function(expr, scache=NULL) {
  if (is.uminus(expr[[2]])) {
    expr[[2]] <- expr[[2]][[2]]
  } else if (is.call(expr[[2]])) {
    subf <- format1(expr[[2]][[1]])
    if (subf == "^") {
      p <- expr[[2]][[3]]
      if (is.numeric(p) && p%%2 == 0)
        expr <- expr[[2]]
    } else if (subf == "exp" || subf == "sqrt") {
      expr <- expr[[2]]
    }
  }
  expr
}
Simplify.sign <- function(expr, scache=NULL) {
  if (is.uminus(expr[[2]])) {
    expr[[2]] <- expr[[2]][[2]]
    expr <- call("-", expr)
  } else if (is.call(expr[[2]])) {
    subf <- format1(expr[[2]][[1]])
    if (subf == "^") {
      p <- expr[[2]][[3]]
      if (is.numeric(p) && p%%2 == 0)
        expr <- 1
    } else if (subf == "exp" || subf == "sqrt") {
      expr <- 1
    }
  }
  expr
}
Simplify.if <- function(expr, scache=NULL) {
  cond <- expr[[2]]
  if ((is.logical(cond) || is.numeric(cond)) && isTRUE(!!cond)) {
    expr <- expr[[3]]
  } else if (length(expr) == 4) {
    if ((is.logical(cond) || is.numeric(cond)) && isTRUE(!cond)) {
      expr <- expr[[4]]
    } else if (identical(expr[[3]], expr[[4]])) {
      expr <- expr[[3]]
    }
  }
  expr
}
Simplify.bessel <- function(expr, scache=NULL) {
  if (length(expr) < 4)
    return(expr)
  cond <- expr[[4]]
  if ((is.logical(cond) || is.numeric(cond)) && isTRUE(!cond)) {
    expr[[4]] <- NULL
  }
  expr
}
`Simplify.=` <- function(expr, scache=NULL) {
  # just strore the rhs in the scache
  if (is.symbol(expr[[2]]) && is.call(expr[[3]])) {
    scache$l[[format1(expr[[3]])]] <- expr[[2]]
  }
  expr
}
`Simplify.{` <- function(expr, scache=NULL) {
  # if the last expression is a constant just return it
  n <- length(expr)
  la <- expr[[n]]
  if (is.conuloch(la)) {
    expr <- la
  }
  expr
}

Numden <- function(expr) {
  # Return a list with "num" as numerator and "den" as denominator sublists.
  # "fa" field is for numeric factors in "num" and "den" subfields.
  # "sminus" is logical for applying or not "-" to the whole expression
  # Each sublist regroups the language expressions which are not products neither
  # divisions. The terms are decomposed in b^p sublists
  #print(expr)
  #browser()
  if (is.uminus(expr)) {
    a=Numden(expr[[2]])
    a$sminus <- !a$sminus
    a
  } else if (is.uplus(expr)) {
    Numden(expr[[2]])
  } else if (is.symbol(expr)) {
    list(num=list(b=list(expr), p=list(1)),
         sminus=FALSE,
         fa=list(num=1, den=1))
  } else if (is.numeric(expr)) {
    sminus <- length(expr) == 1 && expr < 0
    list(fa=list(num=if (sminus) -expr else expr, den=1),
         sminus=sminus)
  } else if (is.call(expr)) {
    if (expr[[1]] == as.symbol("*")) {
      # recursive call
      a=Numden(expr[[2]])
      b=Numden(expr[[3]])
      list(num=list(b=c(a$num$b, b$num$b), p=c(a$num$p, b$num$p)),
           den=list(b=c(a$den$b, b$den$b), p=c(a$den$p, b$den$p)),
           sminus=xor(a$sminus, b$sminus),
           fa=list(num=a$fa$num*b$fa$num, den=a$fa$den*b$fa$den))
    } else if (expr[[1]] == as.symbol("/")) {
      # recursive call
      a=Numden(expr[[2]])
      b=Numden(expr[[3]])
      list(num=list(b=c(a$num$b, b$den$b), p=c(a$num$p, b$den$p)),
           den=list(b=c(a$den$b, b$num$b), p=c(a$den$p, b$num$p)),
           sminus=xor(a$sminus, b$sminus),
           fa=list(num=a$fa$num*b$fa$den, den=a$fa$den*b$fa$num))
    } else if (expr[[1]] == as.symbol("^")) {
      if (is.neg.expr(expr[[3]])) {
        # make the power look positive
        list(den=list(b=list(expr[[2]]), p=list(negate.expr(expr[[3]]))),
             sminus=FALSE,
             fa=list(num=1, den=1)
        )
      } else {
        list(num=list(b=list(expr[[2]]), p=list(expr[[3]])),
             sminus=FALSE,
             fa=list(num=1, den=1)
        )
      }
    } else {
      list(num=list(b=list(expr), p=list(1)),
           sminus=FALSE,
           fa=list(num=1, den=1))
    }
  } else {
    list(num=list(b=list(expr), p=list(1)),
         sminus=FALSE,
         fa=list(num=1, den=1))
  }
}
is.uminus <- function(e) {
  # detect if e is unitary minus, e.g. "-a"
  return(is.call(e) && length(e) == 2 && e[[1]] == as.symbol("-"))
}
is.uplus <- function(e) {
  # detect if e is unitary plus, e.g. "+a"
  return(is.call(e) && length(e) == 2 && e[[1]] == as.symbol("+"))
}
is.unumeric <- function(e) {
  # detect if numeric with optional unitary sign(s)
  return(is.numeric(e) || ((is.uminus(e) || is.uplus(e)) && is.unumeric(e[[2]])))
}
is.conuloch <- function(e) {
  # detect if e is complex, numeric, logical or character
  return(is.numeric(e) || is.logical(e) || is.complex(e) || is.character(e))
}
is.neg.expr <- function(e) {
  # detect if e is a negative expression, i.e. is one of:
  #  - negative real number
  #  - unitary minus (-a)
  return((is.numeric(e) && e < 0) || is.uminus(e))
}
negate.expr <- function(e) {
  # make negative expression looking positive or inverse the difference
  if (is.numeric(e)) 
    -e
  else # e is supposed to be a unitary minus
    e[[2]]
}
is.assign <- function(e) {
  # detect if it is an assignment operator
  is.call(e) && (e[[1]] == as.symbol("<-") || e[[1]] == as.symbol("="))
}
is.subindex <- function(e) {
  # is e a simple subindex expression?
  is.call(e) && any(format1(e[[1]]) == c("$", "[", "[[")) && (is.symbol(e[[2]]) && (is.symbol(e[[3]]) || is.conuloch(e[[3]])))
}
is.numconst <- function(e, val=NULL) {
  res=is.numeric(e) && length(e) == 1L
  if (res && !is.null(val))
    res = res && is.numconst(val) && e == val
  return(res)
}
Lincomb <- function(expr) {
  # decompose expr in a list of product terms (cf Numden)
  # the sign of each term is determined by the nd$sminus logical item.
  if (is.call(expr) && length(expr) == 3) {
    if (expr[[1]] == as.symbol("+")) {
      # recursive call
      c(Lincomb(expr[[2]]), Lincomb(expr[[3]]))
    } else if (expr[[1]] == as.symbol("-")) {
      # recursive call
      a <- Lincomb(expr[[2]])
      b <- Lincomb(expr[[3]])
      # inverse the sign in b terms
      b <- lapply(b, function(it) {it$sminus <- !it$sminus; it})
      c(a, b)
    } else {
      list(Numden(expr))
    }
  } else {
    list(Numden(expr))
  }
}

# return an environement in wich stored subexpressions with
# an index giving the position of each subexpression in the
# whole statement st ("rhs" entry). Index is given as a string i1.i2.i3...
# where the integeres iN refer to st[[i2]][[i3]][[...]]
# "lhs" is index to char mapping (what is defined where)
# "def" is a mapping of lhs (char) to rhs (char)
# "{" where accolade operators are
Leaves <- function(st, ind="1", res=new.env()) {
  if (is.null(res$rhs)) {
    res$rhs <- list()
    res$lhs <- list()
    res$def <- list() # store definitions by asignments
    res[["{"]] <- list()
  }
  if (is.call(st)) {
    if (st[[1]] != as.symbol("<-") && st[[1]] != as.symbol("=")) {
      res$rhs[[ind]] <- format1(st)
      if (st[[1]] == as.symbol("{")) {
        res[["{"]] <- append(res[["{"]], ind)
      }
    } else {
      if (!is.null(res$lhs[[ind]]))
        stop("Re-assignment is not supported yet in caching.")
      if (is.call(st[[2]]))
        stop("Cannot handle yet indexing in left values.")
      lhs <- format1(st[[2]])
      res$lhs[[ind]] <- lhs # we cannot handle yet `[`, `$` etc.
      res$def[[lhs]] <- format1(st[[3]])
      # exclude this assignement from replacements if .eX
      #if (regexpr("\\.+e[0-9]+", lhs) > 0)
      #	return(res)
    }
    args <- as.list(st)[-1]
    l <- lapply(seq_along(args), function(i) Leaves(args[[i]], paste(ind, i+1, sep="."), res))
  }
  return(res)
}

# convert index calculated by Leaves() to a call like st[[i2]][[i3]]...
# the first two chars "1." are striped out
ind2call <- function(ind, st="st")
  if (ind == "1") as.symbol(st) else parse(text=sprintf("%s[[%s]]", st, gsub(".", "]][[", substring(ind, 3), fixed=TRUE)))[[1]]

# replace repeated subexpressions by cached values
# prefix is used to form auxiliary variable
##' @rdname Simplify
Cache <- function(st, env=Leaves(st), prefix="") {
  stch <- if (is.call(st)) format1(st[[1]]) else ""
  env$lhs <- unlist(env$lhs)
  #if (stch == "<-" || stch == "=") {
  #	return(call("<-", st[[2]], Cache(st[[3]], env=env, prefix=paste(".", st[[2]], sep=""))))
  #} else if (stch == "{" || stch == "c") {
  #	return(as.call(c(list(st[[1]]), lapply(as.list(st)[-1], Cache, env=env))))
  #}
  alva <- all.vars(st)
  p <- grep(sprintf("^%s.e[0-9]+", prefix), alva, value=T)
  if (nchar(prefix) == 0 && length(p) > 0) {
    prefix <- max(p)
  }
  ve <- unlist(env$rhs)
  defs <- unlist(env$def)
  tdef <- outer(ve, defs, "==")
  #browser()
  # if the subexpression is in defs, replace it with the symbol in the lhs
  for (ic in seq_len(ncol(tdef))) {
    v <- tdef[,ic]
    nme <- colnames(tdef)[ic]
    idef <- names(which(env$lhs==nme))
    for (i in which(v)) {
      ind <- names(v)[i] # subexpression index in st
      # skip self assignment
      ispl <- strsplit(ind, ".", fixed=TRUE)[[1]]
      indup <- paste(ispl[-length(ispl)], collapse=".")
      stup <- eval(ind2call(indup))
      if ((is.assign(stup) && (format1(stup[[2]]) == nme || natcompare(indup, idef) < 0)))
        next
      ve[i] <- NA
      do.call(`<-`, list(ind2call(ind), quote(as.symbol(nme))))
    }
  }
  
  suppressWarnings(ve <- ve[!is.na(ve)])
  # skip simple subindex
  isi <- sapply(ve, function(e) is.subindex(parse(text=e)[[1]]))
  ve <- ve[!isi]
  
  ta <- table(ve)
  ta <- ta[ta > 1]
  if (length(ta) == 0)
    return(st)
  e <- list() # will store the result code
  alva <- list()
  for (sub in names(sort(ta, decreasing=TRUE))) {
    # get st indexes for this subexpression
    isubs <- names(which(ve == sub))
    for (i in seq_along(isubs)) {
      isub <- isubs[i]
      subst <- ind2call(isub)
      if (i == 1) {
        esubst <- try(eval(subst), silent=TRUE)
        if (inherits(esubst, "try-error"))
          break # was already cached
        # add subexpression to the final code
        ie=length(e)+1
        estr <- sprintf("%s.e%d", prefix, ie)
        esub <- as.symbol(estr)
        e[[ie]] <- call("<-", esub, esubst)
        alva[[estr]] <- all.vars(esubst)
      }
      # replace subexpression in st by .eX
      do.call(`<-`, list(subst, as.symbol("esub")))
    }
  }
  alva[["end"]] <- all.vars(st)
  # where .eX are used? If only once, develop, replace and remove it
  wh <- lapply(seq_along(e), function(i) {
    it=sprintf("%s.e%d", prefix, i)
    which(sapply(alva, function(v) any(it == v)))
  })
  dere <- sapply(wh, function(it) if (length(it) == 1 && names(it) != "end") it[[1]] else 0)
  for (i in which(dere != 0)) {
    idest <- dere[i]
    li <- list()
    li[[sprintf("%s.e%d", prefix, i)]] <- e[[i]][[3]]
    e[[idest]][[3]] <- do.call("substitute", c(e[[idest]][[3]], list(li)))
  }
  e <- e[which(!dere)]
  #browser()
  # place auxiliary vars after the definition of the used vars
  if (stch != "{") {
    l <- c(as.symbol("{"), e, st)
    st <- as.call(l)
  } else {
    n <- length(st)
    res <- c(e, as.list(st)[-c(1, n)])
    alva <- lapply(res, all.vars)
    i <- toporder(alva)
    res <- c(as.symbol("{"), res[i], st[[n]])
    st <- as.call(res)
  }
  return(st)
}

##' @rdname Simplify
deCache <- function(st) {
  # do the job inverse to Cache(), i.e. substitute all auxiliary expressions
  # in the final one
  # NB side effect: all assignement not used in the last operation in {...} are
  # just lost.
  if (!is.call(st)) {
    return(st)
  }
  stch <- format1(st[[1]])
  stl <- as.list(st)
  if (stch == "{") {
    repl <- list()
    for (op in stl[-1]) {
      # gather substitutions
      if (is.assign(op)) {
        repl[[as.character(op[[2]])]] <- do.call("substitute", list(deCache(op[[3]]), repl))
      }
    }
    # the last operation subst
    la <- stl[[length(stl)]]
    if (is.assign(la)) {
      st <- repl[[length(repl)]]
    } else {
      st <- do.call("substitute", list(deCache(la), repl))
    }
  } else {
    # recurrsive call to deCache on all arguments of the statement
    stl <- lapply(stl, deCache)
    st <- as.call(stl)
  }
  return(st)
}

nd2expr <- function(nd, scache, sminus=NULL) {
  # form symbolic products
  # if sminus is not null, use it instead of the nd's one
  if (length(nd) == 0)
    return(0)
  eprod <- list()
  sminus <- (!is.null(sminus) && sminus) || (is.null(sminus) && nd$sminus)
  for (na in c("num", "den")) {
    if (length(nd[[na]]$b) == 0)
      next
    # alphabetic order for bases, symbols first, then calls
    cSeq <- 1:length(nd[[na]]$b)
    #cSeq <- cSeq[order(sapply(nd[[na]]$b,is.call))]
    #for (i in order(sapply(nd[[na]]$b, is.call), sapply(nd[[na]]$b, format1))) {
    for (i in cSeq) {
      #print(order(sapply(nd[[na]]$b, is.call), sapply(nd[[na]]$b, format1)))
      p <- nd[[na]]$p[[i]]
      if (p == 0)
        next
      term <- if (p == 1) nd[[na]]$b[[i]] else Simplify_(call("^", nd[[na]]$b[[i]], p), scache)
      if (term == 0)
        return(if (na == "num") 0 else if (sminus) -Inf else Inf)
      if (is.null(eprod[[na]]))
        eprod[[na]] <- term # start the sequence
      else
        eprod[[na]] <- call("*", eprod[[na]], term)
    }
  }
  expr <- if (is.null(eprod$num)) 1 else eprod$num
  if (!is.null(eprod$den)) {
    expr <- call("/", expr, eprod$den)
  }
  # put numeric factor at first place
  fa=nd$fa
  if (any(fa$num != 1) && any(fa$den != 1)) {
    # add to both num. and denom.
    if (!is.null(eprod$den)) {
      expr[[2]] <- call("*", fa$num, expr[[2]])
      expr[[3]] <- call("*", fa$den, expr[[3]])
    } else {
      expr <- call("/", call("*", fa$num, expr), fa$den)
    }
  } else if (any(fa$num != 1)) {
    if (is.call(expr) && expr[[1]] == as.symbol("/") && expr[[2]] == 1)
      expr[[2]] <- fa$num
    else if (expr == 1)
      expr <- fa$num
    else
      expr <- call("*", fa$num, expr)
  } else if (any(fa$den != 1)) {
    if (is.call(expr) && expr[[1]] == as.symbol("/"))
      expr[[3]] <- call("*", fa$den, expr[[3]])
    else
      expr <- call("/", expr, fa$den)
  }
  expr <- if (sminus) call("-", expr) else expr
  #print(sprintf("nd->%s", format1(expr)))
  return(expr)
}

lc2expr <- function(lc, scache) {
  # form symbolic sum and diff form a list of nds
  # separate in positive and negative
  smin <- sapply(lc, "[[", "sminus")
  epos <- lapply(lc[which(!smin)], nd2expr, scache)
  if (length(epos) > 1) {
    #cat("epos orig=", sapply(epos, format1), sep="\n")
    #cat("epos order=", order(sapply(epos, format1)), sep="\n")
    #cat("order - +=", order(c("-", "+")), sep="\n")
    epos <- epos[order(sapply(epos, format1), decreasing = FALSE)]
    #cat("epos=", sapply(epos, format1), sep="\n")
  }
  eneg <- lapply(lc[which(smin)], nd2expr, scache, sminus=FALSE)
  if (length(eneg) > 1) {
    eneg <- eneg[order(sapply(eneg, format1))]
  }
  if (length(epos) == 0)
    return(if (length(eneg) == 0) 0 else call("-", li2sum(eneg)))
  else
    return(if (length(eneg) == 0) li2sum(epos) else Simplify_(call("-", li2sum(epos), li2sum(eneg)), scache))
}

li2sum <- function(li) {
  # form a long sum of expressions from the list li
  len <- length(li)
  if (len == 0)
    0
  else if (len == 1)
    li[[1]]
  else if (len == 2)
    if (li[[1]] == 0)
      li[[2]]
  else if (li[[2]] == 0)
    li[[1]]
  else
    call("+", li[[1]], li[[2]])
  else
    call("+", li2sum(li[-len]), li[[len]])
}
toporder <- function(l, ind=seq_along(l), vars=sapply(l, `[[`, 1)) {
  # Topological ordering of assignement operators
  # l is a list whose memebers are resulted from all.vars(op)
  # ind is a subindexing vector for l (for recursive call)
  # vars is a vector of variable which are assigned in ops[ind]
  # return a vector of indexes like in order()
  
  # find independet assignements, i.e. whose rhs vars are not in vars
  #cat("ind=", ind, "\n")
  if (length(ind) <= 1) {
    return(ind)
  }
  rhsvar <- lapply(l[ind], `[`, -1)
  indep <- which(!sapply(rhsvar, function(v) any(v %in% vars)))
  #cat("indep=", ind[indep], "\n")
  return(c(ind[indep], toporder(l, ind[-indep], vars[-indep])))
}
natcompare <- function(s1, s2, sep="[^0-9]+") {
  # Compare two strings in natural ordering,
  # i.e. natlower("1.12", "1.2") returns 1 (i.e s1 is greater than s2)
  # while plain "1.12" < "1.2" returns TRUE
  # sep is separator for string splitting
  # By default any non number chain of characters
  # is used as a single separator and thus is exlculed
  # from comparison.
  # The fields after string splitting are compared as numerics
  # Empty string or NA are considered as -Inf, i.e. they are less
  # than any other finite number.
  # Return -1 if s1 is lower s2, 0 if s1 is equal to s2 and 1 otherwise
  # 
  v1 <- as.numeric(strsplit(s1, sep[1])[[1]])
  v1[is.na(v1)] <- -Inf
  v2 <- as.numeric(strsplit(s2, sep[1])[[1]])
  v2[is.na(v2)] <- -Inf
  l1 <- length(v1)
  l2 <- length(v2)
  lmin <- min(l1, l2)
  # complete the shortest vector by -Inf
  v1 <- c(v1, rep(-Inf, l2-lmin))
  v2 <- c(v2, rep(-Inf, l1-lmin))
  m1 <- v1 < v2
  eq <- v1 == v2
  p1 <- v1 > v2
  if (all(m1) || (any(m1) && all(!p1)) || any(head(which(m1), 1) < head(which(p1), 1))) {
    -1
  } else if (all(eq)) {
    0
  } else {
    1
  }
}
simplifications <- new.env()

assign("+", `Simplify.+`, envir=simplifications)
assign("-", `Simplify.-`, envir=simplifications)
# Rearranges the order unfornuately 
assign("*", `Simplify.*`, envir=simplifications)
assign("/", `Simplify./`, envir=simplifications)
assign("(", `Simplify.(`, envir=simplifications)
assign("^", `Simplify.^`, envir=simplifications)
assign("log", `Simplify.log`, envir=simplifications)
assign("logb", `Simplify.log`, envir=simplifications)
assign("sqrt", `Simplify.sqrt`, envir=simplifications)
assign("abs", `Simplify.abs`, envir=simplifications)
assign("sign", `Simplify.sign`, envir=simplifications)
assign("if", `Simplify.if`, envir=simplifications)
assign("besselI", `Simplify.bessel`, envir=simplifications)
assign("besselK", `Simplify.bessel`, envir=simplifications)
#assign("<-", `Simplify.=`, envir=simplifications)
#assign("=", `Simplify.=`, envir=simplifications)
assign("{", `Simplify.{`, envir=simplifications)


# Work and notes start here
library(spida2)
library(effects)
ds <- Arrests
ds$numVar <- rnorm(dim(ds[1]))
#asc <- glm(released ~ I(2*exp(age)*log(age)*age), data=ds, family=binomial)
#asc <- glm(released ~ I(exp(age) + log(age) + age) + I(2*exp(age)*log(age)):numVar + I(2*exp(age)*log(age)*age):numVar, data=ds, family=binomial)
asc <- glm(released ~ (log(age) + sex + colour)^3 + I(age^2) + I(age^2)*sex*colour + I(1.0/age) + I(exp(age) + log(age) + age) + I(exp(age)*log(age)*age):numVar, data=ds, family=binomial)
summary(asc)

# NOTE: Modified Deriv_ starting line 520 to wrap the M and I in calls rather than do the chain rule or ignore them
# NOTE: changed nd2expr line 1569 to NOT sort alphabetically, this way only the variable that is derived with respect to is moved to the front or the back only
# NOTE: changed nd2expr line 1569 to NOT sort by function call, this way only the variable that is derived with respect to is moved to the front or the back only
# NOTE: Commented out the sub expression sub in in Deriv at line 340 ish, not sure this is useful, just looks neater to have {.e1 <- log(age); .e1*exp(age)}
# Just makes it difficult to parse

reArrange <- function(dsSp, oriDexs = oriDexs, newDexs = newDexs)
{
  # Need to shuffle
  # Output always seems first or last
  # Take the value beyond the new spot for the variable up to the end of the previous spot (-1 as the oriDexs have the 1 * as the first element of the list)
  # If the derivate doesn't involve the variable we should still keep the deriviate results
  shuf <- sapply(seq_along(1:length(dsSp)),function(i) {
    if((length(newDexs[[i]]) > 0))
    {
      if(newDexs[[i]] != (oriDexs[[i]]))
      {
        if(newDexs[[i]] > oriDexs[[i]])
        {
          c(dsSp[[i]][newDexs[[i]]],dsSp[[i]][1:max(1,(oriDexs[[i]]-1))]) #dsSp[[i]][oriDexs[[i]]:max(oriDexs[[i]],(length(dsSp[[i]])-1))])
        }
        else
        {
          c(dsSp[[i]][(newDexs[[i]]+1):max(oriDexs[[i]],newDexs[[i]]+1)],dsSp[[i]][newDexs[[i]]]) #,dsSp[[i]][oriDexs[[i]]:max(length(dsSp[[i]]),oriDexs[[i]])])
        }
      }
      else
      {
        dsSp[[i]] <- dsSp[[i]]
      }
    }
    else
    {
      dsSp[[i]] <- dsSp[[i]]
    }
  } )

  rmNA <- lapply(shuf, function(x) x[!is.na(x)])

  strV <- lapply(rmNA, paste, collapse=" * ")
  res <- lapply(seq_along(1:length(strV)), function(i) {
    if(strV[[i]] == "0")
    {
      strV[[i]] = strV[[i]]
    }
    else
    {
      paste("1 * ", strV[[i]], collapse="")
    }
  }
  )

  res
}

LfxDerivate <- function(fit, vName, yName = NULL, factors = getFactorNames(fit), wrap = FALSE, debug = debug)
{
  ts <- colnames(attr(terms(fit), "factors"))
  ts <- strsplit(ts, ":")
  factors <- getFactorNames(fit)
  facl <- lapply(ts, function(x) x %in% factors | grepl(")",
                                                        x))
  ret <- lapply(seq_along(facl), function(i) {
    lapply(seq_along(facl[[i]]), function(j) if (facl[[i]][j])
      paste("M(", ts[[i]][j], ")", sep = "")
      else ts[[i]][j])
  })

  # Implicit Differentation Required
  if(!is.null(yName))
  {
    lhs <- colnames(getData(fit))[1]
    lhsD <- Deriv(lhs,yName)
    lhsD <- paste("(",lhsD,")^(-1)",sep="")
  }
  
  varList <- ret
  oriDexs <- lapply(varList, function(x) grep(vName,unlist(x)))
  varList <- lapply(varList, paste, collapse = " * ")
  varArr <- array(varList,c(1,length(varList)))

  ds <- lapply(varArr,function(x) Deriv(x,vName))
  # https://stackoverflow.com/questions/35347537/using-strsplit-in-r-ignoring-anything-in-parentheses/35347645
  dsSp <- lapply(ds, function (x) gsub("\\) \\* M","\\)  \\*  M", x))
  dsSp <- lapply(dsSp, function (x) unlist(strsplit(x, " \\* |(?>\\(.*?\\).*?\\K(, |$))", perl = TRUE)))
  dsSp <- lapply(dsSp, function(x) unlist(strsplit(x, "  \\*  ")))
  newDexs <- lapply(dsSp, function(x) grep(vName,x))

  dsFormatted <- reArrange(dsSp,oriDexs,newDexs)
  
  # Implicit Differentiation so multiple by inverse of LHS
  if(!is.null(yName))
  {
    dsFormatted <- lapply(dsFormatted,function (x) if( x != "0") paste(x," * ", lhsD, sep="") else x)
  }
  
  cbind(c(varList),c(dsFormatted))
}

vName <- "age"
dsForm <- LfxDerivate(asc, vName)
dsForm

vName <- "sex"
dsForm <- LfxDerivate(asc, vName)
dsForm

library(car)
mod <- lm(log(income) ~  education * type * women * prestige, data = Prestige)

vName <- "education"
dsForm <- LfxDerivate(mod, vName, yName = "income")
dsForm

######################################################################
# Compare to actual usage case from notes
######################################################################

LfxDerivData <- function(fit, vName, data, yName = NULL, factors = getFactorNames(fit), wrap = FALSE, debug = debug)
{
  ts <- colnames(attr(terms(fit), "factors"))
  ts <- strsplit(ts, ":")
  factors <- getFactorNames(fit)
  facl <- lapply(ts, function(x) x %in% factors | grepl(")",
                                                        x))
  ret <- lapply(seq_along(facl), function(i) {
    lapply(seq_along(facl[[i]]), function(j) if (facl[[i]][j])
      paste("M(", ts[[i]][j], ")", sep = "")
      else ts[[i]][j])
  })
  
  # Implicit Differentation Required
  if(!is.null(yName))
  {
    lhs <- colnames(getData(fit))[1]
    lhsD <- Deriv(lhs,yName)
    lhsD <- paste("(",lhsD,")^(-1)",sep="")
  }
  
  varList <- ret
  oriDexs <- lapply(varList, function(x) grep(vName,unlist(x)))
  varList <- lapply(varList, paste, collapse = " * ")
  varArr <- array(varList,c(1,length(varList)))
  
  ds <- lapply(varArr,function(x) Deriv(x,vName))
  # https://stackoverflow.com/questions/35347537/using-strsplit-in-r-ignoring-anything-in-parentheses/35347645
  dsSp <- lapply(ds, function (x) gsub("\\) \\* M","\\)  \\*  M", x))
  dsSp <- lapply(dsSp, function (x) unlist(strsplit(x, " \\* |(?>\\(.*?\\).*?\\K(, |$))", perl = TRUE)))
  dsSp <- lapply(dsSp, function(x) unlist(strsplit(x, "  \\*  ")))
  newDexs <- lapply(dsSp, function(x) grep(vName,x))
  
  dsFormatted <- reArrange(dsSp,oriDexs,newDexs)
  
  # Implicit Differentiation so multiple by inverse of LHS
  if(!is.null(yName))
  {
    dsFormatted <- lapply(dsFormatted,function (x) if( x != "0") paste(x," * ", lhsD, sep="") else x)
  }
  
  evalExpr <- do.call(cbind,lapply(dsFormatted, function (x) eval(parse(text=x),data)))
  # Add 0 column for the intercept 
  L <- cbind(0,evalExpr)
  gg <- getFix(fit)
  rownames(L) <- rownames(data)
  colnames(L) <- names(gg$fixed)
  attr(L, "data") <- data
  L
}

server <- 'blackwell.math.yorku.ca'
server <- '3.83.113.57'

dall <- read.csv(paste0("http://",server,"/data/Smoking3.csv"))
dd <- subset( dall, sex == 'BTSX')   # subset of a data frame (combined sexex)
dd$LifeExp <- dd$lifeexp.Birth # Life expectancy at birth
dd$LE <- dd$LifeExp
dd$smoke <- dd$consumption.cigPC # cigarette consumption per adult per year  
dd$HE <- dd$HealthExpPC.Tot.ppp  # health expenditures per capita in US$ PPP
dd$hiv <- dd$hiv_prev15_49  # prevalence of HIV in population 15 to 49
dd$special <- ifelse(
  dd$country %in% c('Angola','Sierra Leone','Equatorial Guinea'),
  1,
  0)  # indicator variable for 3 outlying countries
head(dd)

fit.hiv2 <- lm( LifeExp ~ log(HE) * (smoke + I(smoke^2)) + hiv+special, dd , 
                na.action = na.exclude) 

pred <- expand.grid(
  HE = c(50,150,500, 1000, 1500, 5000), 
  smoke = seq(10,2000,20), 
  hiv = 0, 
  special = 0)

vName <- "smoke"
autoL <- LfxDerivData(fit.hiv2,vName,pred)
wwA <- wald(fit.hiv2,autoL)
wwA <- as.data.frame(wwA)
head(wwA)

L <- Lfx(fit.hiv2,
         list( 0,  
               0 * M(log(HE)),
               1 ,
               1 * M(I(2*smoke)),
               0 * hiv,
               0 * special,
               1 * M(log(HE)) * 1,
               1 * M(log(HE)) * M(I(2*smoke)) 
         ), pred)

wwO <- wald(fit.hiv2, L)
wwO <- as.data.frame(wwO)
head(wwO)

sum(wwA == wwO) == sum(dim(wwO)[1] * dim(wwO)[2])

# Rough work

#vName <- "age"
#l2 <- sub("list\\(","",l)
#l3 <- sub("\\n\\)","",l2)
#l3 <- unlist(strsplit(l3,",\\n"))
#l4 <- array(l3,c(1,length(l3)))

#l5 <- apply(l4,1,function(x) sub("^[0-9]* \\* ","",x))
#exprs <- apply(l4,2,function(x) parse(text=x))
#l4[11]
#Deriv(l4[11],vName)
#splitLists <- lapply(l3, function(x) strsplit(x, " * ",fixed=T))
#oriDexs <- lapply(splitLists, function(x) grep(vName,unlist(x)))
#ds <- lapply(l4,function(x) Deriv(x,vName))
#dsExSpace <- lapply(ds, function (x) gsub("\\) \\* M","\\)  \\*  M", x))
#dsSp <- lapply(dsExSpace, function(x) unlist(strsplit(x, "  \\*  ")))
#newDexs <- lapply(dsSp, function(x) grep(vName,x))


#cbind(c(l4),c(dsForm))

#fit <- asc
#ts <- colnames(attr(terms(fit), "factors"))
#ts <- strsplit(ts, ":")
#factors <- getFactorNames(asc)
#facl <- lapply(ts, function(x) x %in% factors | grepl(")",
#                                                      x))
#ret <- lapply(seq_along(facl), function(i) {
#  lapply(seq_along(facl[[i]]), function(j) if (facl[[i]][j])
#    paste("M(", ts[[i]][j], ")", sep = "")
#    else ts[[i]][j])
#})

#varList <- ret
#oriDexs <- lapply(varList, function(x) grep(vName,unlist(x)))
#varList <- lapply(varList, paste, collapse = " * ")
#varArr <- array(varList,c(1,length(varList)))

#ds <- lapply(varArr,function(x) Deriv(x,vName))
# https://stackoverflow.com/questions/35347537/using-strsplit-in-r-ignoring-anything-in-parentheses/35347645
#dsSp <- lapply(ds, function (x) gsub("\\) \\* M","\\)  \\*  M", x))
#dsSp <- lapply(dsSp, function (x) unlist(strsplit(x, " \\* |(?>\\(.*?\\).*?\\K(, |$))", perl = TRUE)))
#dsSp <- lapply(dsSp, function(x) unlist(strsplit(x, "  \\*  ")))
#newDexs <- lapply(dsSp, function(x) grep(vName,x))

#dsFormatted <- reArrange2(dsSp,oriDexs,newDexs)
#cbind(c(varList),c(dsFormatted))

