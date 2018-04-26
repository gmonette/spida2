#' Constuctor for Positive-Definite Matrix With Zero Covariances Between
#' Predictor Random Effects
#' 
#' This function is a constructor for the \code{pdInd} class, representing a
#' positive-definite matrix with zero covariances except possibly in the first
#' row and column. If the matrix associated with \code{object} is of dimension
#' $n$, it is represented by $n + (n-1)$ unrestricted parameters representing a
#' lower-triangular log-Cholesky decomposition. The first $n$ parameters are
#' the logs of the diagonal elements of the matrix and the last $n-1$
#' components are the $n-1$ remaining elements of the lower-triangular
#' decomposition corresponding the to the possibly non-zero covariances in the
#' first row.
#' 
#' @param value an optional initialization value, which can be any of the
#' following ...
#' @param form an optional one-sided linear formula specifying the row/column
#' names for the matrix represented by \code{object}.
#' @param nam and optional vector of character strings specifying the
#' row/column names for the matrix represented by \code{object}.
#' @param data and optional data frame i which to evaluate the variables names
#' in \code{value} and \code{form}. ...
#' @param cov optional position in lower triangle of covariances that are
#' estimated and, thus, possibly non-zero. The default is that the covariances
#' in the first column are estimated and possibly non-zero.
#' @param zero optional way of specifying covariances constrained to be equal
#' to zero. Only the lower triangular portion of the \code{zero} is used. The
#' elements that are equal to 0 corresponds to the pattern of elements that are
#' constrained to zero in the covariance matrix.
#' @section CAUTION:
#' cov and zero do not work. Until fixed, \code{pdInd} only creates the
#' default covariance pattern in which the only non-zero covariances are those with the first element.
#' @export
pdInd <-
  function (value = numeric(0), form = NULL, nam = NULL,
            data = sys.parent(), cov = NULL, zero = NULL)
{
 # unchanged
 object <- numeric(0)
 class(object) <- c("pdInd", "pdMat")
 pdConstruct(object, value, form, nam, data, cov, zero)
}
#' Construct pdInd object
#' 
#' This function is a constructor for a \link{list("pdInd")} object.
#' 
#' @param object an object inheriting from the class \code{pdInd}, representing
#' a positive definite matrix with zero covariances except in the first row and
#' column.
#' @param value an option initialization value, which can be any of the
#' following ...
#' @param form an optional one-sided linear formula specifying the row/column
#' names for the matrix represented by \code{object}.
#' @param nam and optional vector of character strings specifying the
#' row/column names for the matrix represented by \code{object}.
#' @param data and optional data frame i which to evaluate the variables names
#' in \code{value} and \code{form}. ...
#' @param cov optional position in lower triangle of covariances that are
#' estimated and, thus, possibly non-zero. The default is that the covariances
#' in the first column are estimated and possibly non-zero.
#' @param zero optional way of specifying covariances constrained to be equal
#' to zero. Only the lower triangular portion of the \code{zero} is used. The
#' elements that are equal to 0 corresponds to the pattern of elements that are
#' constrained to zero in the covariance matrix.
#' @export
pdConstruct.pdInd <-
  function (object, value = numeric(0), form = formula(object),
         nam = Names(object), data = sys.parent(),
         cov = NULL, zero = NULL,...)
{
  # note that pdConstruct.pdMat return an upper-triangular R factor,
  # i.e. chol(value)
  val <- nlme:::pdConstruct.pdMat(object, value = value,
                                  form = form, nam = nam, data = data)
  if (length(val) == 0) {
    class(val) <- c("pdInd", "pdMat")
    return(val)
  }
  # mod 2015 07 04: added arbitrary cov structure of non zero
  # covariance
  isRmat <- function(x) all( x[row(x) > col(x)] == 0)
#  disp(zero)
#  disp(cov)
  if (is.matrix(val)) {
    if(!is.null(zero)) {
      zero <- zero == 0
      zero[col(zero) >= row(zero)] <- FALSE
      cov <- t(zero)
#      disp(cov)
    }
    if(is.null(cov)) cov <- (row(val) == 1) & (col(val) > 1)
#    disp(cov)
    if(isRmat(val) ){
      value <- c(log(diag(val)), val[cov])
          # keeping only the entries that should be non-zero
    } else stop("matrix should be an upper triangular matrix")
    attributes(value) <-
      attributes(val)[names(attributes(val)) != "dim"]
    attr(value,"cov") <- cov
    class(value) <- c("pdInd", "pdMat")
    attr(value,"invert") <- FALSE
    return(value)
  }
  stop("shouldn't get here in pdConstruct.pdInd")
  Ncol <- (length(val) + 1)/2
  if (length(val) != 2*round(Ncol) - 1) {
    stop(gettextf("an object of length %d does not match a pdInd factor (diagonal + covariances with intercept",
                  length(val)), domain = NA)
  }
  class(val) <- c("pdInd", "pdMat")
  val
}

#' Factor of a pdInd object.
#'
#' Function to compute the upper triangular factor of a pdInd object representing the factorization of the inverse variance matrix.
#'
#' Returns a factor for a right log-Cholesky object for positive-definite inverse variance matrix corresponding to a variance matrix with zero covariances except in the first row and column. i.e.
#' $$
#' V^{-1} = R'R
#' $$
#' with $R$ a right-triangular matrix.
#'
#' Then if the upper-diagonal elements of $R$ below the first row are all 0 then the corresponding variance matrix with will have zero covariances except on the first row (and column).
#'
#' @param object a 'pdInd' object from which the right-triangular factor of the variance matrix it represents will be extracted
#' @return the full right-triangular factor, including zeros in the lower triangle, is returned as a vector in column order
#' @examples
#' mat <- pdInd(diag(1:4))
#' pdFactor(mat)



#' Factor of a pdInd object.
#' 
#' Function to compute the upper triangular factor of a pdInd object
#' representing the factorization of the inverse variance matrix.
#' 
#' Returns a factor for a right log-Cholesky object for positive-definite
#' inverse variance matrix corresponding to a variance matrix with zero
#' covariances except in the first row and column. i.e. $$ V^-1 = t(R)R $$ with
#' $R$ a right-triangular matrix.
#' 
#' Then if the upper-diagonal elements of $R$ below the first row are all 0
#' then the corresponding variance matrix with will have zero covariances
#' except on the first row (and column).
#' 
#' @param object a 'pdInd' object from which the right-triangular factor of the
#' variance matrix it represents will be extracted
#' @return the full right-triangular factor, including zeros in the lower
#' triangle, is returned as a vector in column order
#' @examples
#' mat <- pdInd(diag(1:4))
#' pdFactor(mat)
#' @export
pdFactor.pdInd <-
function (object)
{
  invert <- attr(object,"invert")
  cov <- attr(object,"cov")
  object <- as.vector(object)
  Ncov <- sum(cov)
  Ncol <- length(object) - Ncov
# was:
#   L <- matrix(0,Ncol,Ncol)
#   diag(L) <- exp( object[1:Ncol])
#   if ( Ncol > 1 ) L[row(L)>1 & col(L)==1] <-
#     object[(Ncol+1):length(object)]
#   if(invert) c(t(solve(L))) else c(L2R(L))
  R <- matrix(0,Ncol,Ncol)
  diag(R) <- exp( object[1:Ncol])
  if ( Ncol > 1 ) R[cov] <-
    object[(Ncol+1):length(object)]
  if(invert) c(t(solve(f2L(R)))) else c(R)
}
#' pdMatrix method for pdInd objects
#' 
#' This function is the pdMatrix method for pdInd objects
#' 
#' 
#' @param object a pdInd object
#' @param factor should an upper-triangular factor be returned of the variance
#' matrix
#' @return a variance matrix or it upper-triangular factor.
#' @export
pdMatrix.pdInd <-
function (object, factor = FALSE)
{
  if (!isInitialized(object)) {
    stop("cannot extract matrix from an uninitialized object")
  }
  cov <- attr(object,"cov")
  Ncov <- sum(cov)
  Ncol <- length(object) - Ncov
  value <- array(pdFactor(object), c(Ncol, Ncol),
                 attr(object, "Dimnames"))
  ob <- as.vector(object) # subsetting object calls pdMatrix!
  attr(value, "logDet") <- 2*sum(ob[1:Ncol])
  if (factor) value else  crossprod(value)
}



#' solve method for pdInd objects.
#' 
#' This produces a pdInd object corresponding to the inverse of its argument.
#' 
#' 
#' @param a the pdInd object to invert.
#' @param b is unused but copied from \code{solve.pdLogChol}.
#' @return a pdInd object corresponding to the matrix inverse of \code{a}.
#' @export
solve.pdInd <-
function (a, b, ...)
  {
   if (!isInitialized(a)) {
     stop("cannot get the inverse of an uninitialized object")
   }
   attr(a, 'invert') <- !attr(a, 'invert')
   a
#    Ncol <- (length(a) + 1)/2
#    ob <- as.vector(a)
#    if( Ncol == 1) ret <- -ob[1]
#    else ret <-
#      c( -ob[1:Ncol] ,
#        - exp(ob[1])*ob[(Ncol+1):length(ob)]/exp(ob[2:Ncol]))
#    attributes(ret) <- attributes(a)
#    ret
}



#' Turns a factor of a PD matrix to a right-triangular factor
#' 
#' This function takes a factor of a positive-definite matrix ($V = A'A$ ) and
#' returns the equivalent right-triangular factor ($V = R'R$).
#' 
#' 
#' @param x a factor of a positive definite matrix $V$
#' @return a right-triangular factor of the same $V$.
#' @seealso \code{f2R}
#' @export
f2R <- function(A) {
  val <- qr.R(qr(A))
  sgn <- sign(diag(val))
  sgn * val
}



#' Turns a factor of a PD matrix to a left-triangular factor
#' 
#' This function takes a factor of a positive-definite matrix ($V = A'A$ ) and
#' returns the equivalent left-triangular factor ($V = L'L$).
#' 
#' @param A a factor of a positive definite matrix $V$
#' @return a right-triangular factor of the same $V$.
#' @seealso \code{f2L}
#' @export
f2L <- function(A) {
  A[] <- rev(A)
  A[] <- rev(f2R(A))
  sign(diag(A)) * A
}

### Tests
TESTS <- FALSE
if(TESTS) {



  library(nlme)
  library(magrittr)
  V <- diag(3:6)
  V[1,2:4] <- 1:3
  V[2:4,1] <- 1:3
  V
  eigen(V)  # check that it's pd
  # test f2L and f2R

  fc <- chol(V)
  (V - crossprod(fc)) %>% abs %>% max

  faceigen <- function(V) {
      ei <- eigen(V)
      sqrt(ei$values) * t(ei$vectors)
  }

  fe <- faceigen(V)
  (V - crossprod(fe)) %>% abs %>% max

  (V - crossprod(f2L(fe))) %>% abs %>% max
  (V - crossprod(f2R(fe))) %>% abs %>% max
  (V - crossprod(f2L(fc))) %>% abs %>% max
  (V - crossprod(f2R(fc))) %>% abs %>% max

  crossprod(fe) - crossprod(fc)
  crossprod(f2L(fe)) - crossprod(f2L(fc))
  f2L(fe) %>% diag %>% sign
  f2L(fe) - f2L(fc)
  f2R(fe) - f2R(fc)


  f2R(fe)
  f2R(fc)
  f2L(fe)
  # This is the real test because pdInd is initialized as the
  # inverse of a G matrix so we want the result to ...
  #

  f <- pdInd(solve(V))
  pdMatrix(solve(f))-V

  pdMatrix(solve(f),factor=T) %>% crossprod -V

  attributes(f)
  unclass(f)
  V - pdMatrix(f)
  solve(f)
  pdMatrix(f)
  solve(V)
  V-round(pdMatrix(f),6)
  pdMatrix(solve(f))

  solve(V) - pdMatrix(solve(f))
  V - pdMatrix( solve(solve(f)))
  pf <- pdFactor(f)
  pf <- matrix(pf, 4)
  crossprod(pf) - pdMatrix(f)

  # debug(pdMatrix.pdInd)
  pdMatrix(f)
  max( abs(pdMatrix(f)-V)) < 10*.Machine$double.eps
  pdMatrix(f)-V

  summary(f)


  fac <- pdMatrix(f, factor = TRUE)
  fac
  fac2 <- pdFactor(f)
  fac2
  all(c(fac) == fac2)

  max(abs(crossprod( pdMatrix(f, factor = TRUE)) - V)) < 10*.Machine$double.eps

  finv <- solve(f)
  unclass(f)
  unclass(finv)
  f.mat <- pdMatrix(f)
  finv.mat <- pdMatrix(finv)
  zapsmall(f.mat)
  zapsmall(finv.mat)
  zapsmall(f.mat %*% finv.mat)
  diag(f.mat)
  diag(finv.mat)

# test with lme
  # library(spida)
  head(hs)
  hs$sex <- 1*(hs$Sex == 'Female')
fit1 <- lme(mathach ~ ses + sex, hs,
            random = list( school = pdDiag( ~ 1 + ses + sex)))
fit2 <- lme(mathach ~ ses + sex, hs,
            random = list( school = pdSymm( ~ 1 + ses + sex)),
            control = list(msVerbose=T,returnObject=T))
fit3 <- lme(mathach ~ ses + sex, hs,
            random = list( school = pdLogChol( ~ 1 + ses + sex)),
            control = list(msVerbose=T,returnObject=T))
hs$sex <- with(hs, 1*(Sex == "Male"))
fit4 <- lme(mathach ~ ses + sex, hs,
            random = list( school = pdInd( ~ 1 + ses + sex)),
            control = list(msVerbose=T,returnObject=T,msMaxIter=1000))
fit5 <- lme(mathach ~ ses + Sex, hs,
            random = list( school = pdInd( ~ 1 + ses + Sex)),
            control = list(msVerbose=T,returnObject=T,msMaxIter=1000))
summary(fit4)
summary(fit5)
summary(fit2)
ff4 <- (fit4$modelStruct$reStruct$school)
solve(ff4)
summary(ff4)
unclass(ff4)
summary(solve(ff4))

ff3 <- fit3$modelStruct$reStruct$school
ff2 <- fit2$modelStruct$reStruct$school
summary(ff3)
summary(solve(ff3))
summary(ff3)
summary(ff4)
AIC(fit3)
AIC(fit4)
AIC(fit2)
AIC(fit1)

unclass(ff2)
unclass(ff3)

}
