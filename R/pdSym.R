# pdSymm -> pdSym
# to test definition of new pdClass
# This is a copy of pdSymm methods
# MAKE A COPY BEFORE MODIFYING
#
#' Copy of pdSymm to debug pdInd
#' 
#' @export  
pdSym <- function (value = numeric(0), form = NULL, nam = NULL, data = parent.frame()) 
{
  object <- numeric(0)
  class(object) <- c("pdSymm", "pdMat")
  pdConstruct(object, value, form, nam, data)
}
# <bytecode: 0x58284f6ac9f8>
#   <environment: namespace:nlme>
#' @export  
coef.pdSym <- function (object, unconstrained = TRUE, ...) 
{
  if (unconstrained || !isInitialized(object)) 
    NextMethod()
  else {
    val <- as.matrix(object)
    aN <- Names(object)
    aN1 <- paste("cov(", aN, sep = "")
    aN2 <- paste(aN, ")", sep = "")
    aNmat <- t(outer(aN1, aN2, paste, sep = ","))
    aNmat[row(aNmat) == col(aNmat)] <- paste("var(", aN, 
                                             ")", sep = "")
    val <- val[row(val) <= col(val)]
    names(val) <- aNmat[row(aNmat) <= col(aNmat)]
    val
  }
}
#' @export  
Dim.pdSym <- function (object, ...) 
{
  if (isInitialized(object)) {
    val <- round((sqrt(8 * length(object) + 1) - 1)/2)
    c(val, val)
  }
  else {
    NextMethod()
  }
}

#' @export  
logDet.pdSym <- function (object, ...) 
{
  if (!isInitialized(object)) {
    stop("cannot extract the log of the determinant from an uninitialized object")
  }
  attr(pdMatrix(object, factor = TRUE), "logDet")
}

#' @export
pdConstruct.pdSym <- function (object, value = numeric(0), form = formula(object), 
          nam = Names(object), data = parent.frame(), ...) 
{
  val <- NextMethod()
  if (length(val) == 0) {
    class(val) <- c("pdSymm", "pdMat")
    return(val)
  }
  if (is.matrix(val)) {
    vald <- svd(val, nu = 0)
    object <- vald$v %*% (log(vald$d) * t(vald$v))
    value <- object[row(object) <= col(object)]
    attributes(value) <- attributes(val)[names(attributes(val)) != 
                                           "dim"]
    class(value) <- c("pdSymm", "pdMat")
    return(value)
  }
  Ncol <- round((sqrt(8 * length(val) + 1) - 1)/2)
  if (length(val) != round((Ncol * (Ncol + 1))/2)) {
    stop(gettextf("an object of length %d does not match the required parameter size", 
                  length(val)), domain = NA)
  }
  class(val) <- c("pdSymm", "pdMat")
  val
}
#' @export
pdFactor.pdSym <- function (object) 
{
  Ncol <- round((-1 + sqrt(1 + 8 * length(object)))/2)
  .C(matrixLog_pd, Factor = double(Ncol * Ncol), as.integer(Ncol), 
     as.double(object))$Factor
}

#' @export
pdMatrix.pdSym <- function (object, factor = FALSE) 
{
  if (!isInitialized(object)) 
    stop("cannot extract matrix from an uninitialized object")
  if (factor) {
    Ncol <- Dim(object)[2]
    value <- array(pdFactor(object), c(Ncol, Ncol), attr(object, 
                                                         "Dimnames"))
    attr(value, "logDet") <- sum(log(abs(svd.d(value))))
    value
  }
  else {
    NextMethod()
  }
}

#' @export
solve.pdSym <- function (a, b, ...) 
{
  if (!isInitialized(a)) {
    stop("cannot extract the inverse from an uninitialized object")
  }
  coef(a) <- -coef(a, TRUE)
  a
}

#' @export
summary.pdSym <- function (object, structName = "General positive-definite", ...) 
{
  summary.pdMat(object, structName)
}

if(FALSE) {
  # tests
  hs$sex <- 1*(hs$Sex == 'Female')
  fit <- lme(mathach ~ ses + sex, hs, random = list(school = pdSym(~ses + sex)))
}