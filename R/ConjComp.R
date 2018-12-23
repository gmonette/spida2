###
### modified 12/23/2018 2:51 PM 
### - uses svd instead of qr
###

#' Conjugate complement of span(X) in span(Z) with respect to inner product ip
#' 
#' Generates a matrix that is a basis for the conjugate complement, i.e. the orthogonal complement
#' with respect to a norm given by an inner-product matrix so that the norm may be more general
#' than the default euclidean norm.
#' 
#' @param X matrix generating space whose complement is needed
#' @param Z matrix generating space within which complement is neeeded
#' @param ip inner-product matrix
#' @param tol tolerance, default: 1e-07
#' @example {
#' X <- cbind(c(1,1,1,1,1), c(2,-1,-1,0,0), c(3,0,0,1,1))
#' X
#' Z <- cbind(diag(5)[,1:3], c(1,1,1,0,0))
#' Z 
#' # These two matrices illustrate some issues 
#' # related to finding an orthogonal complement:
#' # 1) both matrices are column rank deficient: rank(X) = 2 and rank(Z) = 3 
#' # 2) a subspace of span(X), spanned by c(2,-1,-1,0,0), is orthogonal to span(Z)
#' # Thus the projection of span(X) into span(Z) has rank 1 and its orthogonal
#' # complement in span(Z) has rank 2. 
#' rk <- function(m, tol = 10*.Machine$double.eps) sum(svd(m,0,0)$d > tol)
#' rk(X)
#' rk(Z)
#' CC <- ConjComp2(X, Z)
#' CC
#' rk(CC)
#' CC <- ConjComp2(X, Z, tol = 1e-14)
#' CC
#' }
#' @export
ConjComp <- function( X , Z = diag( nrow(X)) , ip = diag( nrow(X)), tol = 1e-07 ) {
  xq <- qr(t(Z) %*% ip %*% X, tol = tol)
  if ( xq$rank == 0 ) return( Z )
  a <- qr.Q( xq, complete = TRUE )[ , -(1:xq$rank)]
  Z %*% a
}
#' @describeIn ConjComp based on svd
#' @export
ConjComp2 <- function(X, Z = NULL, ip = NULL, tol = .Machine$double.eps, full.rank = FALSE) {
  D <-  as.matrix(X)
  if(!is.null(ip)) D <- ip %*% D
  if(!is.null(Z)) D <- t(Z) %*% D
  sv <- svd(D, nu = nrow(D), nv = 0)
  d <- c(sv$d, rep(0, nrow(D) - length(sv$d)))
  d <- d/d[1]
  disp(d)
  a <- sv$u[, d < tol, drop = FALSE]
  disp(a)
  if(ncol(a) == 0) return(matrix(0, nrow(a), 1))
  if(full.rank) {
    sv <- svd(Z %*% a, nv = 0)
    d <- sv$d
    d <- d/d[1]
    sv$u[, d > tol, drop = FALSE]
  }
  if(is.null(Z)) a
  else Z %*% a
}
ConjComp2(diag(3))

ConjComp2(X,Z, tol = 1e-17)
X
X[,-2]
Z
ConjComp2(X[,-2],Z, tol = 1e-17)
ConjComp2(X[,-c(2,3)],Z, tol = 1e-17)
ConjComp2(X[,-c(2,3)],Z, tol = 1e-14, full = T) %>% svd %>% {.$d}
ConjComp2(X[,-c(2,3)],Z, tol = 1e-17, full = T) 
ConjComp2(X,Z, full.rank = T)
ConjComp2(diag(3)[,-3])
ConjComp2(diag(3)[,-3], cbind(c(1,1,1)))
ConjComp2(X,Z[,-4])
Proj(ConjComp(X,Z)) - Proj(ConjComp2(X,Z))  %>% rk
rk(ConjComp2(X,Z))
rk(ConjComp2(Z,X))
X
sum(diag(Proj(ConjComp(X,Z))))
- Proj(ConjComp2(X,Z))
svd(ConjComp(X,Z)) 
svd(ConjComp2(X,Z))
rk(ConjComp(X,Z)) 
rk(Proj(ConjComp(X,Z))) 
rk(ConjComp2(X,Z))

library(rbenchmark)
install.packages('rbenchmark')

system.time(replicate(10000, ConjComp(X,Z)))
system.time(replicate(10000, ConjComp2(X,Z)))
system.time(replicate(10000, ConjComp2(X,Z, full.rank = T)))


#' @describeIn ConjComp Orthogonal complement of space(X) in span(Z)
#' @export
OrthoComp <-
  function (X, Z , tol = 1e-07) {
    ConjComp(X, Z, tol = tol)
  }
#' @describeIn ConjComp Orthogonal complement of space(X) in span(Z)
#' @export
Proj <-
function(X, tol = 10*.Machine$double.eps) {
  # projection matrix onto span(x)
  sv <- svd(X, nv=0)
  d <- sv$d
  if(d[1] == 0) return(matrix(0, nrow(X), nrow(X)))
  d <- d/d[1]
  u <- sv$u[,d > tol] 
  tcrossprod(u)
}
#' @describeIn ConjComp Projection into orthogonal complement of span(X) 
#' @export
oProj <-
function(X,...) {
  diag(nrow(X)) - Proj(X,...)
}

############
#############

c("tolNorm2", "qr.R", "qr",
                      "useGrad", "maybeGrad") %>% 
  lapply(function(meth) rankMatrix(Hilbert(16),meth))

rk <- function(m, tol =10*.Machine$double.eps) sum(diag(Proj(m, tol = tol)))
rk(Hilbert(18), tol = .001*.Machine$double.eps)
rk(Hilbert(11))

#' X <- cbind(c(1,1,1,1,1), c(2,-1,-1,0,0), c(3,0,0,1,1))
#' X
#' Z <- cbind(diag(5)[,1:3], c(1,1,1,0,0))
#' Z 
#' rk(X)
#' rk(Z)
#' rk.int(X,Z)
#' rk.union(X,Z)

X
Z
rk(X)
rk(Z)
cancor(X,Z,F,F) # perhaps confused by matrices of non-li columns
X %*% cancor(X,Z,F,F)$xcoef
ConjComp2(X,Z)

svd(t(X) %*% Z)
