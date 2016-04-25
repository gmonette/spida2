###
### ConjComp
###

#' Conjugate complement of span(X) in span(Z) with respect to inner product ip
#'
#' @param X
#' @param Z
#' @param ip
#' @param tol
#' @export
ConjComp <- function( X , Z = diag( nrow(X)) , ip = diag( nrow(X)), tol = 1e-07 ) {
  help <- "
  ConjComp returns a basis for the conjugate complement of the
  conjugate projection of X into span(Z) with respect to inner product with
  matrix ip.
  Note: Z is assumed to be of full column rank but not necessarily X.
  "

  xq <- qr(t(Z) %*% ip %*% X, tol = tol)
  if ( xq$rank == 0 ) return( Z )
  a <- qr.Q( xq, complete = TRUE ) [ ,-(1:xq$rank)]
  Z %*% a
}

#' Orthogonal basis for the orthogonal complement of the column space of a matrix.
#'
#' @param X
#' @param Z
#' @param tol
#' @export
OrthoComp <-
  function (X, Z , tol = 1e-07) {
    ConjComp(X, Z, tol = tol)
  }

