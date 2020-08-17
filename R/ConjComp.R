###
### ConjComp
###

#' Conjugate complement of span(X) in span(Z) with respect to inner product ip
#'
#' @param X matrix defining space whose complement is computed.
#'          Not necessarily of full column rank
#' @param Z matrix defining space within which complement is computed.
#'          Should be of full column rank. Default: \code{diag( nrow(X))}
#' @param ip positive definite matrix defining inner product with respect to which
#'           complement is computed. Default: \code{diag( nrow(X))}
#' @param tol tolerance (default 1e-07)
#' @export
ConjComp <- function( X , Z = diag( nrow(X)) , ip = diag( nrow(X)), tol = 1e-07 ) {
  xq <- qr(t(Z) %*% ip %*% X, tol = tol)
  if ( xq$rank == 0 ) return( Z )
  a <- qr.Q( xq, complete = TRUE ) [ ,-(1:xq$rank)]
  Z %*% a
}

#' Orthogonal basis for the orthogonal complement of the column space of a matrix.
#'
#' @param X matrix defining space whose orthogonal complement is computed.
#'          Not necessarily of full column rank
#' @param Z matrix defining space within which orthogonal complement is computed.
#'          Should be of full column rank. Default: \code{diag( nrow(X))}
#' @param tol tolerance (default 1e-07)
#' @export
OrthoComp <-
  function (X, Z = diag(nrow(X)), tol = 1e-07) {
    ConjComp(X, Z, tol = tol)
  }