
ker <- function(L, tol = 1e-14) {
  sv <- svd(t(L), 
             nu = NCOL(L))
  d <- sum(sv$d > tol)
  ret <- sv$u[,-seq_len(d)]
  attr(ret, 'd') <- sv$d
  ret
}



ccor <- function(x,y) {
  # canonical correlations to test whether two
  # matrices span the same or orthogonal space,
  # also the dimension of the intersection
  cancor(x,y,F,F)$cor
}

ob <- function(x,  tol = 1e-14) {
  # orthogonal basis
  sv <- svd(x)
  disp(sum(sv$d > tol))
  ret <- sv$u[,seq_len(sum(sv$d > tol))]
  attr(ret, 'd') <- sv$d
  ret
} 

cc_svd <- function( X , Z = diag( NROW(X)) , ip = diag( NROW(X)), tol = 1e-14 ) {
  #
  # conjugate complement using the SVD
  ret <- Z %*% ker(t(X) %*% ip %*% Z)
  ret
}

cc_qr <- function( X , Z = diag( NROW(X)) , ip = diag( NROW(X)), tol = 1e-14 ) {
  # conjugate complement using QR
  
    # also in package spida:
  # keep versions consistent
  # ConjComp returns a basis for the conjugate complement of the
  # conjugate projection of X into span(Z) with respect to inner product with
  # matrix ip.
  # Note: Z is assumed to be of full column rank but not necessarily X.
  xq <- qr(t(Z) %*% ip %*% X, tol = tol)
  if ( xq$rank == 0 ) return( Z )
  a <- qr.Q( xq, complete = TRUE ) [ ,-(1:xq$rank)]
  Z %*% a
}

cbind(1, 1:5, (1:5)^2, (1:5-10000000000)^2) %>% 
  { ccor(cc_svd(.), cc_qr(.)) -1
  }
cc_svd(Xc)

cc_qr(Xc)

cancor(cbind(1,1:5,(1:5)^2), cbind(11:15,5:1,(1:5-3)^2),F,F)$cor -1
library(car)
library(spida2)
dd <- Prestige 



ob(Xc)




L <- rbind(1, 1:5, 11:15)
ker(L)

ccor(ker(L), t(L))
ConjComp
ccor(ConjComp(t(L)), ker(L))



ConjComp
