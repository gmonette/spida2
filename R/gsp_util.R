##  gsp-util.R
##
##  General polynomial splines: August 8, 2008
##  Revised June 16, 2010 to incorporate degree 0 polynomials
##    and constraining evaluation of spline to 0 with 'intercept'
##   
##  Modified: June 15, 2013 - GM
##  Added specified linear contraints to gsp and 
##  to sc expressed
##  in terms of full polynomial parametrization
##  
##  Added 'PolyShift' function to specify constraints for
##  a periodic spline. See example in the manual page.
##
##  Added 'periodic' argument to gsp
##
##
#' @export
Xmat <- function( x, degree, D = 0, signif = 3) {
  
  # Return design matrix for d-th derivative
  
  if ( length(D) < length( x ) ) D = rep(D, length.out = length(x))
  if ( length(x) < length( D ) ) x = rep(x, length.out = length(D))
  xmat = matrix( x, nrow = length(x), ncol = degree + 1)
  # disp( d)
  # disp( x)
  expvec <- 0:degree
  coeffvec <- rep(1, degree+1)
  expmat <- NULL
  coeffmat <- NULL
  
  for (i in 0:max(D) ) {
    expmat <- rbind(expmat, expvec)
    coeffmat <- rbind(coeffmat, coeffvec)
    coeffvec <- coeffvec * expvec
    expvec <- ifelse(expvec > 0, expvec - 1, 0)
  }
  G = coeffmat[ D + 1, ] * xmat ^ expmat[ D + 1, ]
  
  xlab = signif( x, signif)
  rownames(G) = ifelse( D == 0, paste('f(', xlab ,')', sep = ''),
                        paste( "D",D,"(",xlab,")", sep = ""))
  colnames(G) = paste( "X", 0:(ncol(G)-1), sep = "")
  G
}

#  Xmat( c(1:10,pi), 5,0:4,3)

#' @export
Xf <-  function(   x, knots, degree = 3, D = 0, right = TRUE , signif = 3) {
  
  # With the default, right ,== TRUE, if x is at a knot then it is included in
  # in the lower interval. With right = FALSE, it is included in the higher
  # interval. This is needed when building derivative constraints at the
  # knot
  
  xmat = Xmat ( x, degree, D , signif )
  k = sort(knots)
  g = cut( x, c(-Inf, k, Inf), right = right)
  do.call( 'cbind',
           lapply( seq_along(levels(g)) , function( i ) (g ==levels(g)[i]) *  xmat))
}

# Xf( c(1:10,pi), c(3,6))

# Xf( 3, c(3,6),3, 0, )
# Xf( 3, c(3,6),right = F)

#' @export
Cmat <- function( knots, degree, smooth, lin = NULL, intercept = 0, signif = 3) {
  # add lin: contraints, GM 2013-06-13
  # generates constraint matrix
  dm = max(degree)
  
  # intercept
  
  cmat = NULL
  if( !is.null(intercept))  cmat = rbind( cmat, "Intercept" =
                                            Xf( intercept, knots, dm, D=0 ))
  
  # continuity constraints
  
  for ( i in seq_along(knots) ) {
    k = knots[i]
    sm = smooth[i]
    if ( sm > -1 ) {  # sm = -1 corresponds to discontinuity
      dmat  =Xf( k , knots, dm, D = seq(0,sm) , F ) -   Xf( k ,
                                                            knots, dm, D = seq(0,sm) , T )
      rownames( dmat ) = paste( "C(",signif(k, signif),").",
                                seq(0,sm), sep = '')
      cmat = rbind( cmat,  dmat)
    }
  }
  
  # reduced degree constraints
  
  degree = rep( degree, length.out = length(knots) + 1)
  # disp ( degree )
  for ( i in seq_along( degree)) {
    di = degree[i]
    
    if ( dm > di ) {
      dmat = diag( (length(knots) + 1) * (dm +1)) [  (i - 1)*(dm + 1) +
                                                       1 + seq( di+1,dm), , drop = F]
      rownames( dmat ) = paste( "I.", i,".",seq(di+1,dm),sep = '')
      #disp( cmat)
      #disp( dmat)
      cmat = rbind( cmat, dmat)
    }
  }
  
  # add additional linear constraints
  
  if( ! is.null(lin)) cmat <- rbind(cmat,lin) # GM:2013-06-13
  rk = qr(cmat)$rank
  spline.rank = ncol(cmat) - rk
  attr(cmat,"ranks") = c(npar.full = ncol(cmat), C.n = nrow(cmat),
                         C.rank = rk , spline.rank = spline.rank )
  attr(cmat,"d") = svd(cmat) $ d
  cmat
  
}

# Cmat( c(-.2,.4), c(2,3,2), c(2,2))


#' @export
Emat <- function(  knots, degree, smooth , intercept = FALSE, signif = 3) {
  # estimation matrix:
  # polynomial of min
  # note that my intention of using a Type II estimate is not good
  #    precisely because that means the parametrization would depend
  #    on the data which I would like to avoid
  #    BUT perhaps we can implement later
  if ( length( degree) < length(knots) + 1) degree = rep( degree,
                                                          length.out =  length(knots) + 1)
  dmin = min(degree)
  # if ( dmin < 1) stop("minimum degree must be at least 1")
  dmax = max(degree)
  smin = min(smooth)
  smax = max(smooth)
  imax = length(degree)
  
  zeroi = as.numeric(cut( 0, c(-Inf, sort( knots), Inf)))
  dzero = degree[zeroi]
  
  cmat = Xf( 0, knots, degree = dmax , D = seq( 1, dzero))
  
  if ( imax > zeroi ){
    for ( i in (zeroi+1): imax) {
      d.right = degree[ i ]
      d.left  = degree[ i-1 ]
      k = knots[ i - 1 ]
      sm = smooth[ i - 1 ]
      if (  d.right > sm ) {
        dmat  =Xf( k , knots, dmax, D = seq(sm+1,d.right) , F ) -
          Xf( k , knots, dmax, D = seq(sm+1,d.right) , T )
        rownames( dmat ) = paste( "C(",signif(k, signif),").",
                                  seq(sm+1,d.right), sep = '')
        cmat = rbind( cmat,  dmat)
        
      }
    }
  }
  
  if ( zeroi > 1 ){
    for ( i in zeroi: 2) {
      d.right = degree[ i ]
      d.left  = degree[ i-1 ]
      k = knots[ i - 1 ]
      sm = smooth[ i - 1 ]
      if (  d.left > sm ) {
        dmat = Xf( k , knots, dmax, D = seq(sm+1,d.left) , F ) -
          Xf( k , knots, dmax, D = seq(sm+1,d.left) , T )
        rownames( dmat ) = paste( "C(",signif(k, signif),").",
                                  seq(sm+1,d.left), sep = '')
        cmat = rbind( cmat,  dmat)
        
      }
    }
  }
  cmat
}

#' Orthogonal basis for column space of a matrix
#' 
#' @param X 
#' @param coef 
basis  <- function( X , coef = FALSE ) {
  # returns linear independent columns
  #
  #   X = fr(X) %*% attr(fr(X),"R")
  #
  q <- qr(X)
  sel <- q$pivot[ seq_len( q$rank)]
  ret <- X[ , sel ]
  attr(ret,"cols") <- sel
  if ( coef )attr(ret,"R") <- qr.coef( qr(ret) , mat)
  colnames(ret) <- colnames(X)[ sel ]
  ret
}

#    tt( c(-1,1,2), c(3,3,3,2), c(3,3,2))

#    tmat = rbind( c(1,1,1,1,1), c(1,0,0,0,1), c(0,1,1,1,0), c(1,2,0,0,1))
#    t(tmat)
#    basis(tmat)
#    basis( t(tmat))



#' Spline tranformation matrix
#'
#' @param knots 
#' @param degree 
#' @param smooth 
#' @param intercept 
#' @param signif 
#' 
#' @export
spline.T <-
  function( knots, degree, smooth, lin = NULL, intercept = 0, signif = 3 ) {
    cmat = Cmat( knots, degree, smooth, lin, intercept, signif  )  # constraint matrix
    emat = Emat( knots, degree, smooth, !is.null(intercept), signif  )  # estimation matrix
    #disp( list(C= cmat, E=emat ) )
    nc = nrow( cmat )
    ne = nrow( emat )
    basisT = t( basis( cbind( t(cmat), t(emat) ) ))
    cols = attr(basisT,"cols")
    ncc = sum( cols <= nc )
    Tmat = solve( basisT )
    attr( Tmat, "nC") =  ncc
    attr( Tmat, "CE") = rbind( cmat, emat)
    attr( Tmat, "CEbasis" ) = t(basisT)
    Tmat
  }





#' Spline estimation function
#' 
#' @param knots 
#' @param degree 
#' @param smooth 
#' @param intercept 
#' @param signif 
#' 
#' @export
spline.E <- function( knots, degree, smooth, lin = NULL, intercept = 0, signif = 3 ) {
  cmat = Cmat( knots, degree, smooth, lin, intercept, signif  )  # constraint matrix
  emat = Emat( knots, degree, smooth, !is.null(intercept), signif  )  # estimation matrix
  # disp( list(C= cmat, E=emat ) )
  nc = nrow( cmat )
  ne = nrow( emat )
  basisT = t( basis( cbind( t(cmat), t(emat) ) ))
  cols = attr(basisT,"cols")
  ncc = sum( cols <= nc )
  Tmat = solve( basisT )
  Tmat[, (ncc+1): ncol(Tmat)]
}





#' Projection matrix
#' 
#' Projection matrix into column space of a matrix. Works even if the matrix
#' does not have full column rank.
#' 
#' @param x 
#' @export
Proj.1 <- function( x)  x %*% ginv(x)



#' Projection matrix
#' 
#' Matrix of orthogonal projection into column space of matrix.
#' 
#' @param x 
#' @export
Proj <- function(x) {
  # projection matrix onto span(x)
  u <- svd(x,nv=0)$u
  u %*% t(u)
  
}



#' Test accuracy of projection method.
#' 
#' @param x 
#' @param fun 
#' @export
Proj.test <- function( x, fun = Proj) {
  help = "
  # sample test:
  zz <- matrix( rnorm(120), ncol = 3) %*% diag(c(10000,1,.000001))
  Proj.test( zz, fun = Proj)
  Proj.test( zz, fun = Proj.1)
  "
  # limited testing suggests Proj is generally
  # roughly as good as Proj.1  with errors 1e-17
  # but Proj.1 occasionally has errors ~ 1e-11
  pp <- fun(x)
  ret <- pp%*%pp - pp
  attr(ret,"maxabs") <- max(abs(ret))
  attr(ret,"maxabs")
}




#  round(basis( Proj(tmat)),3)
#  round(basis( Proj(t(tmat))),3)



#' Marked for deletion
#' 
#' @param x 
#' @param knots 
#' @param degree 
#' @param smooth 
#' @param intercept 
#' @param signif
#' @export
gspf.1 <- function( x, knots, degree= 3, smooth = pmax(pmin( degree[-1],
                                                             degree[ - length(degree)]) - 1,0 ), intercept = 0, signif = 3) {
  # generates a 'full' spline matrix
  tmat = spline.T(knots, degree, smooth, intercept = intercept, signif = signif)
  # put E first then C  and intercept last
  nC = attr( tmat, "nC")
  nE = ncol(tmat) - nC
  tmat = tmat[, c(seq(nC+1, ncol(tmat)), 2:nC, 1)]  # first nE columns are E matrix
  # tmat = tmat[ , - ncol(tmat)]   # drop intercept columns
  # full matrix
  X = Xf ( x, knots, max(degree), signif = signif) %*% tmat
  qqr = qr(X)
  Q = qr.Q( qqr )
  R = qr.R( qqr )
  # we need to form the linear combination of
  
  gsp( seq(-2,10), 0, c(2,2), 0)
  ret
}
#' Generate a smoothing spline
#' 
#' @param formula 
#' @param data 
#' 
#' @export
smspline <-
  function (formula, data)
  {
    if (is.vector(formula)) {
      x <- formula
    }
    else {
      if (missing(data)) {
        mf <- model.frame(formula)
      }
      else {
        mf <- model.frame(formula, data)
      }
      if (ncol(mf) != 1)
        stop("formula can have only one variable")
      x <- mf[, 1]
    }
    x.u <- sort(unique(x))
    Zx <- smspline.v(x.u)$Zs
    Zx[match(x, x.u), ]
  }
#' smspline.v from library splines with better ginverse
#' 
#' @param time 
#' 
#' @export
smspline.v <-
  function (time)
  {
    t1 <- sort(unique(time))
    p <- length(t1)
    h <- diff(t1)
    h1 <- h[1:(p - 2)]
    h2 <- h[2:(p - 1)]
    Q <- matrix(0, nr = p, nc = p - 2)
    Q[cbind(1:(p - 2), 1:(p - 2))] <- 1/h1
    Q[cbind(1 + 1:(p - 2), 1:(p - 2))] <- -1/h1 - 1/h2
    Q[cbind(2 + 1:(p - 2), 1:(p - 2))] <- 1/h2
    Gs <- matrix(0, nr = p - 2, nc = p - 2)
    Gs[cbind(1:(p - 2), 1:(p - 2))] <- 1/3 * (h1 + h2)
    Gs[cbind(1 + 1:(p - 3), 1:(p - 3))] <- 1/6 * h2[1:(p - 3)]
    Gs[cbind(1:(p - 3), 1 + 1:(p - 3))] <- 1/6 * h2[1:(p - 3)]
    Gs
    Zus <- t(ginv(Q))  # orig: t(solve(t(Q) %*% Q, t(Q)))
    R <- chol(Gs, pivot = FALSE)
    tol <- max(1e-12, 1e-08 * mean(diag(R)))
    if (sum(abs(diag(R))) < tol)
      stop("singular G matrix")
    Zvs <- Zus %*% t(R)
    list(Xs = cbind(rep(1, p), t1), Zs = Zvs, Q = Q, Gs = Gs,
         R = R)
  }
#' Approximation for semi-parametric splines
#' 
#' @param Z 
#' @param oldtimes 
#' @param newtimes 
#' 
#' @export
approx.Z <-
  function (Z, oldtimes, newtimes)
  {
    oldt.u <- sort(unique(oldtimes))
    if (length(oldt.u) != length(oldtimes) || any(oldt.u != oldtimes)) {
      Z <- Z[match(oldt.u, oldtimes), ]
      oldtimes <- oldt.u
    }
    apply(Z, 2, function(u, oldt, newt) {
      approx(oldt, u, xout = newt)$y
    }, oldt = oldtimes, newt = newtimes)
  }


#' @export
smsp <- function( x, knots ) {
  # generates matrix for smoothing spline over x with knots at knots
  #  Usage:   e.g. from example in smspline
  # lme ( ....., random=list(all= pdIdent(~smsp( x, 0:100) - 1),
  #              plot= pdBlocked(list(~ days,pdIdent(~smsp( x, 0:100) - 1))),
  #              Tree = ~1))
  
  if( max(knots) < max(x) || min(x) < min(knots)) warning( "x beyond range of knots")
  if( length( knots ) < 4 ) warning( "knots should have length at least 4")
  sp = smspline( knots)
  approx.Z( sp, knots, x)
}

#' Transformation of polynomial coefficients to change origin
#'  
#' Useful to equate polynomial coefficients for a periodic spline
#' @param x 
#' @param n
#' @examples
#' coefs <- c(3,2,4)
#' x <- 3
#' n <- length(coefs)
#' a <- 2
#' sum( coefs * x ^ (0:(n-1))) 
#' sum( PolyShift(a,n)%*%coefs * (x -a)^(0:(n-1)))
#' @export
PolyShift <- function( x, n) {
  ret <- matrix(0,n,n)
  pow <- x^(col(ret) - row(ret))
  ret[col(ret)>=row(ret)] <- unlist( sapply(0:(n-1),function(x) choose(x,0:x)))
  ret*pow
}
