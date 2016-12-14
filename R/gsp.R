##
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
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param X %% ~~Describe \code{X} here~~
#' @param coef %% ~~Describe \code{coef} here~~
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
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param knots %% ~~Describe \code{knots} here~~
#' @param degree %% ~~Describe \code{degree} here~~
#' @param smooth %% ~~Describe \code{smooth} here~~
#' @param intercept %% ~~Describe \code{intercept} here~~
#' @param signif %% ~~Describe \code{signif} here~~
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
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param knots %% ~~Describe \code{knots} here~~
#' @param degree %% ~~Describe \code{degree} here~~
#' @param smooth %% ~~Describe \code{smooth} here~~
#' @param intercept %% ~~Describe \code{intercept} here~~
#' @param signif %% ~~Describe \code{signif} here~~
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





#' General regression splines with variable degrees and smoothnes, smoothing
#' splines
#' 
#' These functions implement a general spline with possibly different degrees
#' in each interval and and different orders of smoothness at each knot,
#' including the possibility of allowing a discontinuity at a knot. The
#' function \code{sc} helps in the construction of linear hypothesis matrices
#' to estimate and test levels and derivatives of splines at arbitrary points
#' and the saltus of derivatives that have discontinuities at knots.
#' 
#' Many polynomial regression splines can be generated by 'plus' functions
#' although the resulting basis for the spline may be ill conditioned. For
#' example a 'quadratic spline' (a spline that is quadratic in each interval
#' with and smooth with a first derivative at each knot) with knots at 1 and 3
#' can be fitted with: \preformatted{ plus <- function(x,y) ifelse(x>0,y,0) lm(
#' y ~ x + I(x^2) + plus(x-1,(x-1)^2) + plus(x-3,(x-3)^2)) } All 'standard'
#' polynomial splines with the same degree in each interval and continuity of
#' order one less than the degree at each knot can be constructed in this
#' fashion. A convenient aspect of this parametrization of the spline is that
#' the estimated coefficients have a simple interpretation. The coefficients
#' for 'plus' terms represent the 'saltus' (jump) in the value of a coefficient
#' at the knot. Testing whether the true value of a coefficient is 0 is
#' equivalent to a test for the need to retain the corresponding knot. \cr\cr
#' This approach does not work for some more complex splines with different
#' degrees or different orders of continuity at the knots. An example is the
#' commonly used natural quadratic spline. A natural quadratic spline with
#' knots at -1,0 and 1 (where -1 and 1 are termed 'boundary knots') is linear
#' in the intervals (-Inf,-1) and (1,+Inf), and quadratic in the intervals
#' (-1,0) and (0,1). The spline is smooth of order 1 at each knot.  \cr\cr Many
#' techniques for fitting splines generate a basis for the spline (columns of
#' the design matrix) that has good numerical properties but whose coefficients
#' do not have a simple interpretation. \cr\cr The \code{gsp} function makes it
#' easy to specify a spline with arbitrary degree in each interval and
#' arbitrary smoothness at each knot. The parametrization produces coefficients
#' that have a simple interpretation. For a spline of degree p at x = 0,
#' coefficients correspond to the 1st, 2nd, ... pth derivative at 0. Additional
#' coefficients correspond to each free saltus at each knot. \cr\cr The
#' \code{sc} function generates a matrix to estimate features of a fitted
#' spline that can be expressed as linear combinations of the spline
#' coefficients.  Examples are various derivatives of the spline at any point,
#' left or right derivatives of different orders and the saltus in derivatives
#' at a knot.  \cr\cr A disadvantage of \code{gsp} is that the spline basis may
#' be poorly conditioned.  The impact of this problem can be mitigated by
#' rescaling the x variable so that it has an appropriate range. For example,
#' with a spline whose highest degree is cubic, ensuring that x has a range
#' between -10 and 10 should avert numerical problems. \cr\cr \code{gsp}
#' generates a matrix of regressors for a spline with knots, degree of
#' polynomials in each interval and the degree of smoothness at each knot.
#' Typically, \code{gsp} is used to define a function that is then used in a
#' model equation. See the examples below.
#' 
#' A function to fit a cubic spline with knots at 5 and 10 is generated with:
#' 
#' \preformatted{ sp <- function( x ) gsp( x, c(5,10), c(3,3,3), c(2,2)) }
#' 
#' indicating that a cubic polynomial is used in each of the three intervals
#' and that the second derivative is continuous at each knot.
#' 
#' A 'natural cubic spline' with linear components in each unbounded interval
#' would have the form:
#' 
#' \preformatted{ sp <- function( x ) gsp( x, c(0,5,10,15), c(1,3,3,3,1),
#' c(2,2,2,2)) }
#' 
#' Quadratic and linear splines, respectively: \preformatted{ sp.quad <-
#' function( x ) gsp( x, c(5,10), c(2,2,2), c(1,1)) sp.lin <- function( x )
#' gsp( x, c(5,10), c(1,1,1), c(0,0)) } Where the same degree is used for all
#' intervals and knots, it suffices to give it once: \preformatted{ sp.quad <-
#' function( x ) gsp( x, c(5,10), 2, 1) sp.lin <- function( x ) gsp( x,
#' c(5,10), 1, 0) } An easy way to specify a model in which a knot is dropped
#' is to force a degree of continuity equal to the degree of adjoining
#' polynomials, e.g.  to drop the knot at 10, use: \preformatted{ sp.1 <-
#' function( x ) gsp( x, c(5,10), c(3,3,3), c(2,3)) } This is sometimes easier
#' than laboriously rewriting the spline function for each null hypothesis.
#' 
#' Depending on the maximal degree of the spline, the range of x should not be
#' excessive to avoid numerical problems. The spline matrix generated is 'raw'
#' and values of max(abs(x))^max(degree) may appear in the matrix.  For
#' example, for a cubic spline, it might be desirable to rescale x and/or
#' recenter x so abs(x) < 100 if that is not already the case. Note that the
#' knots need to be correspondingly transformed.
#' 
#' The naming of coefficients should allow them to be easily interpreted.  For
#' example: \preformatted{ > zapsmall( gsp ( 0:10, c(3, 7) , c(2,3,2), c(1,1))
#' )
#' 
#' D1(0) D2(0) C(3).2 C(3).3 C(7).2 f(0) 0 0.0 0.0 0.00000 0.0 f(1) 1 0.5 0.0
#' 0.00000 0.0 f(2) 2 2.0 0.0 0.00000 0.0 f(3) 3 4.5 0.0 0.00000 0.0 f(4) 4 8.0
#' 0.5 0.16667 0.0 f(5) 5 12.5 2.0 1.33333 0.0 f(6) 6 18.0 4.5 4.50000 0.0 f(7)
#' 7 24.5 8.0 10.66667 0.0 f(8) 8 32.0 12.5 20.66667 0.5 f(9) 9 40.5 18.0
#' 34.66667 2.0 f(10) 10 50.0 24.5 52.66667 4.5 }
#' 
#' The coefficient for the first regressor is the first derivative at x = 0;
#' for the second regressor, the second derivative at 0; the third, the saltus
#' (change) in the second derivative at x = 3, the fourth, the saltus in the
#' third derivative at x = 3 and, finally, the saltus in the second derivative
#' at x = 7.
#' 
#' Example: \preformatted{ > sp <- function(x) gsp ( x, c(3, 7) , c(2,3,2),
#' c(1,1)) > zd <- data.frame( x = seq(0,10, .5), y = seq(0,10,.5)^2 + rnorm(
#' 21)) > fit <- lm( y ~ sp( x ), zd) > summary(fit) > Ls <-cbind( 0, sc( sp,
#' c(1,2,3,3,3,5,7,7,7,8), D=3, + type = c(0,0,0,1,2,0,0,1,2,0))) > zapsmall(
#' Ls )
#' 
#' D1(0) D2(0) C(3).2 C(3).3 C(7).2 D3(1) 0 0 0 0 0 0 D3(2) 0 0 0 0 0 0 D3(3-)
#' 0 0 0 0 0 0 D3(3+) 0 0 0 0 1 0 D3(3+)-D3(3-) 0 0 0 0 1 0 D3(5) 0 0 0 0 1 0
#' D3(7-) 0 0 0 0 1 0 D3(7+) 0 0 0 0 0 0 D3(7+)-D3(7-) 0 0 0 0 -1 0 D3(8) 0 0 0
#' 0 0 0
#' 
#' > wald( fit, list( 'third derivatives' = Ls))
#' 
#' numDF denDF F.value p.value third derivatives 1 15 2.013582 0.17634
#' 
#' Estimate Std.Error DF t-value p-value Lower 0.95 Upper 0.95 D3(1) 0.000000
#' 0.000000 15 -1.123777 0.27877 0.000000 0.000000 D3(2) 0.000000 0.000000 15
#' -1.123777 0.27877 0.000000 0.000000 D3(3-) 0.000000 0.000000 15 -1.123777
#' 0.27877 0.000000 0.000000 D3(3+) 0.927625 0.653714 15 1.419008 0.17634
#' -0.465734 2.320984 D3(3+)-D3(3-) 0.927625 0.653714 15 1.419008 0.17634
#' -0.465734 2.320984 D3(5) 0.927625 0.653714 15 1.419008 0.17634 -0.465734
#' 2.320984 D3(7-) 0.927625 0.653714 15 1.419008 0.17634 -0.465734 2.320984
#' D3(7+) 0.000000 0.000000 Inf NaN NaN 0.000000 0.000000 D3(7+)-D3(7-)
#' -0.927625 0.653714 15 -1.419008 0.17634 -2.320984 0.465734 D3(8) 0.000000
#' 0.000000 Inf NaN NaN 0.000000 0.000000
#' 
#' Warning messages: 1: In min(dfs[x != 0]) : no non-missing arguments to min;
#' returning Inf 2: In min(dfs[x != 0]) : no non-missing arguments to min;
#' returning Inf
#' 
#' } Note that some coefficients that are 0 by design may lead to invalid DRs
#' and t-values.
#' 
#' \code{sc} generates a portion of a hypothesis matrix for the coefficients of
#' a general spline constructed with \code{gsp} With: \preformatted{ sc( sp, x,
#' D, type ):
#' 
#' sp is the spline function for which coefficients are required x values at
#' which spline is evaluated D order of derivative: 0 = value of spline, 1 =
#' first derivative, etc.  type at knots: 0 limit from the left, 1 from the
#' right, 2 is saltus (i.e. jump from left to right) }
#' 
#' Warning: \code{sc} will not work correctly if the function defining the
#' spline transforms the variable, e.g. sp <- function(x) gsp( x/100, knot=2 )
#' 
#' Example:
#' 
#' \preformatted{ simd <- data.frame( age = rep(1:50, 2), y =
#' sin(2*pi*(1:100)/5)+ rnorm(100), G = rep( c("male","female"), c(50,50)))
#' 
#' sp <- function(x) gsp( x, knots = c(10,25,40), degree = c(1,2,2,1), smooth =
#' c(1,1,1))
#' 
#' fit <- lm( y ~ sp(age)*G, simd) xyplot( predict(fit) ~ age , simd, groups =
#' G,type = "l") summary(fit) # convenient display
#' 
#' # output:
#' 
#' Call: lm(formula = y ~ sp(age) * G, data = simd)
#' 
#' Residuals: Min 1Q Median 3Q Max -2.5249 -0.7765 -0.0760 0.7882 2.6265
#' 
#' Coefficients: Estimate Std. Error t value Pr(>|t|) (Intercept) 0.733267
#' 0.605086 1.212 0.229 sp(age)D1(0) -0.084219 0.055163 -1.527 0.130
#' sp(age)C(10).2 0.010984 0.006910 1.590 0.115 sp(age)C(25).2 -0.023034
#' 0.012881 -1.788 0.077 . Gmale -0.307665 0.855721 -0.360 0.720
#' sp(age)D1(0):Gmale 0.058384 0.078012 0.748 0.456 sp(age)C(10).2:Gmale
#' -0.010556 0.009773 -1.080 0.283 sp(age)C(25).2:Gmale 0.026410 0.018216 1.450
#' 0.150 ---
#' 
#' Residual standard error: 1.224 on 92 degrees of freedom Multiple R-squared:
#' 0.0814, Adjusted R-squared: 0.0115 F-statistic: 1.165 on 7 and 92 DF,
#' p-value: 0.3308 # end of output
#' 
#' L0 <- list( "hat" = rbind( "females at age=20" = c( 1, sc(sp,20), 0, 0*
#' sc(sp,20)), "males at age=20" = c( 1, sc(sp,20), 1, 1* sc(sp,20))),
#' "male-female" = rbind( "at 20" = c( 0 , 0*sc(sp,20), 1, 1*sc(sp,20)))) wald(
#' fit, L0 )
#' 
#' ...
#' 
#' L1 <- list("D(yhat)/D(age)"= rbind( "female at age = 25" = c(0, sc(sp,25,1),
#' 0, 0*sc(sp,25,1)), "male at x = 25" = c(0, sc(sp,25,1), 0, 1*sc(sp,25,1))))
#' wald( fit, L1) # output: numDF denDF F.value p.value D(yhat)/D(age) 2 92
#' 1.057307 0.35157
#' 
#' Estimate Std.Error DF t-value p-value Lower 0.95 Upper 0.95 female at age =
#' 25 0.080544 0.056974 92 1.413694 0.16083 -0.032612 0.193700 male at x = 25
#' -0.019412 0.056974 92 -0.340712 0.73410 -0.132568 0.093744 }
#' 
#' A periodic spline function can be generated by forcing the coefficients
#' beyond a periodic knot to repeat the pattern in a previous inteval. For
#' example a periodic spline of period 1 can be created as follows:
#' 
#' \preformatted{ x <- seq(-3,3,.01) y <- rnorm(x)
#' 
#' per.sp <- function(x) gsp( x %% 1, knots = 1, degree = c(3, 3), smooth = 1,
#' lin = cbind( diag(4), - PolyShift(1,4))) fit <- lm( y ~ per.sp(x))
#' summary(fit) xyplot( y + predict(fit) ~ x , panel = panel.superpose.2,type =
#' c("p","l")) xyplot( predict(fit) ~ x , type = 'l') }
#' 
#' A periodic spline with additional knots can be created as shown below. Note
#' that constraint matrix 'lin' expresses constraints for the 'full' polynomial
#' parametrization, i.e. polynomials of degree 3 (thus of order 4 when
#' including the constant term) in each of the 4 intervals.  The constraint
#' given in 'lin' forces the coefficients beyond the periodic knot at x = 1, to
#' repeat the polynomial in the interval just to the right of x = 0.
#' 
#' \preformatted{ x <- seq(-3,3,.01) y <- rnorm(x)
#' 
#' per.sp <- function(x) gsp( x %% 1, knots = c(.333,.666,1), degree = 3,
#' smooth = 1, lin = cbind( diag(4),0*diag(4),0*diag(4), - PolyShift(1,4))) fit
#' <- lm( y ~ per.sp(x)) summary(fit) xyplot( y + predict(fit) ~ x , panel =
#' panel.superpose.2,type = c("p","l")) xyplot( predict(fit) ~ x , type = 'l')
#' }
#' 
#' Overview of utility functions: \preformatted{ Xmat = function( x, degree, D
#' = 0, signif = 3) design/estimation matrix for f[D](x) where f(x) is
#' polynomial of degree degree.
#' 
#' Xf = function( x, knots, degree = 3, D = 0, right = TRUE , signif = 3) uses
#' Xmat to form 'full' matrix with blocks determined by knots intervals
#' 
#' Cmat = function( knots, degree, smooth, intercept = 0, signif = 3) linear
#' constraints
#' 
#' Emat = function( knots, degree, smooth , intercept = FALSE, signif = 3)
#' estimates - not necessarily a basis
#' 
#' basis = function( X , coef = FALSE ) selects linear independent subset of
#' columns of X
#' 
#' spline.T = function( knots, degree, smooth, intercept = 0, signif = 3 ) full
#' transformation of Xf to spline basis and constraints
#' 
#' spline.E = function( knots, degree, smooth, intercept = 0, signif = 3 )
#' transformation for spline basis (first r columns of spline.T)
#' 
#' gsp = function( x , knots, degree = 3 , smooth = pmax(pmin( degree[-1],
#' degree[ - length(degree)]) - 1,0 ), intercept = 0, signif = 3)
#' 
#' }
#' 
#' @aliases gsp sc smsp Xf Xmat qs lsp Cmat PolyShift
#' @param x value(s) where spline is evaluated
#' @param knots vector of knots
#' @param degree vector giving the degree of the spline in each interval. Note
#' the number of intervals is equal to the number of knots + 1. A value of 0
#' corresponds to a constant in the interval. If the spline should evaluate to
#' 0 in the interval, use the \code{intercept} argument to specify some value
#' in the interval at which the spline must evaluate to 0.
#' @param smooth vector with the degree of smoothness at each knot (0 =
#' continuity, 1 = smooth with continuous first derivative, 2 = continuous
#' second derivative, etc. The value -1 allows a discontinuity at the knot.
#' @param intercept value(s) of x at which the spline has value 0, i.e. the
#' value(s) of x for which yhat is estimated by the intercept term in the
#' model. The default is 0. If NULL, the spline is not constrained to evaluate
#' to 0 for any x.
#' @param periodic if TRUE generates a period spline on the base interval
#' [0,max(knots)]. A constraint is generated so that the coefficients generate
#' the same values to the right of max(knots) as they do to the right of 0.
#' Note that all knots should be strictly positive.
#' @param lin provides a matrix specifying additional linear contraints on the
#' 'full' parametrization consisting of blocks of polynomials of degree equal
#' to max(degree) in each of the length(knots)+1 intervals of the spline. See
#' below for examples of a spline that is 0 outside of its boundary knots.
#' @param signif number of significant digits used to label coefficients
#' @param exclude number of leading columns to drop from spline matrix: 0:
#' excludes the intercept column, 1: excludes the linear term as well.  Terms
#' that are excluded from the spline matrix can be modeled explicitly.
#' @param sp a spline function defined by \code{gsp}. See the examples below.
#' @param D the degree of a derivative: 0: value of the function, 1: first
#' derivative, 2: second derivative, etc.
#' @param type how a derivative or value of a function is measured at a
#' possible discontinuity at a knot: 0: limit from the left, 1: limit from the
#' right, 2: saltus (limit from the right minus the limit from the left)
#' @param a a 'periodic' knot at which a spline repeats itself
#' @param n the maximal 'order', i.e. maximal degree + 1, of a periodic spline
#' @return \code{gsp} returns a matrix generating a spline. \code{cs},
#' \code{qs} and \code{lsp} return matrices generating cubic, quadratic and
#' linear splines respectively.
#' 
#' \code{smsp}, whose code is adapted from function in the package
#' \code{lmeSplines}, generates a smoothing spline to be used in the random
#' effects portion of a call to \code{lme}.
#' @note %% ~~further notes~~
#' @section Warning: The variables generated by \code{gsp} are designed so the
#' coefficients are interpretable as changes in derivatives at knots. The
#' resulting matrix is not designed to have optimal numerical properties.
#' 
#' The intermediate matrices generated by \code{gsp} will contain \code{x}
#' raised to the power equal to the highest degree in \code{degree}. The values
#' of \code{x} should be scaled to avoid excessively large values in the spline
#' matrix and ensuing numerical problems.
#' @author Monette, G. \email{georges@@yorku.ca}
#' @seealso \code{\link{wald}}
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ## Fitting a quadratic spline
#' simd <- data.frame( age = rep(1:50, 2), y = sin(2*pi*(1:100)/5) + rnorm(100),
#'           G = rep( c('male','female'), c(50,50)))
#' # define a function generating the spline
#' sp <- function(x) gsp( x, knots = c(10,25,40), degree = c(1,2,2,1),
#'           smooth = c(1,1,1))
#' 
#' fit <- lm( y ~ sp(age)*G, simd)
#' 
#' require(lattice)
#' xyplot( predict(fit) ~ age , simd, groups = G,type = 'l')
#' summary(fit)
#' 
#' ## Linear hypotheses
#' L <- list( "Overall test of slopes at 20" = rbind(
#'       "Female slope at age 20" =  c( F20 <- cbind( 0 , sc( sp, 20, D = 1), 0 , 0 * sc( sp, 20, D = 1))),
#'       "Male slope at age 20" =  c( M20 <- cbind( 0 , sc( sp, 20, D = 1), 0 , 1 * sc( sp, 20, D = 1))),
#'       "Difference" = c(M20 - F20))
#'       )
#' wald( fit, L)
#' 
#' ## Right and left second derivatives at knots and saltus
#' L <- list( "Second derivatives and saltus for female curve at knot at 25" =
#'           cbind( 0, sc( sp, c(25,25,25), D = 2, type =c(0,1,2)), 0,0,0,0))
#' L
#' wald( fit, L )
#' 
#' ## Smoothing splines
#'       library(nlme)
#'       data(Spruce)
#'       Spruce$all <- 1
#'       range( Spruce$days)
#'       sp <- function(x) smsp ( x, seq( 150, 675, 5))
#'       spruce.fit1 <- lme(logSize ~ days, data=Spruce,
#'           random=list(all= pdIdent(~sp(days) -1),
#'           plot=~1, Tree=~1))
#'       summary(spruce.fit1)
#'       pred <- expand.grid( days = seq( 160, 670, 10), all = 1)
#'       pred$logSize <- predict( spruce.fit1, newdata = pred, level = 1)
#'       require( lattice )
#'       xyplot( logSize ~ days, pred, type = 'l')
#' 
#' @export
gsp <- function( x , knots, degree = 3 , smooth = pmax(pmin( degree[-1],
                                                             degree[ - length(degree)]) - 1,0 ), lin = NULL, periodic = FALSE, intercept = 0, signif = 3) {
  
  if ( periodic  ) {
    maxd <- max(degree)
    #disp(maxd)
    maxk <- max(knots)
    #disp(maxk)
    mid <- matrix(0,maxd+1,(maxd+1)*(length(knots)-1))
    const.per <- do.call(cbind,list( diag(maxd+1),
                                     mid, -PolyShift(maxk,maxd+1)))
    #disp(const.per)
    lin <- rbind(lin,const.per)
    #disp(lin)
  }
  degree = rep( degree, length = length(knots) + 1)
  smooth = rep( smooth, length = length(knots))
  spline.attr <-   list( knots = knots, degree = degree,
                         smooth = smooth, lin = lin, intercept = intercept, signif = signif)
  if (is.null(x)) return ( spline.attr)
  if( periodic ) x <- x %% maxk
  ret = Xf ( x, knots, max(degree), signif = signif) %*%
    spline.E( knots, degree, smooth, lin = lin, intercept = intercept, signif = signif)
  attr(ret, "spline.attr") <- spline.attr
  ret
}




#' Projection matrix %% ~~function to do ... ~~
#' 
#' Projection matrix into column space of a matrix. Works even if the matrix
#' does not have full column rank.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param x %% ~~Describe \code{x} here~~
#' @export
Proj.1 <- function( x)  x %*% ginv(x)



#' Projection matrix
#' 
#' Matrix of orthogonal projection into column space of matrix.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param x %% ~~Describe \code{x} here~~
#' @export
Proj <- function(x) {
  # projection matrix onto span(x)
  u <- svd(x,nv=0)$u
  u %*% t(u)
  
}



#' Test accuracy of projection method.
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param x %% ~~Describe \code{x} here~~
#' @param fun %% ~~Describe \code{fun} here~~
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
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param x %% ~~Describe \code{x} here~~
#' @param knots %% ~~Describe \code{knots} here~~
#' @param degree %% ~~Describe \code{degree} here~~
#' @param smooth %% ~~Describe \code{smooth} here~~
#' @param intercept %% ~~Describe \code{intercept} here~~
#' @param signif %% ~~Describe \code{signif} here~~
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

#' @export
sc <-
  function( sp, x , D = 0, type = 1) {
    
    a = sp(NULL)
    D = rep(D, length.out = length(x))
    type = rep ( type, length.out = length(x))
    left = Xf( x, knots = a$knots, degree = max( a$degree), D = D, right = T)
    right = Xf( x, knots = a$knots, degree = max( a$degree), D = D, right = F)
    cleft = c(1,0,-1) [ match( type, c(0,1,2)) ]
    cright = c(0,1,1) [ match( type, c(0,1,2)) ]
    raw = left * cleft + right * cright
    nam = rownames( raw )
    nam = sub("^f","g", nam)
    nam0 = sub( "\\)","-)", nam)
    nam1 = sub( "\\)","+)", nam)
    nam2 = paste( nam1 ,"-",nam0,sep = '')
    rownames( raw ) =
      ifelse( match( x , a$knots, 0) > 0,
              cbind( nam0,nam1,nam2) [ cbind( seq_along(type), type+1)],
              ifelse( type != 2, nam, '0'))
    mod = raw %*% spline.E(
      a$knots, a$degree, a$smooth,
      lin = a$lin,
      intercept = a$intercept,
      signif = a$signif)
    mod
  }




#' Generate a smoothing spline
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param formula %% ~~Describe \code{formula} here~~
#' @param data %% ~~Describe \code{data} here~~
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
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param time %% ~~Describe \code{time} here~~
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
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param Z %% ~~Describe \code{Z} here~~
#' @param oldtimes %% ~~Describe \code{oldtimes} here~~
#' @param newtimes %% ~~Describe \code{newtimes} here~~
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



#' @export
PolyShift <- function( x, n) {
  # transformation of polynomial coefficients to change origin 
  # sum( coefs * x ^ 0:n) == sum( Pow(a,n)%*%coefs * (x -a)^0:n
  # Useful to equate polynomial coefficients for a periodic spline
  ret <- matrix(0,n,n)
  pow <- x^(col(ret) - row(ret))
  ret[col(ret)>=row(ret)] <- unlist( sapply(0:(n-1),function(x) choose(x,0:x)))
  ret*pow
}