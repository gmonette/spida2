##
## Copy changes between spida2 and p3d
##
# \name{ell}
# \Rdversion{1.1}
# \alias{ell}
# %\alias{dell}
# \alias{ellplus}
# \alias{elltan}
# \alias{elltanc}
# \alias{ellpt}
# \alias{ellptc}
# \alias{ellbox}
# \alias{ellpts}
# \alias{ellptsc}
# %\alias{cell}
# %- Also NEED an '\alias' for EACH other topic documented here.
# \title{
#   Ellipses in 2D
# }
# \description{
#   Tools to generate 2D data, concentration and confidence ellipses given a center and a 'variance' matrix. Also tools
#   to generate points on a ellipse in a given direction or conjugate to that direction, axes
#   along a vector or conjugate to it, tangent lines at a point or parallel to a vector.
#   %%  ~~ A concise (1-5 lines) description of what the function does. ~~
# }
# \usage{
#   
#   ell( center = c(0,0), shape = diag(2), radius = 1, n = 100)
#   
#   %dell( x, y, radius = 1, n = 100)
#   
#   ellpt( ell, dir = c(0,1) , radius = 1 )
#   
#   ellptc( ell, dir = c(0,1) , radius = 1 )
#   
#   elltan( ell, dir = c(0,1) , radius = 1 , len = 1, v= c(-1,1))
#   
#   elltanc( ell, dir = c(0,1) , radius = 1 , len = 1, v= c(-1,1))
#   
#   ellbox( ell, dir = c(0,1) , radius = 1 )
#   
#   ellpts( ell, dir = c(0,1) , radius = 1 )
#   
#   ellptsc( ell, dir = c(0,1) , radius = 1 )
#   
# }
# %- maybe also 'usage' for other objects documented here.
# \arguments{
#   \item{center}{  the center of the ellipse
#   }
#   \item{shape}{   a 2 x 2 positive semi-definite. It is the variance matrix of a multivariate normal for which the ellipse is a concentration ellipse.
#   }
#   \item{radius}{ of the ellipse or point relative to the unit ellipse. For all functions except those generating tangents \code{radius} can be a vector.
#     For example, \code{ ellpt( ell, c(0,1), radius = c(-1,1))} will generate an axis of the unit ellipse.
#   }
#   \item{n}{  number of points to generate for the ellipse
#     %%     ~~Describe \code{n} here~~
#   }
#   \item{ell}{  is an object of class \code{ell} created by \code{ell}, \code{dell} or \code{cell}
#   }
#   \item{dir}{  vector giving the direction from the center of the ellipse to find a point on the ellipse or a tangent to it, or a direction conjugate to the direction in which ...
#   }
#   \item{len}{  half 'length' of a tangent vector
#   }
#   \item{v}{  tangent vector
#     %%     ~~Describe \code{n} here~~
#   }
# }
# \details{
#   \code{ell} returns a matrix of points on the ellipse, suitable for plotting
#   with \code{lines}.
#   
#   \code{ellpt} returns a point on an ellipse particular direction specified by \code{dir}
#   
#   \code{ellptc} returns a point on an ellipse in a conjugate direction specified by \code{dir}
#   
#   \code{ellpts} returs the nine points of the enclosing parallelogram (+ the centre) with the parallelogram whose
#   with a side parallel to \code{dir}. \code{ellptsc} returns the same points but in a different order
#   
# }
# \value{
#   The functions \code{ell} and \code{dell} return an object of class \code{ell} consisting of
#   matrix whose rows are points on the ellipse and, thus, can be plotted with \code{plot} or \code{lines}.
#   The other functions return a n x 2 matrix of points to plotted.
# }
# %\references{
#   %%% ~put references to the literature/web site here ~
#     %}
# \author{
#   Georges Monette
# }
# %\note{
#   %%%  ~~further notes~~
#     %}
# 
# %% ~Make other sections like Warning with \section{Warning }{....} ~
#   
#   \seealso{
#     \code{\link{cell}}, \code{\link{dell}}, \code{\link{ell.conj}}, ~~~
#   }
# \examples{
#   \dontrun{
#     ##---- Should be DIRECTLY executable !! ----
#     ##-- ==>  Define data, use random,
#     ##--	or do  help(data=index)  for the standard data sets.
#     
#     plot( e1 <- ell(c(1,1), matrix(c(1,.6,.6,1), ncol = 2)) )
#     ellptc ( e1, radius = c(-1,1))
#   }
# }
# 
# Modified by GM 2013-10-27:
#    added ellpts and ellptsc to generate 9 points
#    of tangent parallelogram
#
## Last modified by Georges Monette 2010-12-02
## Consolidated ell.R, dell.R, ell.conj.R, center.R and center.ell.F
## Added function ellpt, ellptc, elltan, elltanc, ellbox to generate points on ellipses,
## axes of ellipses and tangents.

#' Ellipses in 2D
#' 
#' Tools to generate 2D data, concentration and confidence 
#' ellipses given a center and a 'variance' matrix. Also tools
#' to generate points on a ellipse in a given direction or conjugate to that direction, axes
#' along a vector or conjugate to it, tangent lines at a 
#' point or parallel to a vector.
#' 
#' The ellipse is contour of the bivariate normal distribution
#' with variance 'shape'.
#' 
#' @param center (default: c(0,0)) 
#' @param shape variance of bivariate normal. Default: 2 x 2 identity
#' @param radius of ellipse, equivalently: square root of deviance contour
#' @param n number of points on ellipse
#' @param na.rm remove missing data in forming data ellipse, default TRUE
#' 
#' @export
ell <- 
  function( center = c(0,0), shape = diag(2) , radius  = 1, n = 100) {
    fac <- function( x )  {
      # fac(M) is a 'right factor' of a positive semidefinite M
      # i.e. M = t( fac(M) ) %*% fac(M)
      # similar to chol(M) but does not require M to be PD.
      xx <- svd(x)
      t(xx$v) * sqrt(pmax( xx$d,0))
    }
    angles = (0:n) * 2 * pi/n
    if ( length(radius) > 1) {
      ret <- lapply( radius, function(r) rbind(r*cbind( cos(angles), sin(angles)),NA))
      circle <- do.call( rbind, ret)
    }
    else circle = radius * cbind( cos(angles), sin(angles))
    ret <- t( c(center) + t( circle %*% fac(shape)))
    attr(ret,"parms") <- list ( center = rbind( center), shape = shape, radius = radius)
    class(ret) <- "ell"
    ret
  }

#' @export
dell <-
function( x, y, radius = 1, ..., na.rm = TRUE) {
        if ( (is.matrix(x) && (ncol(x) > 1))|| is.data.frame(x)) mat <- as.matrix(x[,1:2])
        else if (is.list(x)) mat <- cbind(x$x, x$y)
        else mat <- cbind( x,y)
        ell( apply(mat,2,mean,na.rm = na.rm), var(mat, na.rm = na.rm), radius = radius, ...)
    }

#' @export
ell.conj <-
function( center, shape, dir, radius = 1, len = 1) {
            # returns conjugate axes or tangent lines to ellipse
            vecs <- uv( shape, dir, radius)
            list( u = list( center-len*vecs$u, center+len*vecs$u),
                  v = list( center-len*vecs$v, center+len*vecs$v),
                  tan1 = list( center + vecs$u - len*vecs$v,center + vecs$u+ len*vecs$v),
                  tan2 = list( center - vecs$u - len*vecs$v,center - vecs$u+ len*vecs$v),
                  tan3 = list( center + vecs$v - len*vecs$u,center + vecs$v+ len*vecs$u),
                  tan4 = list( center - vecs$v - len*vecs$u,center - vecs$v+ len*vecs$u),
                  center = center)
    }


#' @export
center <-
function( obj, ... ) UseMethod("center")


#' @export
center.ell <-
function( obj, ...) attr(obj, 'parms') $ center



#' @export
ConjComp <- function( X , Z = diag( nrow(X)) , ip = diag( nrow(X)), tol = 1e-07 ) {
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

#' @export
uv <- function(object,...) UseMethod('uv')

#' @export
uv.ell <- function( object, u, radius = 1, ...){
        p <- attr(object,"parms")
        uv( p$shape, u=u, radius=radius)
}

#' @export
uv.default <-
function( object, u , radius = 1, ...) {
       # returns 'unit' u and conjugate v
            u <- u / sqrt( sum( u*solve(object,u)))   # 'unit' vector in direction of dir
            v <- c(ConjComp( u, diag(2) , solve(object)))  # conjugate
            v <- v / sqrt( sum( v * solve( object, v)))
       list(u = radius * u, v= radius * v)
    }



#' @export
ellpt <- function(ell, dir = c(0,1) , radius = 1 ) {
   # point on an ellipse in a particular direction
   p <- attr(ell,'parms')
   dir <- cbind(dir)
   dn <- sum(dir* solve(p$radius[1]^2 * p$shape, dir))
   ax <- dir / sqrt(dn)
   #disp(ax)
   #disp( rbind(radius))
   #disp( p$center)
   t( ax %*% rbind(radius) + as.vector(p$center))               # returns a row for plotting
}

#' @export
ellptc <- function(ell, dir = c(0,1) , radius = 1 ) {
   # point on an ellipse in a conjugate direction
   p <- attr(ell,'parms')
   dir <- cbind(dir)
   ax <- Null( solve(p$shape, dir) )
   dn <- sum(ax* solve( p$shape, ax))
   ax <- ax* p$radius[1] /sqrt(dn)
   t( ax %*% rbind(radius) + as.vector(p$center))               # returns a row for plotting
}

#' @export
elltanc <- function( ell, dir = c(0,1), radius = 1, len = 1, v = c(-1,1)) {
       p <- attr(ell,'parms')
       ax <- ellptc( ell, dir = dir, radius = len * v)
       if ( is.null(radius) ) return( t( t(ax) + dir))
       pt <- ellpt( ell, dir = dir, radius = radius)
       t( t(ax) -as.vector(p$center) + as.vector(pt))
}

#' @export
elltan <- function( ell, dir = c(0,1), radius = 1, len = 1, v = c(-1,1)) {
       p <- attr(ell,'parms')
       ax <- ellpt( ell, dir = dir, radius =  len * v)
       pt <- ellptc( ell, dir = dir, radius = radius)
       t( t(ax) -as.vector(p$center) + as.vector(pt))
}

#' @export
ellbox <- function( ell, dir = c(0,1) , radius = 1 ){
      rbind(
 elltan(  ell, dir = dir, radius = radius) , NA,
 elltanc( ell, dir = dir, radius = radius) , NA,
 elltan(  ell, dir  = dir, radius = -radius), NA,
 elltanc( ell, dir = dir, radius = -radius) )
}

#' @export
ellplus <-
function (center = rep(0, 2), shape = diag(2), radius = 1, n = 100,
    angles = (0:n) * 2 * pi/n, fac = chol, ellipse = all, diameters = all,
    box = all, all = FALSE)
{
    help <- "\n        ellplus can produce, in addition to the points of an ellipse, the\n        conjugate axes corresponding to a chol or other decomposition\n        and the surrounding parallelogram.\n        "
    rbindna <- function(x, ...) {
        if (nargs() == 0)
            return(NULL)
        if (nargs() == 1)
            return(x)
        rbind(x, NA, rbindna(...))
    }
    if (missing(ellipse) && missing(diameters) && missing(box))
        all <- TRUE
    circle <- function(angle) cbind(cos(angle), sin(angle))
    Tr <- fac(shape)
    ret <- list(t(c(center) + t(radius * circle(angles) %*% Tr)),
        t(c(center) + t(radius * circle(c(0, pi)) %*% Tr)), t(c(center) +
            t(radius * circle(c(pi/2, 3 * pi/2)) %*% Tr)), t(c(center) +
            t(radius * rbind(c(1, 1), c(-1, 1), c(-1, -1), c(1,
                -1), c(1, 1)) %*% Tr)))
    do.call("rbindna", ret[c(ellipse, diameters, diameters, box)])
}

#' @export
ellpts <- function( ell, dir = c(0,1), radius = 1, len = 1, v = c(-1,1)) {
  rbind(
    ellpt( ell, dir = dir, radius = radius * c(-1,0,1)),
    elltan( ell, dir = dir, radius = radius,
            v = radius * c(-1,0,1)),
    elltan( ell, dir = dir, radius = -radius,
            v = radius * c(-1,0,1)))
}

#' @export
ellptsc <- function( ell, dir = c(0,1), radius = 1, len = 1, v = c(-1,1)) {
  rbind(
    ellptc( ell, dir = dir, radius = radius * c(-1,0,1)),
    elltanc( ell, dir = dir, radius = radius,
            v = radius * c(-1,0,1)),
    elltanc( ell, dir = dir, radius = -radius,
            v = radius * c(-1,0,1)))
}

#' @export
Null <- function (M) 
{
  # from MASS::Null
  tmp <- qr(M)
  set <- if (tmp$rank == 0L) 
    seq_len(ncol(M))
  else -seq_len(tmp$rank)
  qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
}
#' Generate line of osculation between two families of ellipses
#' 
#' @export
osculant <- function(x, ...) {
  UseMethod('osculant')
}
#' @param center1 center of first ellipse
#' @param shape1 shape (variance) of first ellipse
#' @param center2 center of second ellipse
#' @param shape2 shape (variance) of second ellipse
#' @param n number of points to generate, default 100
#' @examples
#' sq <- function(A) {
#'   ev <- eigen(A)
#'   ev$vectors %*% diag(sqrt(ev$values)) %*% t(ev$vectors) 
#' }
#' sqinv <- function(A) {
#'   ev <- eigen(A)
#'   ev$vectors %*% diag(1/sqrt(ev$values)) %*% t(ev$vectors) 
#' }
#' rot <- function(theta, ipm = diag(2)) {
#'   sq(ipm) %*% rbind(c(cos(theta), -sin(theta)), c(sin(theta), cos(theta))) %*% sqinv(ipm)
#' }
#' rotm <- function(mat, theta, ipm = diag(2)){
#'   rot(theta, ipm = ipm) %*% mat %*% rot(-theta, ipm = ipm)
#' }
#' cen1 <- c(0,0)
#' shape1 <- rotm(diag(c(4,1)), pi/4)
#' cen2 <- c(2,2)
#' shape1 <- rotm(diag(c(6,1)), - pi/4)
#' osculant(cen1, shape1, cen2, shape2)
#' 
#' plot(rbind(10*c(1,1),-c(5,5)), asp = 1, xlab = '', ylab = '', type = 'n')
#' (ell1 <- ell(c1 <- c(0,0), G1 <- rotm(diag(c(16,1)), -3*pi/4))) %>% lines(col = 'red')
#' (ell2 <- ell(c2 <- c(6,3), G2 <- rotm(diag(c(16,1)), pi/2))) %>% lines(col = 'red')
#' osculant(ell1, ell2) %>% lines
#' 
#' # (rbind(0, c2) %*% rot(pi/2)) %>% lines
#' rbind(0, c2)  %>% lines
#' # -rbind(0, c2) %*% rot(-pi/2, ipm = solve(G1)) %>% lines
#' rbind(0, c2) %*% solve(G2) %*% G1 %>% lines                       # ok
#' -rbind(c2,0) %*% solve(G1) %*% G2 %>% {. + rbind(c2,c2)} %>% lines # ok
#' -rbind(c2,0) %*% rot(pi/2, ipm = G2) %>% lines
#' -rbind(c2,0) %*% rot(pi/2, ipm = G2) %*% rot(pi/2) %>% lines(col = 'red')
#' -rbind(c2,0) %*% rot(pi/2, ipm = G2) %*% rot(pi/2) %*% rot(pi/2, ipm = solve(G1))%>% lines(col = 'red')
#'
#' 
#' rbind(c1, c2) %>% lines
#' rbind(c1, c2) %*% rot(pi/2) %>% lines
#' rbind(c1, c2) %*% G1 %>% lines
#' rbind(c1, c2) %*% solve(G1) %>% lines(col='blue')

#' ell(cen1, shape1) %>% lines
#' ell(cen2, shape2) %>% lines
#' @rdname osculant 
#' @export
osculant.default <- function(
    center1,
    shape1,
    center2,
    shape2,
    n = 100) {
  w1 <- solve(shape1)
  w2 <- solve(shape2)
  a1 <- w1 %*% center1
  a2 <- w2 %*% center2
  ret <- sapply(seq(0,1, length.out=n),
         function(k){
           solve(k * w1 + (1-k) * w2)%*%(k*a1 + (1-k)*a2)
         }
           )
  t(ret)
}
#' @param ell1 first ellipse, an object of class `ell` created by [ell()]
#' @param ell2 second ellipse
#' @param n number of points to generate, default 100
#' @examples
#' c1 <- c(-5, -1)
#' c2 <- c(5, 0)
#' ell1 <- ell(c1,  cbind(c(3,3.1),c(3.1,4) ))
#' ell2 <- ell(c2,  cbind(c(3,-1),c(-1,4) ))
#' ell3 <- ell(c(5,0),  cbind(c(3,-.1),c(-.1,.1) ))
#' ell4 <- ell(c(5,0),  cbind(c(3,0),c(0,.1) ))
#' ell5 <- ell(c(5,0),  cbind(c(3,0.01),c(0.01,.01) ))
#' 
#' osculant(ell1,  ell2) 
#' perp <- function(x) {
#'   # rotate a line segment 90 deg
#'   x <- rbind(x,0)[1:2,]
#'   v <- c(1,-1) %*% x
#'   vp <- v %*% rbind(c(0,-1),c(1,0))
#'   rbind(vp + x[2, ], x[2,] )
#' }
#'   
#' }UseMethod('perp')
#' perp.vec <- function(x,...){
#'    
#' }
#' 
#' }
#' lseg(c(1,1),c(2,3)) %>% vec
#' plot(c(-5,5), c(-5,5), asp = 1, type = 'n', xlab ='', ylab='')
#' lines(ell1)
#' lines(ell2)
#' osculant(ell1,ell2) %>% lines
#' lines(rbind(c1,c2)) 
#' vec1 <- elltanc(ell2,c1-c2)
#' lines(vec1, col = 'red')
#' vec1
#' assign('vec1',.) %>% lines(col = 'red')
#' perp_dir <-  
#' vec1
#' # ellptsc(ell1, dir = c(1,0)) %>% points(pch = 16, cex = .5)
#' # elltanc(ell2,c(1,0)) %>% lines
#' # elltanc(ell1,c(1,0)) %>% lines
#' elltanc(ell2,elltanc(ell1,c(1,0))[2,]) %>% lines
#' plot(c(-5,5), c(-5,5), asp = 1)
#' lines(ell1) 
#' lines(ell3)
#' osculant(ell1,ell3) %>% lines
#' lines(ell4, col = 'red')
#' osculant(ell1,ell4) %>% lines(col='red')
#' lines(ell5, col = 'blue')
#' osculant(ell1,ell5) %>% lines(col='blue')
#' s11 <- 3.0
#' s12 <- +0.4
#' rotv <- function(v, theta) {
#'   rot <- cbind(c(cos(theta),sin(theta)),c(-sin(theta), cos(theta)))
#'   rot %*% v %*% t(rot)
#' }
#' {
#'   for(theta1 in pi*seq(0:10)/10){
#'     plot(c(-5,5),2* c(-5,5), asp = 1, type = 'n',
#'       xlab = '',ylab = '', sub = 'osculating ellipses')
#'     ell1 <- ell(c(-5,0),  rotv(cbind(c(3,3.1),c(3.1,4) ), theta1), radius =3)
#'     lines(ell1)
#'     for(theta in pi*seq(0:100)/100){
#'       ell6 <- ell(c(10,0),  rotv(cbind(c(s11,s12),c(s12, s12^2/s11 + .1) ), theta))
#'       lines(ell6, col = 'red')
#'       osculant(ell1,ell6,500) %>% lines(col='red')
#'     }
#'     Sys.sleep(1)
#'   }
#' }
#' @rdname osculant
#' @export
osculant.ell <- function(
    ell1,
    ell2,
    n = 100) {
  cen1 <- center(ell1)
  shape1 <- attr(ell1,'parms')$radius * attr(ell1,'parms')$shape
  cen2 <- center(ell2)
  shape2 <- attr(ell2,'parms')$radius * attr(ell2,'parms')$shape
  osculant(as.vector(cen1),shape1,as.vector(cen2),shape2,n=n)
}

