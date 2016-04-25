#' Create 'derivatives' and 'means' of factors to generate, for example,
#' pairwise differences or centres of existing factor for prediction
#'
#' The functions \code{xlevels} and \code{dlevels} are primarily intended to
#' create arguments for \code{expand.grid} to create 'prediction' or 'effect
#' data frames', to generate wald tests and estimates of specific effects in
#' models with interactions,
#'
#' @param f a factor or otherwise a vector interpreted as the levels of a
#' factor
#' @param type one or more types of contrasts ('derivative') or means (convex
#' combinations) of factor levels
#' @param all for some values of \code{type}, indicates whether to use all
#' contrasts. e.g. \code{type = "pairwise"} will produce pairwise comparisons
#' in both directions if \code{all == TRUE}
#' @param sep can add spaces in constructed factor levels
#' @param w weights used for each factor level in creating contrasts
#' differentiating a factor level from others. Weights corresponding to
#' frequencies in the data frame result in effects corresponding to type II
#' effects while equal weights correspond to type III effects for interacting
#' specific effects
#' @export
xlevels <- function(f,   type = c("raw","<mean>"),
                    all = FALSE, sep = '') {
  # produces the equivalent of differentiation for a factor
  # optionally pairwise differences,
  # differences from the mean of other levels
  # or the difference from the weighted mean of other levels

  # BUG?: assumes no hyphens in original factor -> change sep

  if( !is.factor (f))   f <- factor( f, levels = f)
  cmat <- contrasts(f)
  nams <- levels(f)
  n <- length(nams)
  Pmats <- lapply(  type, function(x) Pmat( f, x, all = all, sep = sep))
  Pmat <- do.call( rbind, Pmats)
  #disp(Pmat)
  #disp(cmat)
  Cmat <- Pmat %*% cmat
  if( length(unique(rownames(Cmat))) == 1)  single <- TRUE else single <- FALSE
  if (single) Cmat <- rbind( Cmat,'___' = 0 )     #   CHANGE THIS
  fr <- factor(rownames(Cmat), levels=unique(rownames(Cmat)))
  contrasts(fr,n-1) <- Cmat
  names(fr) <- fr
  if (single) fr <- fr[-length(fr)]
  fr
}


#' @describeIn xlevels equivalent of differentiation for a factor
#' @export
dlevels <- function(f, type = "pairwise", all = FALSE, sep = '') {
  # produces the equivalent of differentiation for a factor
  # optionally pairwise differences,
  # differences from the mean of other levels
  # or the difference from the weighted mean of other levels

  # BUG?: assumes no hyphens in original factor -> change sep

  if( !is.factor (f))   f <- factor( f, levels = f)
  cmat <- contrasts(f)
  nams <- levels(f)
  n <- length(nams)
  Pmats <- lapply(  type, function(x) Pmat( f, x, all = all, sep = sep))
  Pmat <- do.call( rbind, Pmats)
  Cmat <- Pmat %*% cmat
  if( length(unique(rownames(Cmat))) == 1)  single <- TRUE else single <- FALSE
  if (single) Cmat <- rbind( Cmat,'<NULL>' = 0 )     #   CHANGE THIS
  fr <- factor(rownames(Cmat), levels=unique(rownames(Cmat)))
  contrasts(fr,n-1) <- Cmat
  names(fr) <- fr
  if (single) fr <- fr[-length(fr)]
  fr
}


# wtd mean of all but 1
# The common problem in all of the following is the generation of a Pmat matrix
# that combines rows of 'cmat' to achieve what it wants

#' Comparison hypothesis matrix
#'
#' @param f a factor or otherwise a vector interpreted as the levels of a
#' factor
#' @param type one "factor","raw","mean","(mean)","<mean>","II",
#'     "cen","cent","center","centre","<centre>" , "<center>" ,
#'     "(center)","(centre)","III","pairwise",  "<others.m>", "<others.c>",
#'     "diff", "diffmean","diffc","diffcen","diffcentre","diffcentre"
#' @param all for some values of \code{type}, indicates whether to use all
#'      contrasts. e.g. \code{type = "pairwise"} will produce pairwise comparisons
#'      in both directions if \code{all == TRUE}
#' @param sep can add spaces in constructed factor levels
#' @export
Pmat <- function( f ,
                  type = c("factor","raw","mean","(mean)","<mean>","II",
                           "cen","cent","center","centre","<centre>" , "<center>" ,
                           "(center)","(centre)","III","pairwise",  "<others.m>", "<others.c>",
                           "diff", "diffmean","diffc","diffcen","diffcentre","diffcentre"),
                  all = FALSE, sep = '')
{
  if( length(type) > 1) return( do.call( rbind, lapply( type, function( t ) Pmat( f, t))))
  type <- match.arg(type)
  if( ! is.factor(f)) f <- factor( f, levels = unique(f))
  switch( type,
          raw =, factor = {
            ret <- diag(length(levels(f)))
            dimnames(ret) <- list(levels(f),levels(f))
            ret
          },
          cen =, cent =, center =, centre = , "III" = , "(centre)" = , "(center)", "<centre>" = , "<center>"  = {
            nlevs = length( levels(f) )
            matrix(rep(1/nlevs, nlevs), nrow = 1, dimnames = list(type, levels(f)))
          },
          mean =, "(mean)"=,  "<mean>"=, "II" = {
            matrix(table(f) / length(f), nrow = 1, dimnames = list(type, levels(f)))
          },
          pairwise = {
            nams <- levels(f)
            n <- length(nams)
            #if( length(nams) == 2) all = TRUE
            if ( all ) { # all comparisons
              plus <- rep( 1:n, n)
              minus <- rep (1:n, each = n)
              drop <- plus == minus
              plus <- plus[!drop]
              minus <- minus[!drop]
            } else {  # no duplicates
              zm <- matrix(1:n, nrow = n, ncol = n)
              plus <- zm[col(zm) < row(zm)]
              minus <- rep(1:(n - 1), (n - 1):1)
            }
            nrows <- length(plus)
            ret <- matrix( 0, nrows, n)
            ret[ cbind(1:nrows,minus) ] <- -1
            ret[ cbind(1:nrows,plus) ] <- 1
            rownames(ret) <- paste( nams[plus],'-', nams[minus], sep='')
            colnames(ret) <- nams
            ret
          } ,
          diff =, diffmean =, "<others.m>"= {
            ret <- Pmat_diffmat( table(f)/length(f) )
            rownames(ret) <- paste(rownames(ret),"-<others.m>", sep = sep)
            ret
          } ,
          diffcentre =, diffcenter =, diffc =, diffcen = , "<others.c>"= {
            levs <- levels(f)
            w <- rep(1/length(levs), length(levs))
            names(w) <- levs
            ret <- Pmat_diffmat( w )
            rownames(ret) <- paste(rownames(ret),"-<others.c>", sep = sep)
            ret
          })
}
#' @describeIn Pmat a utility function to generate a matrix to apply to contrasts to yield differences from
#'     means of other levels.
#' @export
Pmat_diffmat   <- function( w )  {
  ret <- - outer( 1/(1-w) , w, "*")
  diag(ret) <- 1
  dimnames( ret ) <- list( names(w), names(w))
  ret
}
