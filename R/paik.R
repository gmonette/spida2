#' Paik-Agresti diagrams
#' New version of Paik-Agresti diagrams
#' 
#' @param form formula of the form y ~ x + z + count
#' @param data data frame
#' @param ylab optional parameter
#' @param xlab optional parameter
#' @param cex.factor default 10: scales the size of circles
#' @param ... additional parameters for \code{lattice::xyplot} such
#'        \code{xlim}, \code{ylim}, \code{main}, etc.
#' @examples
#' paik(verdict ~ d.race + v.race | count, death.penalty)
#' paik(verdict ~ d.race + v.race | count, death.penalty, 
#'     ylab = 'verdict', main = 'Paik-Agresti Diagram')
#' paik(Status ~ Gender + Dept | count, Berkeley)
#' gd(lty = 1:4)
#' paik(Status ~ Gender + Dept | count, Berkeley,
#'     zlab = 'Department')
#' @export
paik <- function(form, data, ylab, xlab = xn, zlab = zn,
                 cex.factor = 10, cex.zlab = 1,...){
  vars <- all.vars(form)
  y <- data[[yn <- vars[1]]]
  ynum <- TRUE
  ylab_miss <- missing(ylab)
  if(ylab_miss) ylab <- paste('average of',yn)
  if(!is.numeric(y)) {
    ynum <- FALSE
    y <- as.factor(y)
    if(length(unique(y)) != 2) warning('if y is a factor it should have 2 unique values')
    y2 <- unique(y)[2]
    y <- as.numeric(y) - 1
    if(ylab_miss) ylab <- paste('proportion of', yn,'equal to',y2)
  }
  x <- as.factor(data[[xn <- vars[2]]])
  z <- as.factor(data[[zn <- vars[3]]])
  wt <- if(length(vars) == 4) {
    data[[vars[4]]]
  } else {
    y*0 + 1
  }
  dd <- data.frame(y=y, x = x, z = z)
  dd <- dd[rep(1:nrow(dd), wt),]
  
  fit <- lm(y ~ x * z,  dd)
  fit_marg <- lm(y ~ x, dd)
  dcond <- up(dd, ~ x + z)
  dmarg <- up(dd, ~ x)
  dcond$y <- predict(fit, newdata = dcond)
  dmarg$y <- predict(fit_marg, newdata = dmarg)
  # print(dcond)
  # print(dmarg)
  library(latticeExtra)
  dcond$z <- reorder(dcond$z, -dcond$y)
  print(xyplot(y ~ x, dcond, groups = z, type = 'l', lwd = 2,
         xlab = xlab, ylab = ylab, ...,
         auto.key = list(title = zlab,cex.title=cex.zlab))+
    layer_(panel.grid(v=-1,h=-1))+
    xyplot(y ~ x, dmarg, type = 'b', lwd =3, pch = 16, col = 'black')+
    xyplot(y ~ x, dcond, type = 'p',  pch = 1, cex = cex.factor*dcond$Freq/max(dcond$Freq))
  )
  invisible(list(dcond, dmarg))
}

