#' Draw a brace
#'
#' @param x1, y1, x1, y2 endpoints of brace
#' @param right (default TRUE) direction in which brace points
#' @param rad (default 0.2) radius of curl of brace
#' @param tip the tip of the brace (default FALSE)
#' @author Georges Monette, Michael Friendly
#' @return the coordinates of the point of the brace, invisibly
#' @examples
#' plot(c(0,10), c(0,10), type = 'n')
#' lines(brace(0, 0, 0, 5), lwd = 2)
#' text(brace(0, 0, 0, 5, tip = T),"size", adj = c(-.1,.3))
#' lines(brace(2, 0, 8, 5, right = F), lwd = 2)
#' text(brace_(2, 0, 8, 5, right = F), 'tip')
#' text(brace_(2, 0, 8, 5, right = F), 'tip', adj = c(.5,-.3), col = 'red')
#' plot(c(0,10),c(0,1))
#' brace(5,0,5,1) %>% lines
#' lapply(1:10, function(x) {
#'    lines(brace(x,0,x,.1*x,F))
#'    lines(brace(x,0,x-.2*x,.1*x,F),col = 'red')
#' })
#' @export
brace <- function (x1 = 0, y1 = 0, x2 = 0, y2 = 1, right = TRUE, rad = 0.2, tip = FALSE)
{
#   uin only in Splus
#    uin <- par("uin")
#    ux <- uin[1]
#    uy <- uin[2]
        pin <- par("pin")
        usr <- matrix(par("usr"), ncol=2)
        ux <- pin[1] / diff(usr[,1])
        uy <- pin[2] / diff(usr[,2])

        dx <- x2 - x1
        dy <- y2 - y1
#     atan(x,y) from Splus -> atan(a/y) in R
#   alpha <- atan(ux * dx, uy * dy)
        alpha <- atan((ux * dx) / (uy * dy))
        scale <- sqrt((ux * dx)^2 + (uy * dy)^2)
        if (scale > 5 * rad)
                rad <- rad/scale
        qcirc <- cbind(cos((0:10) * pi/20), sin((0:10) * pi/20))
        qcircr <- cbind(cos((10:0) * pi/20), sin((10:0) * pi/20))
        rot <- function(theta) t(cbind(c(cos(theta), sin(theta)),
                                                c(-sin(theta), cos(theta))))
        seg1 <- t(t(rad * qcirc %*% rot(-pi/2)) + c(0, rad))
        seg4 <- t(t(rad * qcirc) + c(0, 1 - rad))
        seg3 <- t(t((rad * qcircr) %*% rot(pi)) + c(2 * rad, 0.5 +
                                                        rad))
        seg2 <- t(t((rad * qcircr) %*% rot(pi/2)) + c(2 * rad, 0.5 -
                                                        rad))
        bra <- rbind(seg1, seg2, seg3, seg4)
        if (!right)
                bra <- bra %*% diag(c(-1, 1))
        bra <- scale * bra %*% rot(-alpha)
        bra <- bra %*% diag(c(1/ux, 1/uy))
        bra <- t(t(bra) + c(x1, y1))
        if(tip) bra[nrow(bra)/2,,drop = FALSE] else bra
}
#' @describeIn brace default tip = TRUE
#' @export
brace_ <- function(...) brace(..., tip = TRUE)
