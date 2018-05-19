# @param zoom (=FALSE). Should first category, showing proportions of complete data, have its bar on the right scaled down. (NOT YET IMPLEMENTED)


#' Tabulate and visualize missing values
#'
#' Tabulates missing values in variables and patterns of missing values among
#' variables. It produces a missing value plot displaying different patterns of
#' missing values among variables and missing value counts.
#'
#' This function tabulates the missing values in variables and the patterns of
#' missings. The input data should be of class matrix or dataframe. The
#' function also displays a missing value plot. Each row represents a pattern
#' of missings, where the variables with missing values are shown in yellow.
#' The barchart below shows the proportions of missing values by variable. The
#' barchart on the right displays the pattern frequencies. This barchart uses
#' the concept of censored zooming. The topmost bar, representing no missing
#' values in any variable, is likely to be much larger than the others and will
#' then not be drawn to scale, but bordered in red to show this. Ordering of
#' variables and patterns of missings based on counts can be very useful for
#' exploring the structure of missings in a dataset.
#'
#' @param x the dataset to be examined. It should be of class matrix or
#' dataframe.
#' @param sortby a logical value indicating whether the variables, patterns or
#' both are sorted by the numbers of missing values. Possible values are
#' "variable" (sort by variable name), "row" (sort by pattern counts), "column"
#' (sort by variable counts) and "both" (sort by variable and pattern counts).
#' @return
#' \item{finaltable }{A table of class dataframe containing number of missing
#' values in variable and missing pattern}
#' @author Waqas Ahmed Malik, with small modifications.
#' @keywords missing values plot
#' @examples
#' \dontrun{
#' data(hs)
#' zz <- hs
#' missings <- cbind(sample(1:50, 100, replace = T), sample(2:6, 100, replace = T))
#' missings <- unique(missings)
#' nrow(missings)
#' dim(missings)
#' zz[missings] <- rep(NA, nrow(missings))
#' tablemissing(zz, sortby = "variable")
#' }
#' @export
tablemissing <-
  function (x, sortby="variable", bot = 6, mar = c(bot,1,1.5,1),
            adj=c(0,-.1), xpd=TRUE, srt = -60, cex = .8, zoom = FALSE, ...)
  {

    x1 <- as.numeric(apply(x, 2, function(x) length(which(is.na(x)))))
    x1 <- c(x1, nrow(x))

    z1 <- ifelse(is.na(x), 0, 1)
    tab = table(apply(z1, 1, paste, collapse = ","))
    tab = tab[order(names(tab), decreasing=TRUE)]
    tab = data.frame(combination = names(tab), count = as.numeric(tab))
    tabp <- t(apply(tab, 1, function(x) { as.numeric(unlist(strsplit(x, ",", fixed = TRUE))) } ))
    tabp <- as.data.frame(tabp)
    tabp <- rbind(tabp, x1)

    names(tabp) <- c(names(x),"Total")
    row.names(tabp) <- c( seq(1,nrow(tab)), "Total")

    if ( sortby=="variable"){
      tabfinal<-tabp
    }
    if(sortby=="row"){
      tabfinal <- tabp[-nrow(tabp),]
      tabfinal <- tabfinal[order(tabfinal$Total, decreasing = TRUE),]
      tabfinal <- rbind(tabfinal, tabp[nrow(tabp),])
    }

    if (sortby == "column"){
      tabfinal <- tabp[,-ncol(tabp)]
      vals <- unlist(tabfinal[nrow(tabfinal), ])
      tabfinal <- tabfinal[order(vals, decreasing = TRUE)]
      tabfinal <- cbind(tabfinal, Total=tabp$Total)
    }

    if (sortby == "both"){
      tabf <- tabp[-nrow(tabp),]
      tabf <- tabf[order(tabf$Total, decreasing = TRUE),]
      tabf <- rbind(tabf, tabp[nrow(tabp),])

      tabfinal <- tabf[,-ncol(tabf)]
      vals <- unlist(tabfinal[nrow(tabfinal), ])
      tabfinal <- tabfinal[order(vals, decreasing = TRUE)]
      tabfinal <- cbind(tabfinal, Total=tabf$Total)
    }

    finaltable<<- tabfinal
    finaltable

    opar <- par(mar=mar)
    on.exit(par(opar))

    nop=nrow(finaltable)-1
    nov=ncol(finaltable)-1
    width= 100/(nov)
    height=10
    x1=0
    x2=width
    y1=30
    y2=y1+height

    pylim=y1+10*nop

    plot(10,20, type="n", xlim=c(0,120), ylim=c(0,pylim), axes=FALSE, xlab="", ylab="", main="Missing Value Patterns")


    for (i in nop:1){
      for (j in 1:nov){
        if (finaltable[i,j]==0){
          polygon( c(x1,x2,x2,x1), c(y1,y1,y2,y2), col="yellow", border="yellow3")
        }else{
          polygon( c(x1,x2,x2,x1), c(y1,y1,y2,y2), col="blue", border="skyblue")
        }
        x1=x1+width
        x2=x2+width
      }
      x1=0
      x2=width
      y1=y1+height
      y2=y2+height
    }

    bx1=width/4
    bx2=3*bx1
    by1=5
    by3=25
    bsize=20
    for (i in 1:nov){
      m=finaltable[nop+1,i]/finaltable[nop+1,nov+1]*bsize
      p=bsize-m
      by2=by1+p

      polygon( c(bx1,bx2,bx2,bx1), c(by1,by1,by2,by2), col="blue", border=NA)
      polygon( c(bx1,bx2,bx2,bx1), c(by2,by2,by3,by3), col="red", border=NA)
      text(bx1,0, names(finaltable)[i],
           srt = srt, adj= adj, xpd = xpd, cex = cex, ...)
      bx1=bx1+width
      bx2=bx1+width/2
    }


    px1=105
    py1=30
    py2=py1+7
    for ( i in nop:2){
      if ( sum(finaltable[1,1:nov])==nov){
        psize=finaltable[i,nov+1]/(finaltable[nop+1,nov+1]-finaltable[1,nov+1])*20
      }else{
        psize=finaltable[i,nov+1]/(finaltable[nop+1,nov+1])*20
      }
      if (psize<0.2){
        psize=0.2
      }
      px2=px1+psize
      polygon( c(px1,px2,px2,px1), c(py1,py1,py2,py2), col="blue", border=NA)

      py1=py1+10
      py2=py2+10
    }
    psize=finaltable[1,nov+1]/finaltable[nop+1,nov+1]*20
    px2=px1+psize

    if ( sum(finaltable[1,1:nov])==nov){
      polygon( c(px1,px2,px2,px1), c(py1,py1,py2,py2), col="blue", border="red")
    }else{
      polygon( c(px1,px2,px2,px1), c(py1,py1,py2,py2), col="blue", border=NA)
    }

    finaltable
  }


