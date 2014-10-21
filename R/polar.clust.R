#' Convert a hclust object to a polar dendrogram
#' 
#' Any \code{\link{hclust}} object can be converted to a polar dendrogram
#' using this function and then plotted with the standard graphics functions.
#' 
#' @param x \code{\link{hclust}} object.
#' @param labels Chraracter vector of leaf labels.
#' @param col Vector of leaf colours.
#' @param padding Spacing between branches.
#' @param hang Branch lengths of the leaves.
#' @return An object of class \code{polar.clust}
#' @example examples/polar.clust.R
#' @seealso plot.polar.clust
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
#' @import data.table
polar.clust <- function(x, labels=sprintf("#%i", 1:length(x$order)), col, padding=.01, hang=.1*max(x$height)){
    if(padding < 0 || padding > .5){
        warning("Padding must be on the interval [0, .5].")
        padding <- pmin(.5, pmax(0, padding))
    }

    ct <- as.data.table(x$merge)
    setnames(ct, 1:2, c("i", "j"))

    # Branch weights (number of leaves)
    for(k in 1:nrow(ct)){
        ct$ni[k] <- if(ct[k,i] < 0) 1 else sum(ct[i[k], list(ni, nj)])
        ct$nj[k] <- if(ct[k,j] < 0) 1 else sum(ct[j[k], list(ni, nj)])
    }
    
    # Switch branches to keep to smallest one closes to the center
    switch.fun <- function(k, is.left){
        if(ct[k, (ni < nj)] == is.left){
            tmp <- ct$i[k]
            ct$i[k] <<- ct$j[k]
            ct$j[k] <<- tmp
            tmp <- ct$ni[k]
            ct$ni[k] <<- ct$nj[k]
            ct$nj[k] <<- tmp
        }
        if(ct[k,i] > 0) switch.fun(ct[k,i], TRUE)
        if(ct[k,j] > 0) switch.fun(ct[k,j], FALSE)
    }
    switch.fun(nrow(ct), ct[nrow(ct), ni > nj])

    # Branch sectors
    ct$m <- ct$ni/(ct$ni + ct$nj)
    ct$s0[nrow(ct)] <- 0
    ct$s1[nrow(ct)] <- 1
    for(k in nrow(ct):1){
        s <- ct[k, quantile(c(s0, s1), c(m*padding, m-m*padding, m+(1-m)*padding, 1-(1-m)*padding))]
        if(ct[k,i] > 0){
            ct[i[k]]$s0 <- s[1]
            ct[i[k]]$s1 <- s[2]
        }
        if(ct[k,j] > 0){
            ct[j[k]]$s0 <- s[3]
            ct[j[k]]$s1 <- s[4]
        }
    }

    # Branch radius
    I <- ct[,ifelse(i > 0, i, NA)]
    J <- ct[,ifelse(j > 0, j, NA)]
    ct$r <- tail(x$height,1) - x$height
    ct$ri <- ifelse(is.na(I), ct$r + hang, ct[I, r])
    ct$rj <- ifelse(is.na(J), ct$r + hang, ct[J, r])

    # Branch angles
    ct$a <- ct[,s0+m*(s1-s0)]
    ct$ai <- ifelse(is.na(I), ct[,s0+m*(s1-s0)/2], ct[I,s0+m*(s1-s0)])
    ct$aj <- ifelse(is.na(J), ct[,s1-(1-m)*(s1-s0)/2], ct[J,s0+m*(s1-s0)])

    class(ct) <- c("polar.clust", class(ct))
    ct
}

order.polar.clust <- function(x){
    rbind(x[i < 0, list(i=-i, a=ai)],
          x[j < 0, list(i=-j, a=aj)])[order(a), i]
}

#' colors-set or labels to a polar dendrogram
#'
#' @param object The \code{\link{polar.clust}} object.
#' @param value The colors or labels of the cluster leaves. If \code{value}
#'   has the same length as the number of leaves it corresponds to the leaves
#'   (ordered as in the original dataset). Otherwise \code{value} corresponds
#'   to branches (indexed by \code{\link{cutree.polar.clust}}).
#' @return The modified cluster object
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @rdname colors-set
#' @aliases colors<-
#' @export
`colors<-` <- function(object, value){
    UseMethod("colors<-")
}
#' @rdname colors-set
#' @export
`colors<-.polar.clust` <- function(object, value){
    if(is.null(value)){
        object$coli <- object$colj <- NULL
        return(object)
    }
    I <- ifelse(object$i > 0, object$i, NA)
    J <- ifelse(object$j > 0, object$j, NA)

    n <- nrow(object) + 1
    if(length(value) < n)
        value <- value[cutree.polar.clust(object, length(value))]

    if(is.factor(value)){
        require(RColorBrewer)
        value <- brewer.pal(9, "Set1")[as.integer(value)]
    }
    value <- col2rgb(value)/255

    rgbi <- value[,ifelse(object$i < 0, -object$i, NA)]
    rgbj <- value[,ifelse(object$j < 0, -object$j, NA)]
    object$colj <- object$coli <- as.character(NA)
    for(k in 1:nrow(object)){
        if(any(is.na(rgbi[,k]))){
            ik <- I[k]
            rgbi[,k] <- (object$ni[ik]*rgbi[,ik] + object$nj[ik]*rgbj[,ik]) / object[ik,sum(ni,nj)]
        }
        object$coli[k] <- do.call(rgb, as.list(rgbi[,k]))
        if(any(is.na(rgbj[,k]))){
            jk <- J[k]
            rgbj[,k] <- (object$ni[jk]*rgbi[,jk] + object$nj[jk]*rgbj[,jk]) / object[jk,sum(ni,nj)]
        }
        object$colj[k] <- do.call(rgb, as.list(rgbj[,k]))
    }
    object
}
#' @rdname colors-set
#' @export
`labels<-` <- function(object, value){
    UseMethod("labels<-")
}
#' @rdname colors-set
#' @export
`labels<-.polar.clust` <- function(object, value){
    if(is.null(value)){
        if("labi" %in% names(object)) object$labi <- object$labj <- NULL
        return(object)
    }

    if(is.null(object$labi)) object$labi <- NA
    if(is.null(object$labj)) object$labj <- NA
    if(length(value) == nrow(object)+1){
        object$labi <- ifelse(object$i < 0, value[abs(object$i)], object$labi)
        object$labj <- ifelse(object$j < 0, value[abs(object$j)], object$labj)
    } else {
        ii <- get.branch.index(x, k=length(value))
        ii <- ii[order(object$r[ii])]
        newi <- value[match(object$i, ii)]
        newj <- value[match(object$j, ii)]
        object$labi <- ifelse(is.na(newi), object$labi, newi)
        object$labj <- ifelse(is.na(newj), object$labj, newj)
    }
    object
}


#' Plot functions for polar dendrograms
#' 
#' Class specific methods for class \code{\link{polar.clust}}. If you are not
#' happy with the looks of the plot you can build it up manually by calling it
#' with \code{type="n"} and add the annotation functions.
#' 
#' @param x \code{\link{polar.clust}} object.
#' @param type Plot type, see \code{\link{plot}}.
#'   you add whatever layers you want manually.
#' @param lines Whether to plot the trees branches.
#' @param points Whether to plot the trees leaves.
#' @param fill Whether to fill the area under labelled branches, see \code{\link{branches}}.
#' @param labels Whether to plot the leaves' labels.
#' @param label.col Whether to colour labels by leaf colour.
#' @param label.size Relative text size of labels.
#' @param expand Controls how much white space should surround the plot.
#'   Use \code{expand=1} for no white space and larger numbers for more.
#' @param axes Whether to plot axes or not.
#' @param ... Sent to \code{\link{plot}}.
#' @param bg.col Background colour.
#' @param grid.col Radial grid color.
#' @return Nothing, produces a plot.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
plot.polar.clust <- function(x, type=c("b", "l", "p", "n"), sector, fill=TRUE, labels, label.col, label.size=1, las=2, expand, axes=TRUE, ..., bg.col, grid.col="gray90"){
    type <- match.arg(type)
    if(missing(labels)) labels <- "labi" %in% names(x)
    if(missing(expand)) expand <- if(labels) 1.5 else 1.04
    if(length(expand) == 1) expand <- rep(expand, 2)
    xl <- c(-1,1) * max(x$r)*expand[1]
    yl <- c(-1,1) * max(x$r)*expand[2]

    plot(0, 0, type="n", xlab="", ylab="", xlim=xl, ylim=yl, axes=FALSE, ann=FALSE, ...)
    if(!missing(bg.col)){
        do.call(rect, as.list(c(par("usr")[c(1,3,2,4)], border=NA, col=bg.col)))
    }
    if(missing(sector))
        sector <- pi/2+c(.15*axes, 2*pi)

    # Grid circles
    r <- pretty(c(0, par("usr")[2]))
    r <- r[2:(length(r)-1)]
    n.seg <- round(400*r/max(r))
    for(i in seq_along(r)){
        a <- 2*pi*0:n.seg[i]/n.seg[i]
        segments(r[i]*sin(a[-n.seg[i]-1]), r[i]*cos(a[-n.seg[i]-1]),
                 r[i]*sin(a[-1]), r[i]*cos(a[-1]), col=grid.col)
    }

    # Branches
    if(fill)
        branches(x, sector=sector, fill=TRUE, lightness=.7)
    if(missing(type) || type %in% c("b", "l"))
        lines(x, lwd=c(1,10), sector=sector)
    if(missing(type) || type %in% c("b", "p"))
        points(x, sector=sector, pch=20)
    if(labels)
        text(x, sector=sector, cex=label.size, col=label.col, las=las)
    if(axes) axis(1, at=r, pos=0)
}

#' @param lwd Minimum and maximum branch widths.
#' @param col Colour.
#' @param sector The sector of the plot area to use, e.g. set
#'   \code{sector=c(1/2, 2)*pi} to leave one quadrant empty, giving room to put
#'   a \code{\link{legend}}.
#' @rdname plot.polar.clust
#' @export
lines.polar.clust <- function(x, lwd=c(1,10), col, sector=c(0, 2*pi), ...){
    if(missing(col)){
        col <- if("coli" %in% names(x)) c(x$coli, x$colj) else par("fg")
    } else if(is.null(col) || col %in% c(FALSE, NA)){
        col <- par("fg")
    }
    w <- c(x$ni, x$nj)
    my.sin <- function(x) sin(x*diff(range(sector))+sector[1])
    my.cos <- function(x) cos(x*diff(range(sector))+sector[1])
    segments(rep(x[, r*my.sin(a)], 2),
             rep(x[, r*my.cos(a)], 2),
             c(x[, ri*my.sin(ai)], x[, rj*my.sin(aj)]),
             c(x[, ri*my.cos(ai)], x[, rj*my.cos(aj)]),
             lwd = (w-min(w))/diff(range(w))*diff(lwd) + lwd[1],
             col = col,
             ...)
}
#' @param r Radius at which to plot the labels.
#' @param cex Relative size, see \code{\link{par}} for details.
#' @param las Label orientation, see \code{\link{par}}.
#' @rdname plot.polar.clust
#' @export
text.polar.clust <- function(x, r, sector=c(0,2*pi), cex=par("cex"), col, las=par("las"), ...){
    if(missing(col)){
        col <- if("coli" %in% names(x)){
            c(x[!is.na(labi), coli], x[!is.na(labj), colj])
        } else {
            par("fg")
        }
    } else if(is.null(col) || col %in% c(FALSE, NA)){
        col <- par("fg")
    }
    r.max <- if(missing(r)) max(x$ri, x$rj) else r
    ang.i <- x[!is.na(labi), (s0+(s1-s0)*m/2)*diff(range(sector))+sector[1]]
    ang.j <- x[!is.na(labj), (s1-(s1-s0)*(1-m)/2)*diff(range(sector))+sector[1]]
    lab <- rbind(
        x[!is.na(labi), list(
            labels = labi,
            x = 1.02*r.max*sin(ang.i),
            y = 1.02*r.max*cos(ang.i),
            a = ang.i*180/pi
        )],
        x[!is.na(labj), list(
            labels = labj,
            x = 1.02*r.max*sin(ang.j),
            y = 1.02*r.max*cos(ang.j),
            a = ang.j*180/pi
        )]
    )
    lab$srt <- switch(las+1, 
        -lab$a + 180*((lab$a+90) %% 360 > 180), # parallel to axis
        0,                                      # horizontal
        -lab$a + 90 + 180*(lab$a %% 360 > 180), # perpendicular to axis
        90                                      # vertical
    )
    lab$adjx <- switch(las+1,
        0.5,
        as.numeric(lab$a %% 360 > 180),
        as.numeric(lab$a %% 360 > 180),
        as.numeric((lab$a + 90) %% 360 > 180)
    )
    lab$adjy <- switch(las+1,
        as.numeric((lab$a + 90) %% 360 > 180),
        as.numeric((lab$a + 90) %% 360 > 180),
        0.5,
        as.numeric(lab$a + 90 %% 360 > 180)

    )
    if(nrow(lab) > 0) for(i in 1:nrow(lab))
        with(lab[i], text(x, y, labels, srt=srt, adj=c(adjx, adjy), cex=cex, col=col[i], ...))
}
#' @rdname plot.polar.clust
#' @export
points.polar.clust <- function(x, sector=c(0,2*pi), col, ...){
    if(missing(col)){
        col <- if("coli" %in% names(x)) c(x[i < 0, coli], x[j < 0, colj]) else par("fg")
    } else if(is.null(col) || col %in% c(FALSE, NA)){
        col <- par("fg")
    }
    my.sin <- function(x) sin(x*diff(range(sector))+sector[1])
    my.cos <- function(x) cos(x*diff(range(sector))+sector[1])
    pnt <- rbind(
        x[i < 0, list(
            x = ri*my.sin(ai),
            y = ri*my.cos(ai)
        )], 
        x[j < 0, list(
            x = rj*my.sin(aj),
            y = rj*my.cos(aj)
        )]
    )
    with(pnt, points(x=x, y=y, col=col, ...))
}

