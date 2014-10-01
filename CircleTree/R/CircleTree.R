##' Impute missing values with a given fill value
##'
##' @param x Object.
##' @param fill The value that should replace missing values in \code{x}.
##' @return A version of \code{x} with no missing values.
##' @author Christofer \enc{Bäcklin}{Backlin}
na.fill <- function(x, fill){
    x[is.na(x)] <- fill
    x
}


##' Convert a hclust object to a circular dendrogram
##' 
##' Any \code{\link{hclust}} object can be converted to a circular dendrogram
##' using this function and then plotted with the standard graphics functions.
##' 
##' @param x \code{\link{hclust}} object.
##' @param labels Chraracter vector of leaf labels.
##' @param col Vector of leaf colours.
##' @param padding Spacing between branches.
##' @param hang Branch lengths of the leaves.
##' @return An object of class \code{CircleTree}
##' @examples
##' my.data <- iris[-5] + .5*matrix(rnorm(nrow(iris)*(ncol(iris)-1)), nrow(iris))
##' cl <- hclust(dist(my.data))
##' 
##' # Color according to class
##' x <- CircleTree(cl,
##'                 labels=sprintf("%s %i", iris$Species, 1:150),
##'                 col=iris$Species,
##'                 padding=.01)
##' plot(x, label.col=TRUE)
##' 
##' # Color branches
##' x <- CircleTree(cl, labels=sprintf("%s %i", iris$Species, 1:150),
##'                 col=cutree(cl, 3), padding=0, hang=.5)
##' plot(x, labels=FALSE)
##' @seealso plot.CircleTree
##' @author Christofer \enc{Bäcklin}{Backlin}
##' @export
##' @import data.table
CircleTree <- function(x, labels=sprintf("#%i", 1:length(x$order)), col, padding=.01, hang=.1*max(x$height)){
    if(padding < 0 || padding > .5){
        warning("Padding must be on the interval [0, .5].")
        padding <- pmin(.5, pmax(0, padding))
    }

    ct <- data.table(i = ifelse(x$merge[,1] > 0, x$merge[,1], NA),
                     j = ifelse(x$merge[,2] > 0, x$merge[,2], NA))
    if(!missing(labels)){
        ct$labi <- ifelse(x$merge[,1] < 0, labels[abs(x$merge[,1])], NA)
        ct$labj <- ifelse(x$merge[,2] < 0, labels[abs(x$merge[,2])], NA)
    }

    # Branch weights (number of leaves)
    for(k in 1:nrow(ct)){
        ct$ni[k] <- if(is.na(ct$i[k])) 1 else sum(ct[i[k], list(ni, nj)])
        ct$nj[k] <- if(is.na(ct$j[k])) 1 else sum(ct[j[k], list(ni, nj)])
    }
    
    # Switch branches to keep to smallest one closes to the center
    switch.fun <- function(k, is.left){
        if(ct[k, (ni < nj)] == is.left){
            x$merge[k,] <<- rev(x$merge[k,])
            tmp <- ct$i[k]
            ct$i[k] <<- ct$j[k]
            ct$j[k] <<- tmp
            tmp <- ct$labi[k]
            ct$labi[k] <<- ct$labj[k]
            ct$labj[k] <<- tmp
            tmp <- ct$ni[k]
            ct$ni[k] <<- ct$nj[k]
            ct$nj[k] <<- tmp
        }
        if(!is.na(ct$i[k])) switch.fun(ct$i[k], TRUE)
        if(!is.na(ct$j[k])) switch.fun(ct$j[k], FALSE)
    }
    switch.fun(nrow(ct), ct[nrow(ct), ni > nj])

    # Branch lengths
    ct$li <- na.fill(x$height - x$height[ct$i], hang)
    ct$lj <- na.fill(x$height - x$height[ct$j], hang)

    # Colors
    if(missing(col)) col <- "black"
    if(is.factor(col)){
        require(RColorBrewer)
        col <- brewer.pal(9, "Set1")[as.integer(col)]
    }
    col <- col2rgb(col)/255
    n <- length(x$order)
    if(ncol(col) != n)
        col <- col[,rep(1:ncol(col), ceiling(n/ncol(col)))[1:n]]
    rgbi <- col[,ifelse(x$merge[,1] < 0, abs(x$merge[,1]), NA)]
    rgbj <- col[,ifelse(x$merge[,2] < 0, abs(x$merge[,2]), NA)]
    ct$colj <- ct$coli <- as.character(NA)
    for(k in 1:nrow(ct)){
        if(any(is.na(rgbi[,k]))){
            ik <- ct$i[k]
            rgbi[,k] <- (ct$ni[ik]*rgbi[,ik] + ct$nj[ik]*rgbj[,ik]) / ct[ik,sum(ni,nj)]
        }
        ct$coli[k] <- do.call(rgb, as.list(rgbi[,k]))
        if(any(is.na(rgbj[,k]))){
            jk <- ct$j[k]
            rgbj[,k] <- (ct$ni[jk]*rgbi[,jk] + ct$nj[jk]*rgbj[,jk]) / ct[jk,sum(ni,nj)]
        }
        ct$colj[k] <- do.call(rgb, as.list(rgbj[,k]))
    }

    # Branch sectors
    ct$m <- ct$ni/(ct$ni + ct$nj)
    ct$s0[nrow(ct)] <- 0
    ct$s1[nrow(ct)] <- 2*pi
    for(k in nrow(ct):1){
        s <- ct[k, quantile(c(s0, s1), c(m*padding, m-m*padding, m+(1-m)*padding, 1-(1-m)*padding))]
        if(!is.na(ct$i[k])){
            ct[i[k]]$s0 <- s[1]
            ct[i[k]]$s1 <- s[2]
        }
        if(!is.na(ct$j[k])){
            ct[j[k]]$s0 <- s[3]
            ct[j[k]]$s1 <- s[4]
        }
    }

    # Branch radius
    ct$r <- tail(x$height,1) - x$height
    ct$ri <- ifelse(is.na(ct$i), ct$r + hang, ct[ct$i, r])
    ct$rj <- ifelse(is.na(ct$j), ct$r + hang, ct[ct$j, r])

    # Branch angles
    ct$a <- ct[,s0+m*(s1-s0)]
    ct$ai <- ifelse(is.na(ct$i), ct[,s0+m*(s1-s0)/2], ct[ct$i,s0+m*(s1-s0)])
    ct$aj <- ifelse(is.na(ct$j), ct[,s1-(1-m)*(s1-s0)/2], ct[ct$j,s0+m*(s1-s0)])

    class(ct) <- c("CircleTree", class(ct))
    ct
}


##' Plot functions for circular dendrograms
##' 
##' Class specific methods for class \code{\link{CircleTree}}.
##' 
##' @param x \code{\link{CircleTree}} object.
##' @param lines Whether to plot the trees branches.
##' @param points Whether to plot the trees leaves.
##' @param labels Whether to plot the leaves' labels.
##' @param label.size Relative text size of labels.
##' @param label.col Whether to colour labels by leaf colour.
##' @param expand Controls how much white space should surround the plot.
##'   Use \code{expand=1} for no white space and larger numbers for more.
##' @param lwd Minimum and maximum branch widths.
##' @param ... Sent to \code{\link{plot}}.
##' @param bg Background colour.
##' @return Nothing, produces a plot.
##' @seealso CircleTree
##' @author Christofer \enc{Bäcklin}{Backlin}
##' @export
plot.CircleTree <- function(x, lines=TRUE, points=TRUE, labels, label.size=1, label.col=FALSE, expand, lwd=c(1,10), ..., bg){
    if(missing(labels)) labels <- "labi" %in% names(x)
    if(missing(expand)) expand <- if(labels) 2 else 1.04

    lim <- c(-1,1) * max(x$r)*expand
    plot(0, 0, type="n", xlab="", ylab="", xlim=lim, ylim=lim, ...)
    if(!missing(bg)){
        do.call(rect, as.list(c(par("usr")[c(1,3,2,4)], border=NA, col=bg)))
        r <- setdiff(unique(abs(pretty(par("usr")[1:2]))), 0)
        symbols(rep(0, length(r)), circles=r, inches=FALSE, fg="white", add=TRUE)
    }
    if(lines)
        lines(x, lwd)
    if(points)
        points(x, pch=20)
    if(labels)
        labels(x, cex=label.size, col=label.col)
}

##' @param col Colour.
##' @rdname plot.CircleTree
##' @export
lines.CircleTree <- function(x, lwd=c(1,10), col, ...){
    w <- c(x$ni, x$nj)
    segments(rep(x[, r*sin(a)], 2),
             rep(x[, r*cos(a)], 2),
             c(x[, ri*sin(ai)], x[, rj*sin(aj)]),
             c(x[, ri*cos(ai)], x[, rj*cos(aj)]),
             lwd = (w-min(w))/diff(range(w))*diff(lwd) + lwd[1],
             col=if(missing(col)) c(x$coli, x$colj) else col,
             ...)
}
##' @param object Same as \code{x}. Named differently to keep S3 consistency.
##' @param cex Relative size, see \code{\link{par}} for details.
##' @rdname plot.CircleTree
##' @export
labels.CircleTree <- function(object, cex=par("cex"), col=FALSE, ...){
    lab <- rbind(
        object[!is.na(labi), list(
            labels = labi,
            x = 1.02*max(ri)*sin(ai),
            y = 1.02*max(ri)*cos(ai),
            adj = ai > pi,
            srt = 90-(ai %% pi)*180/pi,
            col = if(col) coli else par("fg")
        )],
        object[!is.na(labj), list(
            labels = labj,
            x = 1.02*max(rj)*sin(aj),
            y = 1.02*max(rj)*cos(aj),
            adj = aj > pi,
            srt = 90-(aj %% pi)*180/pi,
            col = if(col) colj else par("fg")
        )]
    )
    for(i in 1:nrow(lab))
        with(lab[i], text(x, y, labels, srt=srt, adj=c(adj, .5), cex=cex, col=col, ...))
}
##' @rdname plot.CircleTree
##' @export
points.CircleTree <- function(x, col, ...){
    pnt <- rbind(
        x[is.na(i), list(
            x = ri*sin(ai),
            y = ri*cos(ai),
            col = coli
        )], 
        x[is.na(j), list(
            x = rj*sin(aj),
            y = rj*cos(aj),
            col = colj
        )]
    )
    points(x=pnt$x, y=pnt$y, col=if(missing(col)) pnt$col else col, ...)
}


#print.CircleTree <- function(x, ...){
#    cat("Circle tree representation of hclust object. Print function not yet implemented.\n")
#}



