#' Cut a tree inte groups
#'
#' A wrapper that allows cutting polar dendrograms with \code{\link{cutree}}.
#'
#' @param tree \code{\link{polar.clust}} object.
#' @param k Number of groups.
#' @param h Height at which to cut the tree.
#' @seealso cutree
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
cutree.polar.clust <- function(tree, k, h){
    if(missing(k) && missing(h) && "labi" %in% names(tree)){
        labs <- rep(NA, nrow(tree)+1)
        lab.fun <- function(i, lab){
            my.lab <- if(is.na(tree[i,labi])) lab else tree[i,labi]
            if(tree[i,i < 0]){
                labs[-tree[i,i]] <<- my.lab
            } else {
                lab.fun(tree[i,i], my.lab)
            }
            my.lab <- if(is.na(tree[i,labj])) lab else tree[i,labj]
            if(tree[i,j < 0]){
                labs[-tree[i,j]] <<- my.lab
            } else {
                lab.fun(tree[i,j], my.lab)
            }
        }
        lab.fun(nrow(tree), "Unlabelled")
        factor(labs)
    } else {
        if(missing(k))
            k <- sum((tree$r < h) != (tree$ri < h)) + sum((tree$r < h) != (tree$rj < h))
        cutree(list(merge = as.matrix(tree[,list(i,j)])), k)
    }
}

#' Extracting line numbers of desired branches
#'
#' Complement to \code{\link{cutree.polar.clust}}.
#' 
#' @param x \code{\link{polar.clust}} object.
#' @param k Number of branches to split \code{x} into.
#' @param h Height at which to split \code{x}.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
get.branch.index <- function(x, k, h){
    if(missing(k) == missing(h))
        stop("You must supply `k` or `h`, but not both.")
    if(!missing(k)){
        branches <- rep(nrow(x), k)
        for(b in 2:k){
            split <- which.min(x[branches[1:(b-1)], pmin(ri,rj)])
            branches[c(split,b)] <- unlist(x[branches[split],list(i,j)])
        }
        branches
    } else {
        c(x[r < h & ri >= h, i], x[r < h & rj >= h, j])
    }
}

#' Color branches of a polar dendrogram
#'
#' This function can either draw a line in the margin outside each branch
#' (\code{fill=FALSE}) or color the area underneath a branch (\code{fill=TRUE}).
#' 
#' @param x \code{\link{polar.clust}} object.
#' @param fill Annotation style, see the details of this function.
#' @param lwd Line width, should be pretty large to be clearly visible.
#' @param sector See \code{\link{plot.polar.clust}}.
#' @param col Branch colors. If omitted it will be taken from \code{x} (if
#'   it has first been specified with \code{\link{colors<-}}).
#' @param alpha Opacity values (0 is transparent, 1 is opaque).
#' @param lightness Color lightness (0 is the original color, 1 is white, and
#'   values in between are interpolated).
#' @param ... Sent to \code{\link{polygon}} or \code{\link{lines}}.
#' @example examples/polar.clust.R
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
branches <- function(x, fill=FALSE, lwd, sector=c(0, 2*pi), col, alpha=1, lightness, ...){
    if("labi" %in% names(x)){
        ii <- c(x[!is.na(labi),i], x[!is.na(labj),j])
        ii <- ii[ii > 0]
        ii <- ii[order(x[ii,r])]
    } else {
        ii <- NULL
    }
    if(missing(lwd)){
        if(!fill){
            lwd <- 20
        } else if(alpha < 1){
            lwd <- NULL # The overlapping border and fill of polygon results in
                        # an ugly inner line if alpha < 1
        } else {
            lwd <- 40 - 4.3*log(nrow(x)) # ~20 for 100 observations and ~10 for 1000
        }
    }
    if(length(ii) > 0){
        get.i <- function(i) rbind(x[i,list(r,a)],
            if(x[i,i] > 0) get.i(x[i,i]) else x[i,list(r=ri, a=ai)])
        get.j <- function(i) rbind(
            if(x[i,j] > 0) get.j(x[i,j]) else x[i,list(r=rj, a=aj)],
            x[i,list(r,a)])
        get.leaves <- function(i) rbind(
            if(x[i,i > 0]) get.leaves(x[i,i]) else x[i,list(r=ri, a=ai)],
            if(x[i,j > 0]) get.leaves(x[i,j]) else x[i,list(r=rj, a=aj)])
        if(missing(col)){
            if("coli" %in% names(x)){
                col <- sweep(
                    sweep(col2rgb(x[ii,coli]), 2, x[ii,ni], "*") +
                    sweep(col2rgb(x[ii,colj]), 2, x[ii,nj], "*"),
                    2, x[ii, ni+nj], "/")/255
            } else {
                require(RColorBrewer)
                col <- col2rgb(brewer.pal(length(ii), "Set1"))/255
            }
        } else {
            col <- col2rgb(col)/255
        }
        if(!missing(lightness))
            col <- col + lightness*(1-col)
        if(!missing(alpha)){
            if(fill && alpha < 1 && !is.null(lwd))
                warning("Specifying both `fill = TRUE` and `alpha < 1` with `lwd` not null is violating the laws of good taste (just look at it!).")
            col <- rbind(col, alpha)
        }
        col <- apply(col, 2, function(x) do.call(rgb, as.list(x)))

        for(i in seq_along(ii)){
            if(fill){
                b <- rbind(get.i(ii[i]), get.leaves(ii[i])[,list(r=max(r), a)], get.j(ii[i]))
                polygon(x = b$r*sin(b$a*diff(sector)+sector[1]),
                        y = b$r*cos(b$a*diff(sector)+sector[1]),
                        border=if(is.null(lwd)) NA else col[i], col=col[i], lwd=lwd, ...)
            } else {
                b <- get.leaves(ii[i])[,list(r=max(r), a)]
                lines(x = b$r*sin(b$a*diff(sector)+sector[1]),
                      y = b$r*cos(b$a*diff(sector)+sector[1]),
                      lend=2, ljoin=1, col=col[i], lwd=lwd, ...)
            }
        }
    }
}

#' Locate branches graphically
#'
#' @param x \code{\link{polar.clust}} object.
#' @return A data.table with relevant attributes.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
locate.branches <- function(x){
    plot(x, type="n", labels=FALSE, axes=FALSE)
    lines(x)
    bb <- data.table()
    l <- locator(1)
    while(!is.null(l)){
        b <- with(l, list(x = x, y = y, radius = sqrt(x^2+y^2), angle = atan(x/y) %% pi + pi*(x<0)))
        b$phase <- b$angle/2/pi
        split.dist <- apply(
            sweep(x[,list(x=r*sin(a*2*pi), y=r*cos(a*2*pi))], 2, unlist(l))^2,
            1, function(x) sqrt(sum(x)))
        b$nearest.split <- which.min(split.dist)
        do.call(points, x[b$nearest.split, list(x=r*sin(a*2*pi), y=r*cos(a*2*pi), col="red", cex=2.2)])
        do.call(text, x[b$nearest.split, list(x=r*sin(a*2*pi), y=r*cos(a*2*pi), labels=nrow(bb)+1, col="red", pos=3, cex=2.2)])
        b[c("left", "right")] <- x[b$nearest.split, list(i, j)]
        h <- x$r[b$nearest.split]
        b$split.h <- mean(c(x[r > h, min(r)], x[r <= h, max(r)]))
        b$split.k <- sum((x$r < h) != (x$ri < h)) + sum((x$r < h) != (x$rj < h))

        bb <- rbind(bb, b)
        l <- locator(1)
    }
    bb
}


#' Annotate a polar dendrogram
#' 
#' @param x \code{\link{polar.clust}} object.
#' @return A modified \code{\link{polar.clust}} object.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
label.branches <- function(x){
    cat("Click a split in the plot to assign labels to its branches. ")
    cat("Assign a blank label to ignore a branch. Right click to stop labeling.\n")

    plot(x, axes=FALSE, sector=c(0, 2*pi))
    title("Interactive labeling in progess", col.main="grey")
    mtext("Right click or ctrl-c to exit", 3, 0, col="grey")
    #lines(x)
    #points(x)
    if(is.null(x$labi)) x$labi <- as.character(NA)
    if(is.null(x$labj)) x$labj <- as.character(NA)

    l <- locator(1)
    count <- 1
    while(!is.null(l)){
        split.dist <- apply(
            sweep(x[,list(x=r*sin(a*2*pi), y=r*cos(a*2*pi))], 2, unlist(l))^2,
            1, function(x) sqrt(sum(x)))
        b <- which.min(split.dist)

        mark <- x[b, list(
            x = r*sin(a*2*pi),    y = r*cos(a*2*pi),
            xi = ri*sin(ai*2*pi), yi = ri*cos(ai*2*pi),
            xj = rj*sin(aj*2*pi), yj = rj*cos(aj*2*pi)
        )]
        points(mark$x, mark$y, col="springgreen3", cex=2.2)
        text(mark$x, mark$y, labels=count, col="springgreen3", pos=3, cex=2.2)
        points(c(mark$xi, mark$xj), c(mark$yi, mark$yj), col=c("red", "blue"), cex=1.5)

        lab <- readline(sprintf("Left label (#%i red): ", count))
        if(lab != ""){
            x$labi[b] <- lab
            text(mark$xi, mark$yi, labels=lab, col="red", pos=3)
        }
        lab <- readline(sprintf("Right label (#%i blue): ", count))
        if(lab != ""){
            x$labj[b] <- lab
            text(mark$xj, mark$yj, labels=lab, col="blue", pos=3)
        }

        cat("Click a new split to continue labeling or right click to exit.\n")
        l <- locator(1)
        count <- count + 1
    }
    dev.off()
    x
}

