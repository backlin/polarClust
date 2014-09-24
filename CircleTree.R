require(data.table)

circletree <- function(cl, labels=sprintf("#%i", 1:length(cl$order)), hang=.1*max(cl$height)){
    hs <- function(i){
        if(i < 0) return(labels[-i])
        if(tail(cl$height, 1) > 0) cl$height <- tail(cl$height, 1) - cl$height
        list(r = ifelse(cl$merge[i,] > 0,
                        cl$height[abs(cl$merge[i,])],
                        if(hang > 0) cl$height[i] + hang else max(cl$height)),
             node = lapply(cl$merge[i,], hs))
    }
    child.r <- function(st){
        if(is.list(st$node[[1]])) st$node[[1]] <- child.r(st$node[[1]])
        if(is.list(st$node[[2]])) st$node[[2]] <- child.r(st$node[[2]])

        st$node.r <- st$r + c(
            if(is.list(st$node[[1]])) sum(st$node[[1]]$node.r) else 0,
            if(is.list(st$node[[2]])) sum(st$node[[2]]$node.r) else 0
        )
        st
    }
    structure(list(hclust = cl, tree = child.r(hs(nrow(cl$merge)))), class="circletree")
}

lines.circletree <- function(x, lwd=c(1,10), padding=.01, labels=TRUE, label.size=1){
    seg.coord <- function(st, x0=0, y0=0, r=0, range=c(0, 2*pi), weight.range){
        if(missing(weight.range)){
            min.f <- function(st)
                if(is.list(st)) min(c(st$node.r, sapply(st$node, min.f))) else Inf
            weight.range <- c(min.f(st), max(st$node.r))
        }
        aa <- quantile(range,
            quantile(c(0, st$node.r[1]/sum(st$node.r), 1),
                c(padding, .25, .5-padding, .5, .5+padding, .75, 1-padding)
            )
        )
        x1 <- st$r*sin(aa[c(2,6)])
        y1 <- st$r*cos(aa[c(2,6)])
        if(is.character(st$node[[1]]))
            text(x1[1], y1[1], st$node[[1]], adj=c(x1[1] < 0, .5), srt=180/pi*atan(y1[1]/x1[1]), xpd=TRUE, cex=label.size)
        if(is.character(st$node[[2]]))
            text(x1[2], y1[2], st$node[[2]], adj=c(x1[2] < 0, .5), srt=180/pi*atan(y1[2]/x1[2]), xpd=TRUE, cex=label.size)
        rbind(
            data.table(x0 = x0, y0 = y0, x1 = x1, y1 = y1,
                       lwd=lwd[1] + diff(lwd)*(st$node.r-weight.range[1])/diff(weight.range)),
            if(is.list(st$node[[1]])){
                seg.coord(st$node[[1]], x0=x1[1], y0=y1[1], st$r[1], aa[c(1,3)], weight.range)
            } else NULL,
            if(is.list(st$node[[2]])){
                seg.coord(st$node[[2]], x0=x1[2], y0=y1[2], st$r[2], aa[c(5,7)], weight.range)
            } else NULL
        )
    }
    do.call(segments, seg.coord(x$tree))
}

plot.circletree <- function(x, lwd, padding, labels=TRUE, label.size=1, ...){
    max.f <- function(st)
        if(is.character(st)) return(-Inf) else max(st$r, sapply(st$node, max.f))
    m <- max.f(x$tree)
    plot(0, 0, type="n", xlab="", ylab="", xlim=c(-m,m), ylim=c(-m,m), ...)
    lines(x, lwd, padding, labels, label.size)
}

print.circletree <- function(x, ...){
    cat("Circle tree representation of hclust object.\n")
    print(x$hclust)
}



cl <- hclust(dist(iris[-5] + 1*matrix(rnorm(nrow(iris)*(ncol(iris)-1)), nrow(iris))))
x <- circletree(cl, labels=paste("Plant", 1:150))

#png("CircleTree0.png", 700, 400, bg="transparent")
par(mfrow=1:2)
plot(cl)
plot(x, lwd=c(1,6), padding=.01, main="Circle Tree", axes=FALSE, bty="n", label.size=.5)
points(0, 0, cex=1.3, pch=19)
#dev.off()

