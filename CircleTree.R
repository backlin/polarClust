
# TODO: Switch branches so that the smallest is close to the center of the
# previous split, e.g. if you a following the left branch put the smallest
# child branch to the left.

require(data.table)
require(RColorBrewer)

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

    # Branch lengths
    ct$li <- na.fill(x$height - x$height[ct$i], hang)
    ct$lj <- na.fill(x$height - x$height[ct$j], hang)


    # Colors
    if(missing(col)) col <- "black"
    if(is.factor(col)) col <- brewer.pal(9, "Set1")[as.integer(col)]
    col <- col2rgb(col)/255
    n <- length(x$order)
    if(ncol(col) != n)
        col <- col[,rep(1:ncol(col), ceiling(n/ncol(col)))[1:n]]
    rgbi <- col[,ifelse(x$merge[,1] < 0, abs(x$merge[,1]), NA)]
    rgbj <- col[,ifelse(x$merge[,2] < 0, abs(x$merge[,2]), NA)]
    ct$coli <- ct$colj <- as.character(NA)
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
labels.CircleTree <- function(object, cex=par("cex"), ...){
    lab <- rbind(
        object[!is.na(labi), list(
            labels = labi,
            x = 1.02*max(ri)*sin(ai),
            y = 1.02*max(ri)*cos(ai),
            adj = ai > pi,
            srt = 90-(ai %% pi)*180/pi
        )],
        object[!is.na(labj), list(
            labels = labj,
            x = 1.02*max(rj)*sin(aj),
            y = 1.02*max(rj)*cos(aj),
            adj = aj > pi,
            srt = 90-(aj %% pi)*180/pi
        )]
    )
    for(i in 1:nrow(lab))
        with(lab[i], text(x, y, labels, srt=srt, adj=c(adj, .5), cex=cex, ...))
}
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

plot.CircleTree <- function(x, points=TRUE, labels, label.size=1, expand, lwd=c(1,10), ...){
    if(missing(labels)) labels <- "labi" %in% names(x)
    if(missing(expand)) expand <- if(labels) 2 else 1.04

    lim <- c(-1,1) * max(x$r)*expand
    plot(0, 0, type="n", xlab="", ylab="", xlim=lim, ylim=lim, ...)
    lines(x, lwd)
    if(points)
        points(x, pch=20)
    if(labels)
        labels(x, cex=label.size)
}

#print.CircleTree <- function(x, ...){
#    cat("Circle tree representation of hclust object. Print function not yet implemented.\n")
#}



