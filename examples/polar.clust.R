my.data <- iris[-5] + .5*matrix(rnorm(nrow(iris)*(ncol(iris)-1)), nrow(iris))
cl <- hclust(dist(my.data))
x <- polar.clust(cl)

plot(x)

# Color according to class
labels(x) <- sprintf("%s %i", iris$Species, 1:150)
colors(x) <- iris$Species
plot(x, expand=1.2, label.size=.3, las=2)

# Color according to branch
colors(x) <- c(brewer.pal(6, "Set1"))
labels(x) <- NULL # Clear the current labels
labels(x) <- paste("Branch", 1:6)
plot(x, las=1)

# Manual compilation of the plot
s <- c(.4,2.2)*pi
colors(x) <- NULL
plot(x, type="n", labels=FALSE, axes=FALSE, expand=1.1)
branches(x, k=6, sector=s, fill=TRUE, lwd=10, lightness=.3)
lines(x, sector=s)
points(x, pch=19, sector=s)
legend("topright", paste("Branch", 1:6), fill=brewer.pal(6, "Set1"), bg="white")

