source("CircleTree.R")

cl <- hclust(dist(iris[-5] + .1*matrix(rnorm(nrow(iris)*(ncol(iris)-1)), nrow(iris))))

# Color according to class
x <- CircleTree(cl,
                labels=sprintf("%s %i", iris$Species, 1:150),
                col=iris$Species,
                padding=.1)
plot(x)

# Color branches
x <- CircleTree(cl, labels=sprintf("%s %i", iris$Species, 1:150),
                col=cutree(cl, 3), padding=.002, hang=.1)
plot(x, labels=FALSE)


# ALL samples
load("demo/all_dist.Rdata")
cl <- hclust(all.dist)
x <- CircleTree(cl)
plot(x)

