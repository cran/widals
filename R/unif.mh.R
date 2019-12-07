unif.mh <-
function(n, center, sd) {
    w <- sd * sqrt(3)
    a <- center - w
    b <- center + w
    cat(a, b, "\n" )
    x <- runif(1, a, b)
    return(x)
}
