Hst.sumup <-
function(Hst.ls, Hs=NULL, Ht=NULL) {
    tau <- length(Hst.ls)
    n <- nrow(Hst.ls[[1]])
    big.sum <- 0
    for(i in 1:tau) {
        if( !is.null(Ht) ) { Ht.mx <- matrix( Ht[i, ], n, ncol(Ht), byrow=TRUE ) } else { Ht.mx <- NULL }
        big.sum <- big.sum + crossprod( cbind( Hs, Ht.mx, Hst.ls[[i]] ) )
    }
    return( big.sum )
}
