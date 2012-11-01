ar.align <-
function(Z, lags, Z.mask=NULL) {
	tau <- nrow(Z)
	n <- ncol(Z)
	low.ndx <- min( min(lags), 0 )
	top.ndx <- max( 0, max(lags) )
	tau0 <- tau - (top.ndx - low.ndx)
	y.sub.ndx <- (1-low.ndx):(tau0-low.ndx)
	Y <- Z[ y.sub.ndx,  , drop=FALSE]
	if( is.null(Z.mask) ) { Y.na <- NULL } else { Y.na <- Z.mask[ y.sub.ndx,  , drop=FALSE ] }
	X <- matrix(0, tau0, n*length(lags))
	for(i in 1:length(lags)) {
		this.lag <- lags[i]
		x.start <- this.lag - low.ndx + 1
		col.rng <- ((i-1)*n+1):(i*n)
		X[ , col.rng ] <- Z[ x.start:(tau0+x.start-1),  , drop=FALSE ]
	}
	return( list(X.sub=X, Y.sub=Y, Y.na=Y.na, y.sub.ndx=y.sub.ndx) )
}
